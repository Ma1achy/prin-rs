//! Append-only checkpoints, so a long run survives being killed.
//!
//! # Why this exists
//!
//! The convergence arm of `wedge_origin` is an ~80 minute indivisible block at 256^2, and it was
//! killed three times, losing everything each time. That is an experiment-design fault, not bad
//! luck: work that cannot be resumed must be sized to fit inside the shortest interruption you
//! expect, and none of these were.
//!
//! # The one thing that would make this dangerous
//!
//! **A checkpoint resumed under different settings is a stale artefact read as current** — the
//! failure this project has hit repeatedly, most expensively with a pre-`dtau`-fix render sitting
//! in a post-fix tree for a day before anyone noticed. So every checkpoint carries a **key**,
//! and [`Ckpt::open`] **refuses to resume** when the key differs rather than mixing the two. The
//! key should be the harness's full config provenance plus whatever else it varies; a key that
//! omits a swept parameter is the `criterion_sweep` filename bug at a new site.
//!
//! # Format
//!
//! `MAGIC | version | key_len | key | (id: u64, len: u32, bytes)*`
//!
//! Append-only and flushed per record, so a kill loses at most the record in flight. A truncated
//! tail — a partial write interrupted mid-record — is **dropped on read**, never returned as a
//! short record: a resumed run must not see half a block as a whole one.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

const MAGIC: &[u8; 8] = b"PRINCKP1";
const VERSION: u32 = 1;

#[derive(Debug)]
pub struct Ckpt {
    w: BufWriter<File>,
}

impl Ckpt {
    /// Open a checkpoint, returning it and everything already recorded.
    ///
    /// Creates the file when absent. **Errors when the stored key differs from `key`** — the
    /// caller should treat that as "these are not the same experiment" and either delete the file
    /// deliberately or use a different path. Silently starting over would be almost as bad as
    /// silently resuming: the run would look fresh and the stale file would still be on disk.
    pub fn open(path: &str, key: &str) -> std::io::Result<(Self, HashMap<u64, Vec<u8>>)> {
        let mut done: HashMap<u64, Vec<u8>> = HashMap::new();

        if Path::new(path).exists() {
            let mut f = File::open(path)?;
            let mut buf = Vec::new();
            f.read_to_end(&mut buf)?;
            let mut o = 0usize;
            let need = |o: usize, n: usize, len: usize| o + n <= len;

            if !need(o, 8 + 4 + 4, buf.len()) || &buf[0..8] != MAGIC {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("{path}: not a checkpoint file"),
                ));
            }
            o = 8;
            let ver = u32::from_le_bytes(buf[o..o + 4].try_into().unwrap());
            o += 4;
            if ver != VERSION {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("{path}: checkpoint version {ver}, expected {VERSION}"),
                ));
            }
            let klen = u32::from_le_bytes(buf[o..o + 4].try_into().unwrap()) as usize;
            o += 4;
            if !need(o, klen, buf.len()) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("{path}: truncated key"),
                ));
            }
            let stored = String::from_utf8_lossy(&buf[o..o + klen]).to_string();
            o += klen;
            if stored != key {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "{path}: checkpoint is for a DIFFERENT configuration and will not be \
                         resumed.\n  stored: {stored}\n  wanted: {key}\n  \
                         Delete the file deliberately, or use another path."
                    ),
                ));
            }
            // Records. A partial tail is dropped rather than returned short.
            while need(o, 12, buf.len()) {
                let id = u64::from_le_bytes(buf[o..o + 8].try_into().unwrap());
                let len = u32::from_le_bytes(buf[o + 8..o + 12].try_into().unwrap()) as usize;
                if !need(o + 12, len, buf.len()) {
                    break;
                }
                done.insert(id, buf[o + 12..o + 12 + len].to_vec());
                o += 12 + len;
            }
            let mut f = OpenOptions::new().append(true).open(path)?;
            f.seek(SeekFrom::End(0))?;
            return Ok((Self { w: BufWriter::new(f) }, done));
        }

        if let Some(p) = Path::new(path).parent() {
            let _ = std::fs::create_dir_all(p);
        }
        let mut f = File::create(path)?;
        f.write_all(MAGIC)?;
        f.write_all(&VERSION.to_le_bytes())?;
        f.write_all(&(key.len() as u32).to_le_bytes())?;
        f.write_all(key.as_bytes())?;
        f.flush()?;
        Ok((Self { w: BufWriter::new(f) }, done))
    }

    /// Record one completed unit of work. Flushed immediately: an unflushed checkpoint is not a
    /// checkpoint, and the cost is one `write` per block rather than per pixel.
    pub fn put(&mut self, id: u64, bytes: &[u8]) -> std::io::Result<()> {
        self.w.write_all(&id.to_le_bytes())?;
        self.w.write_all(&(bytes.len() as u32).to_le_bytes())?;
        self.w.write_all(bytes)?;
        self.w.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> String {
        format!("{}/prin_ckpt_{name}.bin", std::env::temp_dir().display())
    }

    #[test]
    fn a_checkpoint_round_trips_and_resumes() {
        let p = tmp("rt");
        let _ = std::fs::remove_file(&p);
        {
            let (mut c, done) = Ckpt::open(&p, "key-A").unwrap();
            assert!(done.is_empty());
            c.put(7, &[1, 2, 3]).unwrap();
            c.put(9, &[4, 5]).unwrap();
        }
        let (_, done) = Ckpt::open(&p, "key-A").unwrap();
        assert_eq!(done.len(), 2);
        assert_eq!(done[&7], vec![1, 2, 3]);
        assert_eq!(done[&9], vec![4, 5]);
        let _ = std::fs::remove_file(&p);
    }

    /// **The guard that matters.** Resuming under a different configuration would be a stale
    /// artefact read as current -- this project's most expensive recurring failure.
    #[test]
    fn a_different_key_refuses_to_resume() {
        let p = tmp("key");
        let _ = std::fs::remove_file(&p);
        {
            let (mut c, _) = Ckpt::open(&p, "key-A").unwrap();
            c.put(1, &[1]).unwrap();
        }
        let e = Ckpt::open(&p, "key-B").unwrap_err();
        assert!(format!("{e}").contains("DIFFERENT configuration"));
        let _ = std::fs::remove_file(&p);
    }

    /// A kill mid-write leaves a partial record. It must be **dropped**, never returned short:
    /// a resumed run seeing half a block as a whole one is worse than losing the block.
    #[test]
    fn a_truncated_tail_is_dropped_not_returned_short() {
        let p = tmp("trunc");
        let _ = std::fs::remove_file(&p);
        {
            let (mut c, _) = Ckpt::open(&p, "k").unwrap();
            c.put(1, &[9, 9, 9, 9]).unwrap();
            c.put(2, &[8, 8, 8, 8]).unwrap();
        }
        // Chop the last record in half, as an interrupted write would.
        let mut b = std::fs::read(&p).unwrap();
        let n = b.len();
        b.truncate(n - 2);
        std::fs::write(&p, &b).unwrap();

        let (_, done) = Ckpt::open(&p, "k").unwrap();
        assert_eq!(done.len(), 1, "the partial record must not survive");
        assert_eq!(done[&1], vec![9, 9, 9, 9]);
        assert!(!done.contains_key(&2));
        let _ = std::fs::remove_file(&p);
    }
}
