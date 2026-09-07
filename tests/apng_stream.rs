//! **The streaming APNG writer against the whole-slice one.**
//!
//! `apng::Stream` exists because a thousand 1024² frames is 3 GB held only so the encoder can be
//! handed a slice. A writer that produced *different frames* would be a silent corruption — most
//! viewers show a malformed APNG as a still, which is indistinguishable from the failure the
//! duplicate count exists to catch.
//!
//! **The decoded frames are the assertion, and a byte comparison alone would not have caught the
//! bug that was here.** The first cut of `Stream` used `into_stream_writer` — for a *single*
//! image — and concatenated every frame into one: 1390 bytes against 3883, plausible `fcTL` and
//! `fdAT` chunks, and **zero decodable frames**. A smaller file reads as a compression win. So
//! both are asserted: the bytes match, *and* the file decodes to the frames that were pushed.

use prin_rs::output::apng;

fn frames(n: usize, w: usize, h: usize, dupe_at: Option<usize>) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(n);
    for i in 0..n {
        if Some(i) == dupe_at {
            out.push(out[i - 1].clone());
            continue;
        }
        let mut f = vec![0u8; w * h * 3];
        for (k, b) in f.iter_mut().enumerate() {
            *b = ((k * 7 + i * 53) % 251) as u8;
        }
        out.push(f);
    }
    out
}

/// Decode every frame of an APNG back to RGB8.
fn decode(path: &std::path::Path) -> Vec<Vec<u8>> {
    let dec = png::Decoder::new(std::fs::File::open(path).unwrap());
    let mut rdr = dec.read_info().unwrap();
    let mut out = Vec::new();
    let mut buf = vec![0u8; rdr.output_buffer_size()];
    while let Ok(info) = rdr.next_frame(&mut buf) {
        out.push(buf[..info.buffer_size()].to_vec());
    }
    out
}

#[test]
fn the_streamed_animation_decodes_to_the_same_frames() {
    let (w, h) = (17usize, 11usize);
    let fs = frames(6, w, h, None);
    let dir = std::env::temp_dir().join("prin_apng_stream_test");
    let _ = std::fs::create_dir_all(&dir);
    let a = dir.join("whole.png");
    let b = dir.join("stream.png");

    apng::write(a.to_str().unwrap(), w, h, &fs, 1, 30).unwrap();
    let mut s = apng::Stream::begin(b.to_str().unwrap(), w, h, fs.len(), 1, 30).unwrap();
    for f in &fs {
        s.push(f).unwrap();
    }
    s.finish().unwrap();

    let (da, db) = (decode(&a), decode(&b));
    assert_eq!(da.len(), fs.len(), "the whole-slice file decoded {} frames", da.len());
    assert_eq!(db.len(), fs.len(), "the streamed file decoded {} frames", db.len());
    assert_eq!(da, db, "the streamed animation decodes to different frames");
    // And both decode back to what was written, or "they agree" could be two identical failures.
    assert_eq!(db, fs, "the streamed frames are not the frames that were pushed");

    // The control: the comparison must be able to fail. A different frame set has to decode
    // differently, or "identical" is a statement about the decoder ignoring its input.
    let other = frames(6, w, h, Some(3));
    let c = dir.join("other.png");
    apng::write(c.to_str().unwrap(), w, h, &other, 1, 30).unwrap();
    assert_ne!(da, decode(&c), "two different animations decode identically");

    // The two paths differ only in who holds the frames, so the files are byte-identical too.
    // Asserted BESIDE the decode, never instead of it: the bug this test found produced a file
    // that was smaller and structurally plausible, which a size or a byte check alone would have
    // read as an improvement.
    assert_eq!(
        std::fs::read(&a).unwrap(),
        std::fs::read(&b).unwrap(),
        "the streamed file differs byte for byte from the whole-slice one"
    );
}

#[test]
fn the_stream_counts_adjacent_duplicates_the_same_way() {
    let (w, h) = (8usize, 8usize);
    let fs = frames(5, w, h, Some(2));
    let dir = std::env::temp_dir().join("prin_apng_stream_test");
    let _ = std::fs::create_dir_all(&dir);
    let p = dir.join("dup.png");

    let mut s = apng::Stream::begin(p.to_str().unwrap(), w, h, fs.len(), 1, 30).unwrap();
    for f in &fs {
        s.push(f).unwrap();
    }
    let streamed = s.finish().unwrap();
    assert_eq!(streamed, apng::adjacent_duplicates(&fs), "the two duplicate counts disagree");
    // And it is not trivially zero on both sides: this fixture holds exactly one repeat.
    assert_eq!(streamed, 1, "the fixture no longer contains a duplicate, so the test is vacuous");
}
