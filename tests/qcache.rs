//! Tests for the `PRQC` reader.
//!
//! The module doc has said *"a reader never guesses"* since the writer was written, and for a
//! fortnight there was no reader: the committed `.qcache` files could be produced and not
//! consumed, so `total_substeps` -- the one machine-independent cost column in this project --
//! sat on disk unreadable, and the cost ledger could not be built. These tests are about the two
//! ways a table reader silently returns the wrong number.

use prin_rs::output::qcache;

fn read(stem: &str) -> Option<qcache::QuadRows> {
    let path = format!("results/payload/{stem}_t13_L6.qcache");
    let f = std::fs::File::open(path).ok()?;
    Some(qcache::read(&mut std::io::BufReader::new(f)).expect("read qcache"))
}

/// **A space-containing value on a shared header line truncates at the first space.**
///
/// `fcache` documents this and reads such fields line-wise; `qcache` writes `region={} chart={}
/// ...` on one line and three region names carry a space, so a whitespace-token read of
/// `deep interior` returns `deep`. It went unnoticed for the only reason such things do: nothing
/// read the file.
///
/// The negative control is the point. Asserting `region == "deep interior"` alone passes on any
/// parse that happens to work today; asserting that the *naive* parse gives something different
/// is what fires if the careful one is ever simplified back.
#[test]
fn a_region_name_with_a_space_survives_the_header() {
    let Some(qr) = read("deep_interior") else {
        eprintln!("no committed deep_interior qcache; nothing to test");
        return;
    };
    assert_eq!(qr.region, "deep interior", "the space must survive");

    let naive = qr
        .header
        .split_whitespace()
        .find_map(|t| t.strip_prefix("region="))
        .unwrap_or_default();
    assert_eq!(naive, "deep", "the naive parse is expected to truncate -- this is the control");
    assert_ne!(
        naive, qr.region,
        "if these agree the control has gone vacuous and this test no longer discriminates"
    );
}

/// **Columns are addressed by name, never by position.** `FIELDS` has been appended to once
/// already (v2 added the relative mask and the gradient), so a positional read of a v1 file under
/// v2's indices returns a different column with no error anywhere.
#[test]
fn columns_are_found_by_name_and_absence_is_not_zero() {
    let Some(qr) = read("near-field") else {
        eprintln!("no committed near-field qcache; nothing to test");
        return;
    };
    let c = qr.col("total_substeps").expect("v2 carries total_substeps");
    assert_eq!(qr.rows[0][c], qr.get(0, "total_substeps"));
    assert!(qr.get(0, "total_substeps") > 0.0, "the root quad integrated something");

    // An absent column is `NaN`, not `0.0`. A sum over `NaN` is loud; a sum over zero is a
    // plausible wrong answer, and a cost table reading zero would say the tree was free.
    assert!(qr.col("no_such_field").is_none());
    assert!(qr.get(0, "no_such_field").is_nan(), "an absent field must not read as zero");
}

/// The level column indexes a complete tree: `(4^(L+1) - 1)/3` quads, `4^l` at level `l`.
/// Cheap, and it is what says the key reconstruction is not off by a level.
#[test]
fn the_cache_holds_a_complete_tree() {
    let Some(qr) = read("far") else {
        eprintln!("no committed far qcache; nothing to test");
        return;
    };
    let mut n = vec![0usize; qr.levels as usize + 1];
    for i in 0..qr.rows.len() {
        n[qr.key(i).0 as usize] += 1;
    }
    for (l, &c) in n.iter().enumerate() {
        assert_eq!(c, 4usize.pow(l as u32), "level {l} should hold 4^{l} quads");
    }
    assert_eq!(qr.rows.len(), (4usize.pow(qr.levels + 1) - 1) / 3);
}

/// **A truncated file is refused, not read short.** The header names the fields and the body
/// declares its own count; if they ever disagree the file is damaged or from a build whose header
/// and body diverged, and reading the shorter of the two shifts every column by the difference.
#[test]
fn a_damaged_file_is_refused_rather_than_read_short() {
    let Some(bytes) = std::fs::read("results/payload/far_t13_L6.qcache").ok() else {
        eprintln!("no committed far qcache; nothing to test");
        return;
    };
    // Whole file reads.
    assert!(qcache::read(&mut bytes.as_slice()).is_ok());
    // Body cut in half: `read_exact` on the value block must fail rather than return the prefix.
    let cut = bytes.len() / 2;
    assert!(
        qcache::read(&mut bytes[..cut].as_ref()).is_err(),
        "a truncated body must be refused"
    );
    // Wrong magic.
    let mut wrong = bytes.clone();
    wrong[0] = b'X';
    assert!(qcache::read(&mut wrong.as_slice()).is_err(), "a foreign file must be refused");
}
