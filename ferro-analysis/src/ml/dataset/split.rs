//! Frame-level train / valid / test split of one system.

use ferro_core::error::ChemError;

use super::merge::shuffle_order;

/// Frame-level train/valid/test split.
///
/// Membership is drawn from a shuffled order — taking the tail as a test set
/// would hand it the end of the trajectory, which is one contiguous stretch of
/// a single state. The indices of each part are then sorted back into
/// trajectory order, so only membership is random and the output stays
/// byte-identical across runs with the same seed.
///
/// Note what this cannot fix: frames of one MD run are correlated, so a
/// frame-level test set still shares its neighbourhood with the training set
/// and reads optimistic. A split across whole systems is the honest estimate;
/// this one is the convenient one.
#[derive(Clone, Copy, Debug)]
pub struct Split {
    valid: f64,
    test: f64,
    seed: u64,
}

impl Split {
    /// Parses `--ratio`: `train:valid:test`, or `train:test` with two fields.
    ///
    /// The numbers are weights, not fractions — `8:1:1` and `80:10:10` are the
    /// same split — so no one has to make them sum to one. Two fields mean
    /// train:test because that is the pair NEP asks for (`train.xyz` +
    /// `test.xyz`) and the order the phrase "train/test split" already implies;
    /// every run prints the three parts by name, so a misread shows up on the
    /// first line of output rather than in a silently mislabelled dataset.
    pub fn parse(spec: &str, seed: u64) -> ferro_core::Result<Self> {
        let invalid = |msg: String| ChemError::ValidationError(msg);
        let fields: Vec<&str> = spec.split(':').map(|f| f.trim()).collect();
        let bad = || {
            invalid(format!(
                "--ratio must look like 8:1:1 (train:valid:test) or 9:1 \
                 (train:test), got {spec:?}"
            ))
        };
        if !(2..=3).contains(&fields.len()) {
            return Err(bad());
        }
        let mut w = Vec::with_capacity(3);
        for f in &fields {
            let v: f64 = f.parse().map_err(|_| bad())?;
            if !v.is_finite() || v < 0.0 {
                return Err(invalid(format!("--ratio takes non-negative numbers, got {f:?}")));
            }
            w.push(v);
        }
        let (train, valid, test) = match w.len() {
            2 => (w[0], 0.0, w[1]),
            _ => (w[0], w[1], w[2]),
        };
        let total = train + valid + test;
        if total <= 0.0 {
            return Err(invalid("--ratio is all zeros, so nothing would be written".into()));
        }
        if valid <= 0.0 && test <= 0.0 {
            return Err(invalid(format!("--ratio {spec:?} holds out nothing; omit --ratio instead")));
        }
        Ok(Split { valid: valid / total, test: test / total, seed })
    }

    /// The split off state, used when `--ratio` is absent.
    pub fn off(seed: u64) -> Self {
        Split { valid: 0.0, test: 0.0, seed }
    }

    pub fn is_off(&self) -> bool {
        self.valid <= 0.0 && self.test <= 0.0
    }

    /// `[train, valid, test]` frame indices of a system of `n` frames, each in
    /// ascending order. `who` names the system in the error.
    pub fn parts(&self, n: usize, who: &str) -> ferro_core::Result<[Vec<usize>; 3]> {
        if self.is_off() {
            return Ok([(0..n).collect(), Vec::new(), Vec::new()]);
        }
        // A ratio rounds up to at least one frame: a ratio given and 0 frames
        // taken is silently no validation set
        let take = |r: f64| -> usize {
            if r <= 0.0 { 0 } else { ((n as f64 * r).round() as usize).max(1) }
        };
        let (n_valid, n_test) = (take(self.valid), take(self.test));
        if n_valid + n_test >= n {
            return Err(ChemError::ValidationError(format!(
                "{who}: {n} frame(s) cannot give {n_valid} validation + {n_test} test \
                 and still leave a training set"
            )));
        }
        let order = shuffle_order(n, self.seed);
        let mut valid: Vec<usize> = order[..n_valid].to_vec();
        let mut test: Vec<usize> = order[n_valid..n_valid + n_test].to_vec();
        let mut train: Vec<usize> = order[n_valid + n_test..].to_vec();
        for v in [&mut train, &mut valid, &mut test] {
            v.sort_unstable();
        }
        Ok([train, valid, test])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::merge::DEFAULT_SEED;

    fn split_of(spec: &str) -> Split {
        Split::parse(spec, DEFAULT_SEED).unwrap()
    }

    #[test]
    fn a_split_partitions_every_frame_exactly_once() {
        let [tr, va, te] = split_of("7:2:1").parts(20, "x").unwrap();
        assert_eq!((tr.len(), va.len(), te.len()), (14, 4, 2));
        let mut all: Vec<usize> = tr.iter().chain(&va).chain(&te).copied().collect();
        all.sort_unstable();
        assert_eq!(all, (0..20).collect::<Vec<_>>());
    }

    #[test]
    fn each_part_keeps_trajectory_order() {
        // membership is drawn at random, but each part is in frame order, so the
        // products are byte-identical under one seed
        let [tr, va, te] = split_of("2:1:1").parts(40, "x").unwrap();
        for part in [&tr, &va, &te] {
            assert!(part.windows(2).all(|w| w[0] < w[1]), "{part:?}");
        }
    }

    #[test]
    fn a_split_is_not_the_tail_of_the_trajectory() {
        // cutting off the tail would make the test set the last few frames; the
        // draw must shuffle first
        let [_, _, te] = split_of("3:0:1").parts(40, "x").unwrap();
        assert!(te.iter().any(|&i| i < 30), "test set looks like a tail: {te:?}");
    }

    #[test]
    fn a_ratio_too_small_to_reach_one_frame_still_gets_one() {
        // a ratio given and 0 frames taken is silently no validation set
        let [tr, va, _] = split_of("99:1:0").parts(20, "x").unwrap();
        assert_eq!(va.len(), 1);
        assert_eq!(tr.len(), 19);
    }

    #[test]
    fn a_split_that_leaves_no_training_frames_is_an_error() {
        let e = split_of("0.001:1:1").parts(4, "sysA").unwrap_err();
        let msg = e.to_string();
        assert!(msg.contains("sysA"), "{msg}");
        // a continuation written as `\\` leaves a literal backslash and an indent in the message
        assert!(!msg.contains('\\') && !msg.contains("  "), "the continuation is broken: {msg:?}");
    }

    #[test]
    fn a_ratio_is_weights_not_fractions() {
        // 8:1:1 and 80:10:10 are one split: nobody has to make them sum to 1
        let a = split_of("8:1:1").parts(100, "x").unwrap();
        let b = split_of("80:10:10").parts(100, "x").unwrap();
        assert_eq!(a[0].len(), b[0].len());
        assert_eq!((a[0].len(), a[1].len(), a[2].len()), (80, 10, 10));
    }

    #[test]
    fn two_fields_mean_train_and_test() {
        // the 1 of 9:1 is test, not valid; every run prints the parts by name, so a
        // misreading shows on the first line of output
        let [tr, va, te] = split_of("9:1").parts(100, "x").unwrap();
        assert_eq!((tr.len(), va.len(), te.len()), (90, 0, 10));
    }

    #[test]
    fn a_malformed_ratio_is_rejected_with_the_expected_shape() {
        for bad in ["", "8", "8:1:1:1", "8:x:1", "8:-1:1", "0:0:0", "10:0:0"] {
            let e = Split::parse(bad, DEFAULT_SEED).unwrap_err();
            let msg = e.to_string();
            assert!(
                msg.contains("--ratio"),
                "{bad:?} produced an unhelpful message: {msg}"
            );
        }
        assert!(Split::parse("8:1:1", DEFAULT_SEED).is_ok());
        assert!(Split::parse(" 8 : 1 : 1 ", DEFAULT_SEED).is_ok());
        assert!(Split::parse("0.8:0.1:0.1", DEFAULT_SEED).is_ok());
    }
}
