use std::usize;

pub struct GrayFlips {
    state: Vec<usize>,
}

impl Iterator for GrayFlips {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let k = self.state[0];
        if k == self.state.len() - 1 {
            return None;
        }
        self.state[k] = self.state[k + 1];
        self.state[k + 1] = k + 1;
        if k != 0 {
            self.state[0] = 0;
        }
        return Some(k);
    }
}

impl GrayFlips {
    pub fn new(n: usize) -> GrayFlips {
        GrayFlips {
            state: (0..=n).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn gray_flips_ok() {
        let n: usize = 8;
        let k = 2_usize.pow(n as u32);
        let result: Vec<usize> = GrayFlips::new(n)
            .scan(0, |state, x| {
                *state ^= 1 << x;
                Some(*state)
            })
            .sorted()
            .collect();
        let repr = result.iter().map(|x| x.to_string()).join(",");
        assert!(itertools::equal(result, 1..k), "{} != 1..{}", repr, k);
    }
}
