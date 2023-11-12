#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub interval_min: f64,
    pub interval_max: f64,
}

impl Interval {
    pub fn new(min: f64, max: f64) -> Self {
        Interval {
            interval_min: min,
            interval_max: max,
        }
    }

    pub fn new_from_intervals(a: Interval, b: Interval) -> Self {
        Interval {
            interval_min: a.interval_min.min(b.interval_min),
            interval_max: a.interval_max.max(b.interval_max),
        }
    }

    pub fn size(&self) -> f64 {
        self.interval_max - self.interval_min
    }

    pub fn expand(&self, delta: f64) -> Interval {
        let padding = delta / 2.0;
        Interval::new(self.interval_min - padding, self.interval_max + padding)
    }

    pub fn contains(&self, x: f64) -> bool {
        self.interval_min <= x && x <= self.interval_max
    }
}

#[cfg(test)]
mod tests {
    use crate::render::Interval;
    use assert_approx_eq::assert_approx_eq;

    // Test creating a new interval, taking floats as an argument
    #[test_log::test]
    fn test_new() {
        let p: Interval = Interval::new(0.1, 0.2);
        assert_approx_eq!(p.interval_min, 0.1);
        assert_approx_eq!(p.interval_max, 0.2);

        let p: Interval = Interval::new(0.3, -0.4);
        assert_approx_eq!(p.interval_min, 0.3);
        assert_approx_eq!(p.interval_max, -0.4);
    }

    #[test_log::test]
    fn test_new_from_intervals() {
        let a: Interval = Interval::new(0.1, 0.2);
        let b: Interval = Interval::new(0.3, -0.4);

        let int: Interval = Interval::new_from_intervals(a, b);

        assert_approx_eq!(int.interval_min, 0.1);
        assert_approx_eq!(int.interval_max, 0.2);
    }

    #[test_log::test]
    fn test_size() {
        let p: Interval = Interval::new(0.1, 0.2);
        assert_approx_eq!(p.size(), 0.1);
    }

    #[test_log::test]
    fn test_expand() {
        let p: Interval = Interval::new(0.1, 0.2);
        assert_approx_eq!(p.expand(0.1).interval_max, 0.25);
        assert_approx_eq!(p.expand(0.1).interval_min, 0.05);
    }

    #[test_log::test]
    fn test_contains() {
        let p: Interval = Interval::new(1.0, 2.0);
        assert_eq!(p.contains(1.5), true);
        assert_eq!(p.contains(2.5), false);
    }
}
