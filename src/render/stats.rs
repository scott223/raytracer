use super::Color;

pub struct Stats {
    colors: Vec<Color>,
    weights: Vec<f64>
}

impl Stats {
    pub fn new(capacity: usize) -> Self {
        Stats {
            colors: Vec::with_capacity(capacity),
            weights: Vec::with_capacity(capacity),
        }
    }

    pub fn push(&mut self, color: Color, weight: f64) {
        self.colors.push(color);
        self.weights.push(weight);
    }

    pub fn mean(&self) -> Color {
        let colors: Color = self.colors.iter()
            .zip(self.weights.iter())
            .map(|(c,w)| *c * *w)
            .sum();

        let weights: f64 = self.weights.iter().sum();

        colors / weights
    }

    pub fn variance(&self) -> f64 {
        if self.colors.len() == 1 {
            0.0
        } else {
            let mean: f64 = self.mean().illuminance_approx();
            let nominator: f64 = self.weights.iter()
                .zip(self.colors.iter())
                .map(|(w, c)| *w * (c.illuminance_approx() - mean).powf(2.0))
                .sum();

            let sum_of_weights: f64 = self.weights.iter().sum();
            let denominator: f64 = sum_of_weights;
            
            nominator / denominator
        }
    }

    pub fn interval(&self) -> f64 {
        1.96 * (self.variance().sqrt() / (self.colors.len() as f64).sqrt())
    }
}
