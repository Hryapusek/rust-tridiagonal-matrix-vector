use std::ops::{Add, Div, Sub};

trait Stepping<Number>
where
    Number: Div<i32, Output = Number> + std::ops::Sub<Output = Number>,
{
    fn steps_count(&self) -> usize;

    /// **NOTE:** `i > 0` and `i <= self.steps_count()`
    fn step(&self, i: usize) -> Number;

    /// **NOTE:** `i >= 0` and `i <= self.steps_count()`
    fn cross_step(&self, i: usize) -> Number {
        assert!(i <= self.steps_count());
        if i == 0 {
            self.step(1) / 2
        } else if i == self.steps_count() {
            self.step(i) / 2
        } else {
            (self.step(i + 1) - self.step(i)) / 2
        }
    }
}

struct IntervalSplitter<Number> {
    points: Vec<Number>,
    /// points.size() - 1
    steps: Vec<Number>,
    /// points.size()
    cross_steps: Vec<Number>,
}

impl<Number> IntervalSplitter<Number>
where
    Number: Div<i32, Output = Number> + Sub<Output = Number>,
    Number: Copy,
    Number: Add<Output = Number>,
{
    fn new(points: Vec<Number>) -> IntervalSplitter<Number> {
        let mut steps: Vec<Number> = vec![];
        for i in 1..points.len() {
            steps.push(points[i] - points[i - 1]);
        }

        let mut cross_steps: Vec<Number> = vec![];
        for i in 0..points.len() {
            if i == 0 {
                cross_steps.push(steps[i] / 2);
            } else if i == points.len() - 1 {
                cross_steps.push(steps[i - 1] / 2);
            } else {
                cross_steps.push((steps[i - 1] + steps[i]) / 2);
            }
        }

        IntervalSplitter {
            points,
            steps,
            cross_steps,
        }
    }

    fn points(&self) -> &Vec<Number> {
        &self.points
    }

    fn point(&self, i: usize) -> Number {
        self.points[i]
    }
}

impl<Number> Stepping<Number> for IntervalSplitter<Number>
where
    Number: Div<i32, Output = Number> + std::ops::Sub<Output = Number> + Copy,
{
    fn steps_count(&self) -> usize {
        self.points.len() - 1
    }

    fn step(&self, i: usize) -> Number {
        self.steps[i - 1]
    }

    fn cross_step(&self, i: usize) -> Number {
        self.cross_steps[i - 1]
    }
}
