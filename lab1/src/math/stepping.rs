pub trait NumberTrait:
    std::ops::Div<Self, Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Add<Output = Self>
    + std::ops::Mul<Output = Self>
    + num_traits::Pow<u16, Output = Self>
    + Sized
    + Copy
    + std::convert::From<i32>
    + std::convert::From<u16>
{
}

impl NumberTrait for i32 {}
impl NumberTrait for f64 {}

pub trait Stepping<Number>
where
    Number: NumberTrait,
{
    /// **NOTE:** `i > 0` and `i <= self.steps_count()`
    fn step(&self, i: usize) -> Number;

    /// **NOTE:** `i >= 0` and `i <= self.steps_count()`
    fn cross_step(&self, i: usize) -> Number {
        assert!(i < self.points().len());
        if i == 0 {
            self.step(1) / Number::from(2)
        } else if i == self.points().len() - 1 {
            self.step(i) / Number::from(2)
        } else {
            (self.step(i + 1) - self.step(i)) / Number::from(2)
        }
    }

    fn points(&self) -> &Vec<Number>;

    fn point(&self, i: usize) -> Number;

    /// Returns the middle point between `i` and `i+1`
    fn middle_point(&self, i: usize) -> Number;
}

pub struct IntervalSplitter<Number> {
    points: Vec<Number>,
    /// points.size() - 1
    steps: Vec<Number>,
    /// points.size()
    cross_steps: Vec<Number>,
}

impl<Number> IntervalSplitter<Number>
where
    Number: NumberTrait,
{
    pub fn new(points: Vec<Number>) -> IntervalSplitter<Number> {
        let mut steps: Vec<Number> = vec![];
        for i in 1..points.len() {
            steps.push(points[i] - points[i - 1]);
        }

        let mut cross_steps: Vec<Number> = vec![];
        for i in 0..points.len() {
            if i == 0 {
                cross_steps.push(steps[i] / Number::from(2));
            } else if i == points.len() - 1 {
                cross_steps.push(steps[i - 1] / Number::from(2));
            } else {
                cross_steps.push((steps[i - 1] + steps[i]) / Number::from(2));
            }
        }

        IntervalSplitter {
            points,
            steps,
            cross_steps,
        }
    }
}

impl<Number> Stepping<Number> for IntervalSplitter<Number>
where
    Number: NumberTrait,
{
    fn step(&self, i: usize) -> Number {
        self.steps[i - 1]
    }

    fn cross_step(&self, i: usize) -> Number {
        self.cross_steps[i]
    }

    fn points(&self) -> &Vec<Number> {
        &self.points
    }

    fn point(&self, i: usize) -> Number {
        self.points[i]
    }

    fn middle_point(&self, i: usize) -> Number {
        assert!(i < self.points.len() - 1);
        (self.points[i] + self.points[i + 1]) / Number::from(2)
    }
}
