pub trait CoeffCalculator<Number> {
    fn calc_a(&self, i: usize) -> Number;
    fn calc_b(&self, i: usize) -> Number;
    fn calc_c(&self, i: usize) -> Number;
    fn calc_g(&self, i: usize) -> Number;
    fn size(&self) -> usize;
}

pub struct LambdaFunction<Number, Func: Fn(Number) -> Number> {
    ffunc: Func,
    phatnom: std::marker::PhantomData<Number>,
}

impl<Number, Func: Fn(Number) -> Number> Function<Number> for LambdaFunction<Number, Func> {
    fn calc(&self, x: Number) -> Number {
        (self.ffunc)(x)
    }
}

impl<Number, Func: Fn(Number) -> Number> std::convert::From<Func> for LambdaFunction<Number, Func> {
    fn from(ffunc: Func) -> Self {
        LambdaFunction {
            ffunc,
            phatnom: std::marker::PhantomData,
        }
    }
}

pub trait Function<Number> {
    fn calc(&self, x: Number) -> Number;
}

/// This module will contain calculator that works
/// when we have **1st condition** at the left and **3rd condition** at the right
pub mod first_third_calculator {
    use std::marker::PhantomData;

    use super::Function;
    use crate::math::coeff_calculator::CoeffCalculator;
    use crate::math::stepping::{NumberTrait, Stepping};

    pub struct FirstThirdCalculator<
        SteppingObject: Stepping<Number>,
        KFunctionType: Function<Number>,
        QFunctionType: Function<Number>,
        FunctionType: Function<Number>,
        Number: NumberTrait,
    > {
        stepping: SteppingObject,
        kfunc: KFunctionType,
        qfunc: QFunctionType,
        ffunc: FunctionType,
        y1: Number,
        y2: Number,
        n: u16,
    }

    impl<
            SteppingObject: Stepping<Number>,
            KFunctionType: Function<Number>,
            QFunctionType: Function<Number>,
            FunctionType: Function<Number>,
            Number: NumberTrait,
        > FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number>
    {
        pub fn new(
            stepping: SteppingObject,
            kfunc: KFunctionType,
            qfunc: QFunctionType,
            ffunc: FunctionType,
            y1: Number,
            y2: Number,
            n: u16,
        ) -> FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number>
        {
            FirstThirdCalculator {
                stepping,
                kfunc,
                qfunc,
                ffunc,
                y1,
                y2,
                n,
            }
        }
    }

    impl<
            SteppingObject: Stepping<Number>,
            KFunctionType: Function<Number>,
            QFunctionType: Function<Number>,
            FunctionType: Function<Number>,
            Number: NumberTrait,
        > CoeffCalculator<Number>
        for FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number>
    {
        fn calc_a(&self, i: usize) -> Number {
            if i == 0 {
                panic!("i == 0")
            } else if i < self.stepping.points().len() {
                Number::from(-1)
                    * self.stepping.middle_point(i - 1).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i - 1))
                    / self.stepping.step(i)
            } else {
                panic!("i > self.stepping.steps_count()")
            }
        }

        fn calc_b(&self, i: usize) -> Number {
            if i == 0 {
                self.stepping.middle_point(i).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i))
                    / self.stepping.step(i + 1)
            } else if i < self.stepping.points().len() - 1 {
                self.stepping.middle_point(i).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i))
                    / self.stepping.step(i + 1)
            } else {
                panic!("i >= self.stepping.points().len()")
            }
        }

        fn calc_c(&self, i: usize) -> Number {
            if i == 0 {
                self.stepping.middle_point(i).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i))
                    / self.stepping.step(i + 1)
                    + Number::from(1) / Number::from(self.n + 1)
                        * self.stepping.cross_step(i)
                        * self.stepping.middle_point(i).pow(self.n)
                        * self.qfunc.calc(self.stepping.point(i))
            } else if i < self.stepping.points().len() - 1 {
                self.stepping.middle_point(i - 1).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i - 1))
                    / self.stepping.step(i)
                    + self.stepping.middle_point(i).pow(self.n)
                        * self.kfunc.calc(self.stepping.middle_point(i))
                        / self.stepping.step(i + 1)
                    + self.stepping.cross_step(i)
                        * self.stepping.middle_point(i).pow(self.n)
                        * self.qfunc.calc(self.stepping.point(i))
            } else if i == self.stepping.points().len() - 1 {
                self.stepping.middle_point(i - 1).pow(self.n)
                    * self.kfunc.calc(self.stepping.middle_point(i - 1))
                    / self.stepping.step(i)
                    + self.stepping.cross_step(i)
                        * self.stepping.middle_point(i-1).pow(self.n)
                        * self.qfunc.calc(self.stepping.point(i))
            } else {
                panic!("i > self.stepping.points().len()")
            }
        }

        fn calc_g(&self, i: usize) -> Number {
            if i == 0 {
                self.y1
            } else if i < self.stepping.points().len() - 1 {
                self.stepping.cross_step(i)
                    * self.stepping.point(i).pow(self.n)
                    * self.ffunc.calc(self.stepping.point(i))
            } else if i == self.stepping.points().len() - 1 {
                self.stepping.cross_step(i)
                    * self.stepping.point(i).pow(self.n)
                    * self.ffunc.calc(self.stepping.point(i))
                    + self.stepping.point(i).pow(self.n) * self.y2
            } else {
                panic!("i >= self.stepping.points().len()")
            }
        }

        fn size(&self) -> usize {
            self.stepping.points().len()
        }
    }
}

pub mod matrix_building {
    use crate::CoeffCalculator;
    use nalgebra::{DMatrix, DVector};

    /// Builds the tridiagonal matrix and the right-hand side vector
    /// using the provided coefficient calculator.
    pub fn build_tridiagonal_matrix<Number, Calculator>(
        calculator: &Calculator,
    ) -> (DMatrix<Number>, DVector<Number>)
    where
        Number: nalgebra::RealField, // Ensures we can work with nalgebra numbers
        Calculator: CoeffCalculator<Number>,
    {
        let n = calculator.size();
        // Create an NxN matrix initialized to zeros
        let mut matrix = DMatrix::zeros(n, n);
        // Create a vector for the right-hand side of the equation
        let mut rhs_vector = DVector::zeros(n);

        // Fill the matrix with values from the coefficient calculator
        for i in 0..n {
            if i > 0 {
                matrix[(i, i - 1)] = calculator.calc_a(i); // Sub-diagonal (below the main diagonal)
            }
            matrix[(i, i)] = calculator.calc_c(i); // Main diagonal
            if i < n - 1 {
                matrix[(i, i + 1)] = calculator.calc_b(i); // Super-diagonal (above the main diagonal)
            }

            // Fill the right-hand side vector (the g vector)
            rhs_vector[i] = calculator.calc_g(i);
        }

        (matrix, rhs_vector)
    }
}
