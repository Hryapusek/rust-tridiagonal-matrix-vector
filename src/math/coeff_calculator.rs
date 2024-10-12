trait CoeffCalculator<Number> {
    fn calc_a(&self, i: usize) -> Number;
    fn calc_b(&self, i: usize) -> Number;
    fn calc_c(&self, i: usize) -> Number;
    fn calc_g(&self, i: usize) -> Number;
}

trait Function<Number> {
    fn calc(&self, x: Number) -> Number;
}

trait KFunction<Number>: Function<Number> {}

trait QFunction<Number>: Function<Number> {}

/// This module will contain calculator that works
/// when we have **1st condition** at the left and **3rd condition** at the right
mod first_third_calculator {
    use std::marker::PhantomData;

    use super::{KFunction, QFunction, Function};
    use crate::math::coeff_calculator::CoeffCalculator;
    use crate::math::stepping::{NumberTrait, Stepping};

    struct FirstThirdCalculator<
        SteppingObject: Stepping<Number>,
        KFunctionType: KFunction<Number>,
        QFunctionType: QFunction<Number>,
        FunctionType: Function<Number>,
        Number: NumberTrait,
    > {
        stepping: SteppingObject,
        kfunc: KFunctionType,
        qfunc: QFunctionType,
        ffunc: FunctionType,
        y1: Number,
        hi2: Number,
        y2: Number,
        phatnom: PhantomData<Number>,
    }

    impl<
            SteppingObject: Stepping<Number>,
            KFunctionType: KFunction<Number>,
            QFunctionType: QFunction<Number>,
            FunctionType: Function<Number>,
            Number: NumberTrait,
        > FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number>
    {
        fn new(
            stepping: SteppingObject,
            kfunc: KFunctionType,
            qfunc: QFunctionType,
            ffunc: FunctionType,
            y1: Number,
            hi2: Number,
            y2: Number,
        ) -> FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number> {
            FirstThirdCalculator {
                stepping,
                kfunc,
                qfunc,
                ffunc,
                y1,
                hi2,
                y2,
                phatnom: PhantomData,
            }
        }
    }

    impl<
            SteppingObject: Stepping<Number>,
            KFunctionType: KFunction<Number>,
            QFunctionType: QFunction<Number>,
            FunctionType: Function<Number>,
            Number: NumberTrait,
        > CoeffCalculator<Number>
        for FirstThirdCalculator<SteppingObject, KFunctionType, QFunctionType, FunctionType, Number>
    {
        fn calc_a(&self, i: usize) -> Number {
            if i == 0 {
                panic!("i == 0")
            } else if i < self.stepping.points().len() {
                Number::from(-1) * self.kfunc.calc(self.stepping.middle_point(i - 1))
                    / self.stepping.step(i)
            } else {
                panic!("i > self.stepping.steps_count()")
            }
        }

        fn calc_b(&self, i: usize) -> Number {
            if i == 0 {
                Number::from(0)
            } else if i < self.stepping.points().len() - 1 {
                Number::from(-1) * self.kfunc.calc(self.stepping.middle_point(i))
                    / self.stepping.step(i + 1)
            } else {
                panic!("i >= self.stepping.points().len()")
            }
        }

        fn calc_c(&self, i: usize) -> Number {
            if i == 0 {
                Number::from(0)
            } else if i < self.stepping.points().len() - 1 {
                self.kfunc.calc(self.stepping.middle_point(i - 1)) / self.stepping.step(i)
                    + self.kfunc.calc(self.stepping.middle_point(i)) / self.stepping.step(i + 1)
                    + self.qfunc.calc(self.stepping.point(i)) * self.stepping.cross_step(i)
            } else if i == self.stepping.points().len() - 1 {
                self.kfunc.calc(self.stepping.middle_point(i - 1)) / self.stepping.step(i)
                    + self.qfunc.calc(self.stepping.point(i)) * self.stepping.cross_step(i)
                    + self.hi2
            } else {
                panic!("i > self.stepping.points().len()")
            }
        }

        fn calc_g(&self, i: usize) -> Number {
            if i == 0 {
                self.y1
            } else if i < self.stepping.points().len() - 1 {
                self.stepping.cross_step(i) * self.ffunc.calc(self.stepping.point(i))
            } else if i == self.stepping.points().len() - 1 {
                self.stepping.cross_step(i) * self.ffunc.calc(self.stepping.point(i)) + self.y2
            } else {
                panic!("i >= self.stepping.points().len()")
            }
        }
    }
}
