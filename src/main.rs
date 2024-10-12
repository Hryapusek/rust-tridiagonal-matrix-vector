mod math;

use math::{coeff_calculator::*, solver};
use math::stepping::IntervalSplitter;

fn main() {
    let points = (0..10).map(|i| (i as f64) / 10.0).collect::<Vec<f64>>();

    let splitter = IntervalSplitter::new(points);

    let kfunc = LambdaFunction::from(|x| x);
    let qfunc = LambdaFunction::from(|x| x);
    let ffunc = LambdaFunction::from(|x| x);
    let y1 = 0.0;
    let hi2 = 0.0;
    let y2 = 0.0;

    let coeff_calculator_v =
        math::coeff_calculator::first_third_calculator::FirstThirdCalculator::new(
            splitter,
            kfunc,
            qfunc,
            ffunc,
            y1,
            hi2,
            y2,
        );
    
    let (A, g) = math::coeff_calculator::matrix_building::build_tridiagonal_matrix(&coeff_calculator_v);

    println!("A: {}", A);
    println!("g: {}", g);

    let v = solver::solve(&A, &g);

    println!("v: {}", v);

}
