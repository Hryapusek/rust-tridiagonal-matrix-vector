mod math;

use math::{coeff_calculator::*, solver};
use math::stepping::IntervalSplitter;

fn main() {
    // Set radial points (ranging from 0 to 1, simulating the radius of the cylinder)
    let points = (1..=11).map(|i| (i as f64)).collect::<Vec<f64>>();

    let splitter = IntervalSplitter::new(points);

    // k(r) = 1, q(r) = 0, f(r) = 0 for a simple cylindrical heat conduction case
    let kfunc = LambdaFunction::from(|_r| 2.0);
    let qfunc = LambdaFunction::from(|_r| 3.0);
    let ffunc = LambdaFunction::from(|r| -4.0 / r + 54.0 + 6.0 * r);

    // Boundary conditions: T(0) = 0 (center), T(1) = 1 (outer radius)
    let y1 = 20.0;
    let hi2 = 5.0;
    let y2 = 204.0;

    // The 'n' value for cylindrical symmetry (n=1 for 2D axisymmetric cylindrical case)
    let n = 1;

    let coeff_calculator_v =
        math::coeff_calculator::first_third_calculator::FirstThirdCalculator::new(
            splitter,
            kfunc,
            qfunc,
            ffunc,
            y1,
            hi2,
            y2,
            n,
        );

    // Build the matrix and right-hand side vector
    let (A, g) = math::coeff_calculator::matrix_building::build_tridiagonal_matrix(&coeff_calculator_v);

    println!("Matrix A:\n{}", A);
    println!("Vector g:\n{}", g);

    // Solve the system to get the temperature distribution
    let v = solver::solve(&A, &g);

    println!("Solution vector v (temperature distribution):\n{}", v);
}

