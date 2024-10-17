mod math;

use math::{coeff_calculator::*, solver};
use math::stepping::{IntervalSplitter, Stepping};
use nalgebra::DVector;

fn generate_points(left: f64, right: f64, step_count: i32) -> Vec<f64> {
    let step_size = (right - left) / step_count as f64;
    let mut points: Vec<f64> = Vec::<f64>::new();

    for i in 0..=step_count {
        points.push(left + i as f64 * step_size);
    }

    points
}

fn exercise_accuracy() {
    const LEFT: f64 = 1.0;
    const RIGHT: f64 = 11.0;

    // k(r) = 1, q(r) = 0, f(r) = 0 for a simple cylindrical heat conduction case
    let kfunc = LambdaFunction::from(|_r: f64| 2.0);
    let qfunc = LambdaFunction::from(|_r: f64| 3.0);
    let ffunc = LambdaFunction::from(|r: f64| -4.0 / r + 54.0 + 6.0 * r);

    // Boundary conditions: T(0) = 0 (center), T(1) = 1 (outer radius)
    let y1 = 20.0;
    let hi2 = 5.0;
    let y2 = 204.0;

    // The 'n' value for cylindrical symmetry (n=1 for 2D axisymmetric cylindrical case)
    let n = 1;
    let original_function = |r: f64| y1 + (r-1.0) * 2.0;

    let mut step_count_vec = Vec::<i32>::new();
    let mut sum_inaccuracy_vec = Vec::<f64>::new();
    let mut max_inaccuracy_vec = Vec::<f64>::new();
    let mut inaccuracy_in_first_half_vec = Vec::<f64>::new();
    let mut inaccuracy_in_second_half_vec = Vec::<f64>::new();
    

    for step_count in [2, 4, 8, 16, 32, 64, 128].iter() {
        let points = generate_points(LEFT, RIGHT, *step_count);
        let splitter = IntervalSplitter::new(points.clone());

        let coeff_calculator_v =
            math::coeff_calculator::first_third_calculator::FirstThirdCalculator::new(
                splitter,
                &kfunc,
                &qfunc,
                &ffunc,
                y1,
                hi2,
                y2,
                n,
            );
        
        let (A, g) = math::coeff_calculator::matrix_building::build_tridiagonal_matrix(&coeff_calculator_v);

        let calculated_v = solver::solve(&A, &g);
        let expected_v: DVector<f64> = DVector::from_vec(points.iter().map(|x| original_function(*x)).collect());
        let accuracy = &calculated_v - &expected_v;
        let sum_inaccuracy = accuracy.fold(0.0, |acc, x| acc + x.abs());

        step_count_vec.push(*step_count);

        sum_inaccuracy_vec.push(sum_inaccuracy);

        max_inaccuracy_vec.push(accuracy.iter().fold(0.0, |max, x| x.abs().max(max)));

        let first_half_inaccuracy: f64 = accuracy
            .iter()
            .take(accuracy.len() / 2)
            .fold(0.0, |acc, x| acc + x.abs());

        let second_half_inaccuracy: f64 = accuracy
            .iter()
            .skip(accuracy.len() / 2)
            .fold(0.0, |acc, x| acc + x.abs());
        inaccuracy_in_first_half_vec.push(first_half_inaccuracy);
        inaccuracy_in_second_half_vec.push(second_half_inaccuracy);
    }
    let digits_after_dot = 16;
    println!("{}", "-".repeat(109));
    println!(
        "| {:>5} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$} |",
        "steps",
        "sum inaccuracy",
        "max inaccuracy",
        "first half inaccuracy",
        "second half inaccuracy",
        width = digits_after_dot + 6
    );
    println!("{}", "-".repeat(109));
    for i in 0..step_count_vec.len() {
        println!(
            "| {:>5} | {:width$e} | {:width$e} | {:width$e} | {:width$e} |",
            step_count_vec[i],
            sum_inaccuracy_vec[i],
            max_inaccuracy_vec[i],
            inaccuracy_in_first_half_vec[i],
            inaccuracy_in_second_half_vec[i],
            width = digits_after_dot + 6
        );
    }
    println!("{}", "-".repeat(109));
    

}

fn main() {
    exercise_accuracy();
}

