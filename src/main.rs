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

fn exercise_accuracy_base_example() {
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
    let mut avg_inaccuracy_vec = Vec::<f64>::new();
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
        let avg_inaccuracy = accuracy.fold(0.0, |acc, x| acc + x.abs()) / (accuracy.len() as f64 + 1.0);

        step_count_vec.push(*step_count);

        avg_inaccuracy_vec.push(avg_inaccuracy);

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
    let mut rel_inaccuracy: Vec<f64> = vec![];
    for i in 1..avg_inaccuracy_vec.len() {
        rel_inaccuracy.push(avg_inaccuracy_vec[i-1] / avg_inaccuracy_vec[i]);
    }
    let digits_after_dot = 16;
    println!("{}", "-".repeat(136));
    println!(
        "| {:>5} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$}   |",
        "steps",
        "avg inaccuracy",
        "max inaccuracy",
        "first half inaccuracy",
        "second half inaccuracy",
        "rel inaccuracy",
        width = digits_after_dot + 6
    );
    println!("{}", "-".repeat(136));
    for i in 0..step_count_vec.len() {
        println!(
            "| {:>5} | {:width$e} | {:width$e} | {:width$e} | {:width$e} | {:width$.width$} |",
            step_count_vec[i],
            avg_inaccuracy_vec[i],
            max_inaccuracy_vec[i],
            inaccuracy_in_first_half_vec[i],
            inaccuracy_in_second_half_vec[i],
            if i > 0 { rel_inaccuracy[i - 1] } else { 0.0 },
            width = digits_after_dot + 6
        );
    }
    println!("{}", "-".repeat(136));
}

fn exercise_accuracy_kr_high_accuracy() {
    const LEFT: f64 = 1.0;
    const RIGHT: f64 = 10.0;
    
    // k(r) = 1, q(r) = 0, f(r) = 0 for a simple cylindrical heat conduction case
    let kfunc = LambdaFunction::from(|r: f64| 2.0 * r);
    let qfunc = LambdaFunction::from(|_r: f64| 2.0);
    let ffunc = LambdaFunction::from(|r: f64| 60.0 * r - 120.0);

    // Boundary conditions: T(0) = 0 (center), T(1) = 1 (outer radius)
    let y1 = 30.0;
    let hi2 = 5.0;
    let y2 = 2100.0;

    // The 'n' value for cylindrical symmetry (n=1 for 2D axisymmetric cylindrical case)
    let n = 1;
    let original_function = |r: f64| 30.0 * r;

    let mut step_count_vec = Vec::<i32>::new();
    let mut avg_inaccuracy_vec = Vec::<f64>::new();
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
        let avg_inaccuracy = accuracy.fold(0.0, |acc, x| acc + x.abs()) / (accuracy.len() as f64 + 1.0);

        step_count_vec.push(*step_count);

        avg_inaccuracy_vec.push(avg_inaccuracy);

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
    let mut rel_inaccuracy: Vec<f64> = vec![];
    for i in 1..avg_inaccuracy_vec.len() {
        rel_inaccuracy.push(avg_inaccuracy_vec[i-1] / avg_inaccuracy_vec[i]);
    }
    let digits_after_dot = 16;
    println!("{}", "-".repeat(136));
    println!(
        "| {:>5} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$}   |",
        "steps",
        "avg inaccuracy",
        "max inaccuracy",
        "first half inaccuracy",
        "second half inaccuracy",
        "rel inaccuracy",
        width = digits_after_dot + 6
    );
    println!("{}", "-".repeat(136));
    for i in 0..step_count_vec.len() {
        println!(
            "| {:>5} | {:width$e} | {:width$e} | {:width$e} | {:width$e} | {:width$.width$} |",
            step_count_vec[i],
            avg_inaccuracy_vec[i],
            max_inaccuracy_vec[i],
            inaccuracy_in_first_half_vec[i],
            inaccuracy_in_second_half_vec[i],
            if i > 0 { rel_inaccuracy[i - 1] } else { 0.0 },
            width = digits_after_dot + 6
        );
    }
    println!("{}", "-".repeat(136));
}

fn exercise_accuracy_kr_no_accuracy() {
    const LEFT: f64 = 1.0;
    const RIGHT: f64 = 10.0;
    
    // k(r) = 1, q(r) = 0, f(r) = 0 for a simple cylindrical heat conduction case
    let kfunc = LambdaFunction::from(|r: f64| 2.0 * r * r);
    let qfunc = LambdaFunction::from(|r: f64| 2.0 * r * r);
    let ffunc = LambdaFunction::from(|r: f64| 60.0 * r.powf(4.0) - 480.0 * r * r);

    // Boundary conditions: T(0) = 0 (center), T(1) = 1 (outer radius)
    let y1 = 30.0;
    let hi2 = 5.0;
    let y2 = 135000.0;

    // The 'n' value for cylindrical symmetry (n=1 for 2D axisymmetric cylindrical case)
    let n = 1;
    let original_function = |r: f64| 30.0 * r * r;

    let mut step_count_vec = Vec::<i32>::new();
    let mut avg_inaccuracy_vec = Vec::<f64>::new();
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
        let avg_inaccuracy = accuracy.fold(0.0, |acc, x| acc + x.abs()) / (accuracy.len() as f64 + 1.0);

        step_count_vec.push(*step_count);

        avg_inaccuracy_vec.push(avg_inaccuracy);

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
    let mut rel_inaccuracy: Vec<f64> = vec![];
    for i in 1..avg_inaccuracy_vec.len() {
        rel_inaccuracy.push(avg_inaccuracy_vec[i-1] / avg_inaccuracy_vec[i]);
    }
    let digits_after_dot = 16;
    println!("{}", "-".repeat(136));
    println!(
        "| {:>5} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$} | {:width$.width$}   |",
        "steps",
        "avg inaccuracy",
        "max inaccuracy",
        "first half inaccuracy",
        "second half inaccuracy",
        "rel inaccuracy",
        width = digits_after_dot + 6
    );
    println!("{}", "-".repeat(136));
    for i in 0..step_count_vec.len() {
        println!(
            "| {:>5} | {:width$e} | {:width$e} | {:width$e} | {:width$e} | {:width$.width$} |",
            step_count_vec[i],
            avg_inaccuracy_vec[i],
            max_inaccuracy_vec[i],
            inaccuracy_in_first_half_vec[i],
            inaccuracy_in_second_half_vec[i],
            if i > 0 { rel_inaccuracy[i - 1] } else { 0.0 },
            width = digits_after_dot + 6
        );
    }
    println!("{}", "-".repeat(136));
}

fn main() {
    exercise_accuracy_base_example();
    // exercise_accuracy_kr_high_accuracy();
    // exercise_accuracy_kr_no_accuracy();
}

