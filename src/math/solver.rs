use core::panic;

use nalgebra::DMatrix;
use nalgebra::DVector;

/// This type allow us to specify the exact type of the number we are using in the solver
#[allow(unused)]
pub type Number = f64;

/// This function checks if the matrix is tri-diagonally dominant.
///
/// Dominant tridiagonal matrices are the ones that have only
/// non-zero elements on the main diagonal and upper and lower diagonals
#[allow(unused)]
fn check_if_three_diagonals(a: &DMatrix<Number>) -> () {
    if a.nrows() != a.ncols() {
        panic!("The matrix is not square");
    }
    for row in 0..a.nrows() {
        for col in 0..a.ncols() {
            if !(row.abs_diff(col) <= 1) && (a[(row, col)] != 0.0) {
                panic!("Matrix is not diagonally dominant");
            }
        }
    }
}

/// This function calculates v-coefficients
#[allow(unused)]
fn calculate_v_coefficients(a: &DMatrix<Number>) -> DVector<Number> {
    let mut v_arr = DVector::<Number>::from_vec(vec![0.0; a.nrows()]);
    v_arr[0] = a[(0, 1)] / (-a[(0, 0)]);
    for i in 1..a.nrows() - 1 {
        v_arr[i] = a[(i, i + 1)] / (-a[(i, i)] - a[(i, i - 1)] * v_arr[i - 1]);
    }
    v_arr[a.nrows() - 1] = 0.0;
    return v_arr;
}

/// This function calculates u-coefficients.
/// It depends on v-coefficients
#[allow(unused)]
fn calculate_u_coefficients(
    a: &DMatrix<Number>,
    b: &DVector<Number>,
    v_arr: &DVector<Number>,
) -> DVector<Number> {
    let mut u_arr = DVector::<Number>::from_vec(vec![0.0; a.nrows()]);
    u_arr[0] = -b[0] / (-a[(0, 0)]);
    for i in 1..a.nrows() - 1 {
        u_arr[i] =
            (a[(i, i - 1)] * u_arr[i - 1] - b[i]) / (-a[(i, i)] - a[(i, i - 1)] * v_arr[i - 1]);
    }
    u_arr[a.nrows() - 1] = (a[(a.nrows() - 1, a.nrows() - 2)] * u_arr[a.nrows() - 2]
        - b[a.nrows() - 1])
        / (-a[(a.nrows() - 1, a.nrows() - 1)]
            - a[(a.nrows() - 1, a.nrows() - 2)] * v_arr[a.nrows() - 2]);
    return u_arr;
}

/// This function checks if two vectors are almost equal
///
/// We need it because nalgebra::DVector does not have almost_equal function.
/// Generally 0.99999999999 and 1.00000000001 are not equal completely, so we need to use almost_equal
#[allow(unused)]
pub fn vectors_almost_equal(
    vec1: &DVector<Number>,
    vec2: &DVector<Number>,
    epsilon: Number,
) -> bool {
    if vec1.len() != vec2.len() {
        return false;
    }

    for i in 0..vec1.len() {
        if (vec1[i] - vec2[i]).abs() > epsilon {
            return false;
        }
    }

    true
}

/// This function solves the system of linear equations
/// Ax = b
/// It uses tridiagonal matrix method
/// Make sure the matrix is tri-diagonally dominant
#[allow(unused)]
pub fn solve(a: &DMatrix<Number>, b: &DVector<Number>) -> DVector<Number> {
    check_if_three_diagonals(a);
    let mut v_arr = calculate_v_coefficients(a);
    let mut u_arr = calculate_u_coefficients(a, b, &v_arr);

    let mut x_arr = DVector::<Number>::zeros(a.ncols());
    x_arr[a.nrows() - 1] = u_arr[a.nrows() - 1];
    for i in (0..a.nrows() - 1).rev() {
        x_arr[i] = u_arr[i] + v_arr[i] * x_arr[i + 1];
    }

    return x_arr;
}

/// Here are just some tests
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use nalgebra::{DMatrix, RowDVector};

    /// This test checks if the matrix is tri-diagonally dominant
    #[test]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__correct_matrix() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![2., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    /// This test checks if the matrix is not tri-diagonally dominant and should panic because matrix is not even square
    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_1() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 3., 0., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    /// This test checks if the matrix is not tri-diagonally dominant and should panic
    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_2() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![3., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 1., 3., 1., 2.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    /// This test checks if the matrix is not tri-diagonally dominant and should panic
    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_3() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![3., 1., 2., 1., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 1.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    #[test]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__correct_matrix_2() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    #[test]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__correct_matrix_3() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![3., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 4., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_4() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., -1.]),
            RowDVector::from_vec(vec![3., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    /// Here are some tests directly for the solver
    mod solver {
        use super::*;
        use rand::Rng;

        /// This function calculates b to prepare data for text
        #[allow(unused)]
        fn calculate_b(a: &DMatrix<Number>, x_arr: &DVector<Number>) -> DVector<Number> {
            return a * x_arr;
        }

        /// Here is we print the equation that about to be solved
        #[allow(unused)]
        fn print_equation(a: &DMatrix<Number>, x: &DVector<Number>, b: &DVector<Number>) {
            print!("A:{}", a);
            println!("*");
            print!("x:{}", x);
            println!("=");
            print!("b:{}", b);
        }

        /// This test is used to avoid repeating of code and call the solver with different input parameters
        #[allow(unused)]
        fn base_solve_test(a: &DMatrix<Number>, example_x: &DVector<Number>) {
            let example_b = calculate_b(&a, &example_x);
            print_equation(&a, &example_x, &example_b);
            let solve_x = solve(&a, &example_b);
            println!("Calculated X: {}", solve_x.to_string());
            println!("Expected X: {}", example_x.to_string());
            println!("----------------------------");
            assert!(vectors_almost_equal(&solve_x, example_x, 0.001));
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_1() {
            let a = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
                RowDVector::from_vec(vec![2., 1., 2., 0., 0.]),
                RowDVector::from_vec(vec![0., 3., 1., 2., 0.]),
                RowDVector::from_vec(vec![0., 0., 4., 1., 2.]),
                RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 3., 2., 5., 4.]);
            base_solve_test(&a, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_2() {
            let a = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![4., 2., 0., 0.]),
                RowDVector::from_vec(vec![1., 3., 1., 0.]),
                RowDVector::from_vec(vec![0., 1., 3., 2.]),
                RowDVector::from_vec(vec![0., 0., 1., 4.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 2., 3., 4.]);
            base_solve_test(&a, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_3() {
            let a = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![2., 1., 0., 0., 0.]),
                RowDVector::from_vec(vec![1., 3., 1., 0., 0.]),
                RowDVector::from_vec(vec![0., 1., 3., 1., 0.]),
                RowDVector::from_vec(vec![0., 0., 1., 4., 1.]),
                RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![5., 4., 3., 2., 1.]);
            base_solve_test(&a, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_4() {
            let a = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![10., 2., 0.]),
                RowDVector::from_vec(vec![3., 8., 1.]),
                RowDVector::from_vec(vec![0., 4., 7.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![2., 1., 3.]);
            base_solve_test(&a, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_5() {
            let a = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![5., 1., 0.]),
                RowDVector::from_vec(vec![1., 6., 2.]),
                RowDVector::from_vec(vec![0., 2., 5.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 1., 1.]);
            base_solve_test(&a, &example_x);
        }

        /// Here we run solve tests with random matrices
        #[test]
        fn solve_random_matrix_multiple_tests() {
            println!("----------------------------");
            println!("Running random matrix tests...");
            let mut rng = rand::thread_rng();

            // Define the number of test runs
            let test_runs = 10; // You can adjust the number of iterations for larger tests

            for _ in 0..test_runs {
                let n = rng.gen_range(2..10); // Matrix size, adjust for larger matrices

                // Generate a random tridiagonal matrix
                let mut matrix = DMatrix::<Number>::zeros(n, n);
                for i in 0..n {
                    // Diagonal element
                    matrix[(i, i)] = rng.gen_range(1.0..10.0);

                    // Subdiagonal element (if applicable)
                    if i > 0 {
                        matrix[(i, i - 1)] = rng.gen_range(0.1..1.0);
                    }

                    // Superdiagonal element (if applicable)
                    if i < n - 1 {
                        matrix[(i, i + 1)] = rng.gen_range(0.1..1.0);
                    }
                }

                // Generate a random x vector
                let example_x = DVector::<Number>::from_fn(n, |_, _| rng.gen_range(1.0..10.0));

                // Call the base_solve_test function to perform the test
                base_solve_test(&matrix, &example_x);
            }
        }
    }
}
