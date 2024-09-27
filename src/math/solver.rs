use core::panic;

use nalgebra::DMatrix;
use nalgebra::DVector;

#[allow(unused)]
pub type Number = f64;

#[allow(unused)]
fn check_if_three_diagonals(a: &DMatrix<Number>) -> () {
    for row in 0..a.nrows() {
        for col in 0..a.ncols() {
            if !(row.abs_diff(col) <= 1) && (a[(row, col)] != 0.0) {
                panic!("Matrix is not diagonally dominant");
            }
        }
    }
}

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

#[allow(unused)]
pub fn solve(a: &DMatrix<Number>, b: &DVector<Number>) -> DVector<Number> {
    check_if_three_diagonals(a);
    if a.nrows() != a.ncols() {
        panic!("The matrix is not square");
    }
    let mut v_arr = calculate_v_coefficients(a);
    let mut u_arr = calculate_u_coefficients(a, b, &v_arr);

    let mut x_arr = DVector::<Number>::zeros(a.ncols());
    x_arr[a.nrows() - 1] = u_arr[a.nrows() - 1];
    for i in (0..a.nrows() - 1).rev() {
        x_arr[i] = u_arr[i] + v_arr[i] * x_arr[i + 1];
    }

    return x_arr;
}

mod tests {
    use super::*;
    #[allow(unused_imports)]
    use nalgebra::{DMatrix, RowDVector};

    #[test]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__correct_matrix() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

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

    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_2() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![1., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    #[test]
    #[should_panic]
    #[allow(non_snake_case)]
    fn check_if_three_diagonals__incorrect_matrix_3() {
        let example_matrix = DMatrix::<Number>::from_rows(&[
            RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
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
            RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
            RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
        ]);
        check_if_three_diagonals(&example_matrix);
    }

    mod solver {
        use super::*;

        #[allow(unused)]
        fn calculate_b(a: &DMatrix<Number>, x_arr: &DVector<Number>) -> DVector<Number> {
            return a * x_arr;
        }

        #[allow(unused)]
        fn print_equation(a: &DMatrix<Number>, x: &DVector<Number>, b: &DVector<Number>) {
            print!("A:{}", a);
            println!("*");
            print!("x:{}", x);
            println!("=");
            print!("b:{}", b);
        }

        #[allow(unused)]
        fn base_solve_test(example_matrix: &DMatrix<Number>, example_x: &DVector<Number>) {
            let example_b = calculate_b(&example_matrix, &example_x);
            print_equation(&example_matrix, &example_x, &example_b);
            let solve_x = solve(&example_matrix, &example_b);
            println!("Result: {}", solve_x.to_string());
            println!("Expected: {}", example_x.to_string());
            assert!(vectors_almost_equal(&solve_x, example_x, 0.001));
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_1() {
            let example_matrix = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
                RowDVector::from_vec(vec![2., 1., 2., 0., 0.]),
                RowDVector::from_vec(vec![0., 3., 1., 2., 0.]),
                RowDVector::from_vec(vec![0., 0., 4., 1., 2.]),
                RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 3., 2., 5., 4.]);
            base_solve_test(&example_matrix, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_2() {
            let example_matrix = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![4., 2., 0., 0.]),
                RowDVector::from_vec(vec![1., 3., 1., 0.]),
                RowDVector::from_vec(vec![0., 1., 3., 2.]),
                RowDVector::from_vec(vec![0., 0., 1., 4.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 2., 3., 4.]);
            base_solve_test(&example_matrix, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_3() {
            let example_matrix = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![2., 1., 0., 0., 0.]),
                RowDVector::from_vec(vec![1., 3., 1., 0., 0.]),
                RowDVector::from_vec(vec![0., 1., 3., 1., 0.]),
                RowDVector::from_vec(vec![0., 0., 1., 4., 1.]),
                RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![5., 4., 3., 2., 1.]);
            base_solve_test(&example_matrix, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_4() {
            let example_matrix = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![10., 2., 0.]),
                RowDVector::from_vec(vec![3., 8., 1.]),
                RowDVector::from_vec(vec![0., 4., 7.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![2., 1., 3.]);
            base_solve_test(&example_matrix, &example_x);
        }

        #[test]
        #[allow(non_snake_case)]
        fn solve__correct_matrix_5() {
            let example_matrix = DMatrix::<Number>::from_rows(&[
                RowDVector::from_vec(vec![5., 1., 0.]),
                RowDVector::from_vec(vec![1., 6., 2.]),
                RowDVector::from_vec(vec![0., 2., 5.]),
            ]);
            let example_x = DVector::<Number>::from_vec(vec![1., 1., 1.]);
            base_solve_test(&example_matrix, &example_x);
        }
    }
}
