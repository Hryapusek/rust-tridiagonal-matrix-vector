use core::panic;

use nalgebra::DMatrix;
use nalgebra::DVector;

pub type Number = f64;

fn check_if_three_diagonals(A: &DMatrix<Number>) -> () {
    for row in 0..A.nrows() {
        for col in 0..A.ncols() {
            if !(row.abs_diff(col) <= 1) && (A[(row, col)] != 0.0) {
                panic!("Matrix is not diagonally dominant");
            }
        }
    }
}

fn calculate_V_coefficients(A: &DMatrix<Number>) -> DVector<Number> {
    let mut v_arr = DVector::<Number>::from_vec(vec![0.0; A.nrows()]);
    v_arr[0] = A[(0, 1)] / (-A[(0, 0)]);
    for i in 1..A.nrows() - 1 {
        v_arr[i] = A[(i, i + 1)] / (-A[(i, i)] - A[(i, i - 1)] * v_arr[i - 1]);
    }
    v_arr[A.nrows() - 1] = 0.0;
    return v_arr;
}

fn calculate_U_coefficients(
    A: &DMatrix<Number>,
    b: &DVector<Number>,
    v_arr: &DVector<Number>,
) -> DVector<Number> {
    let mut u_arr = DVector::<Number>::from_vec(vec![0.0; A.nrows()]);
    u_arr[0] = -b[0] / (-A[(0, 0)]);
    for i in 1..A.nrows() - 1 {
        u_arr[i] =
            (A[(i, i - 1)] * u_arr[i - 1] - b[i]) / (-A[(i, i)] - A[(i, i - 1)] * v_arr[i - 1]);
    }
    u_arr[A.nrows() - 1] = (A[(A.nrows() - 1, A.nrows() - 2)] * u_arr[A.nrows() - 2]
        - b[A.nrows() - 1])
        / (-A[(A.nrows() - 1, A.nrows() - 1)]
            - A[(A.nrows() - 1, A.nrows() - 2)] * v_arr[A.nrows() - 2]);
    return u_arr;
}

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

pub fn solve(A: &DMatrix<Number>, b: &DVector<Number>) -> DVector<Number> {
    check_if_three_diagonals(A);
    if A.nrows() != A.ncols() {
        panic!("The matrix is not square");
    }
    let mut v_arr = calculate_V_coefficients(A);
    let mut u_arr = calculate_U_coefficients(A, b, &v_arr);

    let mut x_arr = DVector::<Number>::zeros(A.ncols());
    x_arr[A.nrows() - 1] = u_arr[A.nrows() - 1];
    for i in (0..A.nrows() - 1).rev() {
        x_arr[i] = u_arr[i] + v_arr[i] * x_arr[i + 1];
    }

    return x_arr;
}

mod tests {
    use super::check_if_three_diagonals;
    use super::*;
    use nalgebra::{DMatrix, RowDVector};

    #[test]
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

        fn calculate_b(A: &DMatrix<Number>, x_arr: &DVector<Number>) -> DVector<Number> {
            return A * x_arr;
        }

        fn print_equation(A: &DMatrix<Number>, x: &DVector<Number>, b: &DVector<Number>) {
            for i in 0..A.nrows() {
                print!("| ");
                for j in 0..A.ncols() {
                    print!("{} ", A[(i, j)]);
                }
                print!("| * | {} ", x[i]);
                print!("= | {} \n", b[i]);
            }
        }

        fn base_solve_test(example_matrix: &DMatrix<Number>, example_x: &DVector<Number>) {
            let example_b = calculate_b(&example_matrix, &example_x);
            print_equation(&example_matrix, &example_x, &example_b);
            let solve_x = solve(&example_matrix, &example_b);
            println!("Result: {}", solve_x.to_string());
            println!("Expected: {}", example_x.to_string());
            assert!(vectors_almost_equal(&solve_x, example_x, 0.001));
        }

        #[test]
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
