use core::panic;

use nalgebra::DMatrix;
use nalgebra::DVector;

fn check_if_three_diagonals(A: &DMatrix<f64>) -> () {
  for row in 0..A.nrows() {
    for col in 0..A.ncols() {
      println!("-----------");
      println!("Row: {}", row);
      println!("Col: {}", col);
      println!("Value: {}", A[(row, col)]);
      println!("Condition 1: {}", row.abs_diff(col) <= 1);
      println!("Condition 2: {}", A[(row, col)] != 0.0);
      if !(row.abs_diff(col) <= 1) && (A[(row, col)] != 0.0) {
        panic!("Matrix is not diagonally dominant");
      }
    }
  }
}

pub fn solve(A: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
  check_if_three_diagonals(A);
  return DVector::<f64>::zeros(A.ncols());
}

mod tests {
  use nalgebra::{DMatrix, RowDVector};
  use super::check_if_three_diagonals;

  #[test]
  fn check_if_three_diagonals__correct_matrix() {
    let example_matrix = DMatrix::<f64>::from_rows(&[
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
    let example_matrix = DMatrix::<f64>::from_rows(&[
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
    let example_matrix = DMatrix::<f64>::from_rows(&[
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
    let example_matrix = DMatrix::<f64>::from_rows(&[
      RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
      RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
      RowDVector::from_vec(vec![0., 0., 1., 2., 1.]),
      RowDVector::from_vec(vec![0., 0., 0., 1., 2.]),
    ]);
    check_if_three_diagonals(&example_matrix);
  }

  #[test]
  fn check_if_three_diagonals__correct_matrix_2() {
    let example_matrix = DMatrix::<f64>::from_rows(&[
      RowDVector::from_vec(vec![1., 2., 0., 0., 0.]),
      RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
      RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
      RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
    ]);
    check_if_three_diagonals(&example_matrix);
  }

  #[test]
  #[should_panic]
  fn check_if_three_diagonals__incorrect_matrix_4() {
    let example_matrix = DMatrix::<f64>::from_rows(&[
      RowDVector::from_vec(vec![1., 2., 0., 0., -1.]),
      RowDVector::from_vec(vec![0., 1., 2., 0., 0.]),
      RowDVector::from_vec(vec![0., 0., 1., 2., 0.]),
      RowDVector::from_vec(vec![0., 0., 1., 1., 2.]),
    ]);
    check_if_three_diagonals(&example_matrix);
  }
}
