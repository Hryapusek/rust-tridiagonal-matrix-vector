use crate::math::solver::vectors_almost_equal;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;

trait TridiagonalMatrixSize {
    fn size(&self) -> usize;
}
trait TridiagonalNumberType<NumberType>:
    std::convert::From<i32>
    + Clone
    + std::cmp::PartialEq
    + std::cmp::PartialOrd
    + std::ops::Add<Output = NumberType>
    + std::ops::Sub<Output = NumberType>
    + std::ops::Mul<Output = NumberType>
    + std::ops::Div<Output = NumberType>
    + std::ops::AddAssign
    + std::ops::SubAssign
    + std::ops::MulAssign
    + std::ops::DivAssign
    + std::ops::Neg
    + std::ops::Add<NumberType, Output = NumberType>
    + std::ops::Sub<NumberType, Output = NumberType>
    + std::ops::Mul<NumberType, Output = NumberType>
    + std::ops::Div<NumberType, Output = NumberType>
    + std::ops::AddAssign<NumberType>
    + std::ops::SubAssign<NumberType>
    + std::ops::MulAssign<NumberType>
    + std::ops::DivAssign<NumberType>
    + std::ops::Neg<Output = NumberType>
{
}

impl TridiagonalNumberType<i32> for i32 {}
impl TridiagonalNumberType<f64> for f64 {}

pub struct TridiagonalSparseMatrix<NumberType: TridiagonalNumberType<NumberType>> {
    pub subdiagonal: Vec<NumberType>,
    pub maindiagonal: Vec<NumberType>,
    pub updiagonal: Vec<NumberType>,
    pub zero: NumberType,
}

impl<'a, NumberType> Mul<&Vec<NumberType>> for &'a TridiagonalSparseMatrix<NumberType>
where
    NumberType: TridiagonalNumberType<NumberType>,
{
    type Output = Vec<NumberType>;

    fn mul(self, vec: &Vec<NumberType>) -> Self::Output {
        let n = self.size();
        let mut result = vec![self.zero.clone(); n];

        for i in 0..n {
            result[i] = self.maindiagonal[i].clone() * vec[i].clone();
            if i > 0 {
                result[i] =
                    result[i].clone() + self.subdiagonal[i - 1].clone() * vec[i - 1].clone();
            }
            if i < n - 1 {
                result[i] = result[i].clone() + self.updiagonal[i].clone() * vec[i + 1].clone();
            }
        }

        result
    }
}

impl<NumberType: TridiagonalNumberType<NumberType>> TridiagonalSparseMatrix<NumberType> {
    fn new(n: usize) -> TridiagonalSparseMatrix<NumberType> {
        TridiagonalSparseMatrix {
            subdiagonal: vec![NumberType::from(0); n - 1],
            maindiagonal: vec![NumberType::from(0); n],
            updiagonal: vec![NumberType::from(0); n - 1],
            zero: NumberType::from(0),
        }
    }

    #[allow(unused)]
    fn check_if_three_diagonals(&self) {
        let n = self.size();
        if n == 0 {
            panic!("The matrix has zero size");
        }

        // Iterate over each row
        for i in 0..n {
            let main_value = self.maindiagonal[i].clone();

            // Check if the main diagonal element is zero (it shouldn't be)
            if main_value == self.zero {
                panic!("Matrix is not diagonally dominant: zero on the main diagonal");
            }

            // Check the condition that the sum of the subdiagonal and updiagonal elements
            // is less than or equal to the main diagonal element
            let mut sum_of_neighbors = self.zero.clone();

            if i > 0 {
                let sub_value = self.subdiagonal[i - 1].clone();
                sum_of_neighbors = sum_of_neighbors + sub_value;
            }

            if i < n - 1 {
                let up_value = self.updiagonal[i].clone();
                sum_of_neighbors = sum_of_neighbors + up_value;
            }

            if sum_of_neighbors > main_value {
                panic!(
                    "Matrix is not diagonally dominant: the sum of subdiagonal and updiagonal elements at row {} is greater than the main diagonal element",
                    i
                );
            }
        }
    }

    /// This function calculates the v-coefficients for the sparse matrix
    fn calculate_v_coefficients(&self) -> Vec<NumberType> {
        let n = self.size();
        let mut v_arr = vec![self.zero.clone(); n];
        if n > 1 {
            v_arr[0] = self.updiagonal[0].clone() / (-self.maindiagonal[0].clone());
            for i in 1..n - 1 {
                let denom = -self.maindiagonal[i].clone()
                    - self.subdiagonal[i - 1].clone() * v_arr[i - 1].clone();
                v_arr[i] = self.updiagonal[i].clone() / denom;
            }
            v_arr[n - 1] = self.zero.clone();
        }
        v_arr
    }

    /// This function calculates the u-coefficients for the sparse matrix
    fn calculate_u_coefficients(
        &self,
        b: &Vec<NumberType>,
        v_arr: &Vec<NumberType>,
    ) -> Vec<NumberType> {
        let n = self.size();
        let mut u_arr = vec![self.zero.clone(); n];
        if n > 0 {
            u_arr[0] = -b[0].clone() / (-self.maindiagonal[0].clone());
            for i in 1..n - 1 {
                let denom = -self.maindiagonal[i].clone()
                    - self.subdiagonal[i - 1].clone() * v_arr[i - 1].clone();
                u_arr[i] =
                    (self.subdiagonal[i - 1].clone() * u_arr[i - 1].clone() - b[i].clone()) / denom;
            }
            let denom_last = -self.maindiagonal[n - 1].clone()
                - self.subdiagonal[n - 2].clone() * v_arr[n - 2].clone();
            u_arr[n - 1] = (self.subdiagonal[n - 2].clone() * u_arr[n - 2].clone()
                - b[n - 1].clone())
                / denom_last;
        }
        u_arr
    }

    /// This function solves the system of linear equations Ax = b
    /// using the tridiagonal matrix algorithm (Thomas algorithm).
    pub fn solve(&self, b: &Vec<NumberType>) -> Vec<NumberType> {
        // Check if the matrix is diagonally dominant
        self.check_if_three_diagonals();

        // Calculate the v and u coefficients
        let v_arr = self.calculate_v_coefficients();
        let u_arr = self.calculate_u_coefficients(b, &v_arr);

        // Back substitution to solve for x
        let n = self.size();
        let mut x_arr = vec![self.zero.clone(); n];
        if n > 0 {
            x_arr[n - 1] = u_arr[n - 1].clone();
            for i in (0..n - 1).rev() {
                x_arr[i] = u_arr[i].clone() + v_arr[i].clone() * x_arr[i + 1].clone();
            }
        }

        x_arr
    }
}

impl<'a, NumberType: TridiagonalNumberType<NumberType>> Index<(usize, usize)>
    for TridiagonalSparseMatrix<NumberType>
{
    type Output = NumberType;
    fn index(&self, idx: (usize, usize)) -> &NumberType {
        let (x, y) = idx;
        match (x, y) {
            (x, y) if x == y => self.maindiagonal.get(x).unwrap(),
            (x, y) if x == y - 1 => self.subdiagonal.get(x).unwrap(),
            (x, y) if x == y + 1 => self.updiagonal.get(x).unwrap(),
            _ => &self.zero,
        }
    }
}

impl<'a, NumberType: TridiagonalNumberType<NumberType>> IndexMut<(usize, usize)>
    for TridiagonalSparseMatrix<NumberType>
{
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut NumberType {
        let (x, y) = idx;
        match (x, y) {
            (x, y) if x == y => self.maindiagonal.get_mut(x).unwrap(),
            (x, y) if x == y - 1 => self.subdiagonal.get_mut(x).unwrap(),
            (x, y) if x == y + 1 => self.updiagonal.get_mut(x).unwrap(),
            _ => panic!("Index out of bounds"),
        }
    }
}

impl<'a, NumberType: TridiagonalNumberType<NumberType>> TridiagonalMatrixSize
    for TridiagonalSparseMatrix<NumberType>
{
    fn size(&self) -> usize {
        self.maindiagonal.len()
    }
}

#[cfg(test)]
mod tests {
    use crate::math::solver::solve;

    use super::*;
    use nalgebra::{DMatrix, DVector};
    use rand::Rng;
    use std::ops::Mul;

    // A helper function to create a simple tridiagonal matrix for testing
    fn create_test_matrix(
        main_diagonal: Vec<f64>,
        sub_diagonal: Vec<f64>,
        up_diagonal: Vec<f64>,
    ) -> TridiagonalSparseMatrix<f64> {
        TridiagonalSparseMatrix {
            maindiagonal: main_diagonal,
            subdiagonal: sub_diagonal,
            updiagonal: up_diagonal,
            zero: 0.0,
        }
    }

    #[test]
    fn test_valid_tridiagonal_matrix() {
        // A valid tridiagonal matrix where the main diagonal elements are larger
        // than the sum of their neighbors
        let main_diagonal = vec![4.0, 5.0, 6.0];
        let sub_diagonal = vec![1.0, 1.0];
        let up_diagonal = vec![1.0, 1.0];

        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);
        matrix.check_if_three_diagonals(); // Should not panic
    }

    #[test]
    #[should_panic(expected = "Matrix is not diagonally dominant")]
    fn test_invalid_tridiagonal_matrix() {
        // An invalid tridiagonal matrix where the main diagonal elements are smaller
        // than the sum of their neighbors
        let main_diagonal = vec![2.0, 2.0, 2.0];
        let sub_diagonal = vec![1.5, 1.5];
        let up_diagonal = vec![1.5, 1.5];

        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);
        matrix.check_if_three_diagonals(); // Should panic
    }

    #[test]
    fn test_single_element_matrix() {
        // A matrix with a single element, which should be considered diagonally dominant
        let main_diagonal = vec![1.0];
        let sub_diagonal = vec![];
        let up_diagonal = vec![];

        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);
        matrix.check_if_three_diagonals(); // Should not panic
    }

    #[test]
    #[should_panic]
    fn test_zero_size_matrix() {
        // A matrix with no elements, which should panic
        let main_diagonal = vec![];
        let sub_diagonal = vec![];
        let up_diagonal = vec![];

        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);
        matrix.check_if_three_diagonals(); // Should panic
    }

    #[test]
    fn test_basic_multiplication() {
        // A simple 3x3 tridiagonal matrix
        let main_diagonal = vec![4.0, 5.0, 6.0];
        let sub_diagonal = vec![1.0, 1.0];
        let up_diagonal = vec![1.0, 1.0];
        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);

        // Vector to multiply with
        let vector: Vec<f64> = vec![1.0, 2.0, 3.0];

        // Expected result is calculated as:
        // [4*1 + 1*2, 1*1 + 5*2 + 1*3, 1*2 + 6*3] = [6, 12, 20]
        let expected_result = vec![6.0, 14.0, 20.0];

        // Perform the multiplication
        let result = matrix.mul(&vector);

        // Check if the result matches the expected values
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_single_element_multiplication() {
        // A 1x1 tridiagonal matrix
        let main_diagonal = vec![2.0];
        let sub_diagonal = vec![];
        let up_diagonal = vec![];
        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);

        // Vector to multiply with
        let vector = vec![3.0];

        // Expected result: [2 * 3] = [6]
        let expected_result = vec![6.0];

        // Perform the multiplication
        let result = matrix.mul(&vector);

        // Check if the result matches the expected values
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_zero_size_multiplication() {
        // A matrix with zero size
        let main_diagonal = vec![];
        let sub_diagonal = vec![];
        let up_diagonal = vec![];
        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);

        // Empty vector to multiply with
        let vector: Vec<f64> = vec![];

        // Expected result: empty vector
        let expected_result: Vec<f64> = vec![];

        // Perform the multiplication
        let result = matrix.mul(&vector);

        // Check if the result matches the expected values
        assert_eq!(result, expected_result);
    }

    #[test]
    fn test_larger_matrix_multiplication() {
        // A 4x4 tridiagonal matrix
        let main_diagonal = vec![4.0, 5.0, 6.0, 7.0];
        let sub_diagonal = vec![1.0, 1.0, 1.0];
        let up_diagonal = vec![1.0, 1.0, 1.0];
        let matrix = create_test_matrix(main_diagonal, sub_diagonal, up_diagonal);

        // Vector to multiply with
        let vector = vec![1.0, 2.0, 3.0, 4.0];

        // Expected result is calculated as:
        // [4*1 + 1*2, 1*1 + 5*2 + 1*3, 1*2 + 6*3 + 1*4, 1*3 + 7*4]
        let expected_result = vec![6.0, 14.0, 24.0, 31.0];

        // Perform the multiplication
        let result = matrix.mul(&vector);

        // Check if the result matches the expected values
        assert_eq!(result, expected_result);
    }

    /// Function to calculate b = A * x for the sparse matrix
    #[allow(unused)]
    fn calculate_b_sparse(a: &TridiagonalSparseMatrix<f64>, x_arr: &Vec<f64>) -> Vec<f64> {
        a * x_arr
    }

    /// Print the equation for the sparse matrix
    #[allow(unused)]
    fn print_equation_sparse(a: &TridiagonalSparseMatrix<f64>, x: &Vec<f64>, b: &Vec<f64>) {
        println!("Tridiagonal Matrix A:");
        println!("Main Diagonal: {:?}", a.maindiagonal);
        println!("Sub Diagonal: {:?}", a.subdiagonal);
        println!("Up Diagonal: {:?}", a.updiagonal);
        println!("x: {:?}", x);
        println!("b: {:?}", b);
    }

    /// Base test function to solve the sparse matrix system
    #[allow(unused)]
    fn base_solve_test_sparse(a: &TridiagonalSparseMatrix<f64>, example_x: &Vec<f64>) {
        let example_b = calculate_b_sparse(a, example_x);
        print_equation_sparse(a, example_x, &example_b);
        let solve_x = a.solve(&example_b);
        println!("Calculated X: {:?}", solve_x);
        println!("Expected X: {:?}", example_x);
        println!("----------------------------");
        assert!(vectors_almost_equal(
            &DVector::from_vec(solve_x),
            &DVector::from_vec(example_x.clone()),
            0.001
        ));
    }

    #[test]
    #[allow(non_snake_case)]
    fn solve__correct_matrix_1_sparse() {
        // A sparse tridiagonal matrix with size 5x5
        let main_diagonal = vec![5.0, 5.0, 5.0, 5.0, 5.0];
        let sub_diagonal = vec![1.0, 1.0, 1.0, 1.0];
        let up_diagonal = vec![1.0, 1.0, 1.0, 1.0];

        let matrix = TridiagonalSparseMatrix {
            maindiagonal: main_diagonal,
            subdiagonal: sub_diagonal,
            updiagonal: up_diagonal,
            zero: 0.0,
        };

        // Example solution vector x
        let example_x = vec![1.0, 3.0, 2.0, 5.0, 4.0];

        // Test the solver
        base_solve_test_sparse(&matrix, &example_x);
    }

    /// Here we run solve tests with random tridiagonal sparse matrices
    #[test]
    fn solve_random_matrix_multiple_tests_sparse() {
        println!("----------------------------");
        println!("Running random sparse matrix tests...");
        let mut rng = rand::thread_rng();

        // Define the number of test runs
        let test_runs = 10; // You can adjust the number of iterations for larger tests

        for _ in 0..test_runs {
            let n = rng.gen_range(2..10); // Matrix size, adjust for larger matrices

            // Generate random main, sub, and superdiagonals for the tridiagonal matrix
            let mut main_diagonal = vec![0.0; n];
            let mut sub_diagonal = vec![0.0; n - 1];
            let mut up_diagonal = vec![0.0; n - 1];

            for i in 0..n {
                // Generate the main diagonal element
                main_diagonal[i] = rng.gen_range(2.1..10.0);

                // Generate the subdiagonal element if applicable
                if i > 0 {
                    sub_diagonal[i - 1] = rng.gen_range(0.1..1.0);
                }

                // Generate the updiagonal element if applicable
                if i < n - 1 {
                    up_diagonal[i] = rng.gen_range(0.1..1.0);
                }
            }

            // Create the sparse matrix
            let matrix = TridiagonalSparseMatrix {
                maindiagonal: main_diagonal,
                subdiagonal: sub_diagonal,
                updiagonal: up_diagonal,
                zero: 0.0,
            };

            // Generate a random example x vector
            let example_x: Vec<f64> = (0..n).map(|_| rng.gen_range(1.0..10.0)).collect();

            // Perform the test using the base_solve_test_sparse function
            base_solve_test_sparse(&matrix, &example_x);
        }
    }

    use std::time::Instant;

    /// Test with a large matrix to compare the performance between dense and sparse representations
    #[test]
    fn solve_large_matrix_performance_test() {
        println!("----------------------------");
        println!("Running large matrix performance test...");

        let mut rng = rand::thread_rng();

        // Define the size of the matrix (large size for performance testing)
        let n = 100000; // Adjust the size for even larger tests if desired

        // Generate random main, sub, and superdiagonals for the tridiagonal matrix
        let mut main_diagonal = vec![0.0; n];
        let mut sub_diagonal = vec![0.0; n - 1];
        let mut up_diagonal = vec![0.0; n - 1];

        for i in 0..n {
            // Generate the main diagonal element
            main_diagonal[i] = rng.gen_range(2.1..10.0);

            // Generate the subdiagonal element if applicable
            if i > 0 {
                sub_diagonal[i - 1] = rng.gen_range(0.1..1.0);
            }

            // Generate the updiagonal element if applicable
            if i < n - 1 {
                up_diagonal[i] = rng.gen_range(0.1..1.0);
            }
        }

        // Create the sparse matrix
        let sparse_matrix = TridiagonalSparseMatrix {
            maindiagonal: main_diagonal.clone(),
            subdiagonal: sub_diagonal.clone(),
            updiagonal: up_diagonal.clone(),
            zero: 0.0,
        };

        // Generate a random example x vector
        let example_x: Vec<f64> = (0..n).map(|_| rng.gen_range(1.0..10.0)).collect();

        // Calculate b = A * x using the sparse matrix
        let example_b_sparse = calculate_b_sparse(&sparse_matrix, &example_x);

        // Measure the time taken to solve using the sparse matrix approach
        let start_sparse = Instant::now();
        let solve_x_sparse = sparse_matrix.solve(&example_b_sparse);
        let duration_sparse = start_sparse.elapsed();

        println!("Solution using sparse matrix(first five elements): {:?}", solve_x_sparse.iter().take(5usize).collect::<Vec<_>>());

        println!(
            "Sparse matrix solution time: {:.2?} for size {}",
            duration_sparse, n
        );
    }
}
