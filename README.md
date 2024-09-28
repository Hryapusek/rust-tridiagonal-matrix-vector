# Task 1 on Voskoboynikov math - Tridiagonal matrix method
## Description
The goal of tridiagonal method is to solve equation like `Ax = b`, where `A` is a *tridiagonal matrix*, `b` is a *vector*, and `x` is a desired *solution vector*.
```
A:
  │ 1 2 0 0 0 │
  │ 2 1 2 0 0 │
  │ 0 3 1 2 0 │
  │ 0 0 4 1 2 │
  │ 0 0 0 1 2 │

*
x:
  │ x1 │
  │ x2 │
  │ x3 │
  │ x4 │
  │ x5 │

=
b:
  │  7 │
  │  9 │
  │ 21 │
  │ 21 │
  │ 13 │
```
## Build & Run
In rust you dont have to worry about anything. Just install rust on your computer, go to the folder with project and type
```sh
cargo test -- --color always
```
If you want to see more output - use this:
```sh
cargo test -- --nocapture --color always
```
The second command will print the equations in the terminal.

You will see several examples of equations that are getting solved using tridiagonal matrix method.

The main function is placed in `src/math/solver.rs` with name `solve`
```rust
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
```
First we check if the matrix is tridiagonal and square. If not we panic. Then we calculate the v coefficients. Then we calculate the u coefficients. Then we solve the system of equations using the tridiagonal method. You can read about tridiagonal methods in [here](https://dzen.ru/a/YDWWQaMy3XNzjwmJ)

> If we go a little deeper in the file - we will find different tests for this and other private functions. You can check commentaries across this file `src/math/solver.rs` to get more information
