# Task 1 on Voskoboynikov math - Tridigonal matrix method

## Build & Run
In rust you dont have to worry about nothing. Just install rust on your computer, go to the folder with project and type
```sh
cargo test -- --color always
```
If you want to see more output - use this:
```sh
cargo test -- --nocapture --color always
```
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

If we go a little deeper in the file - we will find different tests for this and other private functions.