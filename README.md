# Task 1 on Voskoboynikov math
## Description
Readme in other branches are not updated. This is readme for CP3 branch. Coeff calculator implemented for other branches, yet it only tested for CP3.

## Goal of the task
The goal is to write program which would solve given derivative equation of second order. Main task is to find out how to calculate `a` `b` `c` and `g` coefficients. Part for `i=1..N-1` is the same for each variant. More interesting becoming when we get to the first and last indexes.

## Code
Main part to change for your variant is in two places. First is main.rs - you should put your own function to check accuracy. Second is coeff_calculator.rs - you should put your own function to calculate `a` `b` `c` and `g` coefficients.

There is also part of the code in sparse_matrixes.rs which is needed for Sergei Petrovich. Once ive come to him - he told me that code should not fall if we input huge matrixes. In main.rs i use matrixes from library, yet sparse matrix is works and also fully implements calculations that are needed for this task. You can show tests to the teacher to convince that code works.

## Coefficients
Few words about coefficients. For i = 1..N-1 coefficients are the same. For i = 0 coefficient is special and for i = N-1 coefficient is special too. If you have condition of the first type at the left or at the right - the corresponding side will contain two of the coefficients equal to `one` and `zero`, g will be equal to `v1` or `v2` depending on the side.

>Example:
> If we have condition of the first type at the left - we would result in coefficients `b = 0, c = 1, g = v1`.
> If we have condition of the third type at the right - we would result in more complex coefficients. You can open `docs/report.odt` to see what coefficients would look like.

## Accuracy
Read the `report.odt` to see how accuracy works and how to calculate it.

## How to create function to test
First you should choose `k` `q` coefficients. It could be just numbers or function. Then you should choose left and right threasholds. After that you should calculate coefficients for your threasholds functions. If you have third condition - you can choose `hi_1` randomly. Then you should put all your coefficients and threasholds in the functoin to find what `f` look like. Then you should put everything to the `main.rs`. Don't forget that we `compare` our `theoretical` u and our `calculated` u.

## Tridiagonal matrix method
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

### Main code
To run the main code and see the table - just type
```sh
cargo run
```

### Tridiagonal matrix method
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

First we check if the matrix is tridiagonal and square. If not we panic. Then we calculate the v coefficients. Then we calculate the u coefficients. Then we solve the system of equations using the tridiagonal method. You can read about tridiagonal methods in [here](https://dzen.ru/a/YDWWQaMy3XNzjwmJ)

> If we go a little deeper in the file - we will find different tests for this and other private functions. You can check commentaries across this file `src/math/solver.rs` to get more information

## Note
The code is written by Abraamyan Alexander at 13.11.2024. Code looking f... awful but still works and more than enough for polytech. At this moment i am working and dont really got time to make it better.
