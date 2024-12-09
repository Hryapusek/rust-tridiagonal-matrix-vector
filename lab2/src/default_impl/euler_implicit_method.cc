#include <default_impl/euler_implicit_method.hpp>

#include <vector>
#include <stdexcept>

#include <contract/contract.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

// Helper function to solve a tridiagonal system using the Thomas algorithm
void solve_tridiagonal(
    const Eigen::VectorXd& a,   // Lower diagonal (size n-1)
    const Eigen::VectorXd& b,   // Main diagonal (size n)
    const Eigen::VectorXd& c,   // Upper diagonal (size n-1)
    const Eigen::VectorXd& d,   // Right-hand side (size n)
    Eigen::VectorXd& x          // Solution (size n)
)
{
    size_t n = b.size();
    
    Eigen::VectorXd cp(n);  // Temporary vector for forward substitution
    Eigen::VectorXd dp(n);  // Temporary vector for the right-hand side adjustment

    // Forward elimination
    cp(0) = c(0) / b(0);
    dp(0) = d(0) / b(0);

    for (size_t i = 1; i < n - 1; ++i) {  // Loop from 1 to n-2 for the main part
        double m = 1.0 / (b(i) - a(i - 1) * cp(i - 1));  // a(i-1) for lower diagonal
        cp(i) = c(i) * m;
        dp(i) = (d(i) - a(i - 1) * dp(i - 1)) * m;
    }

    // Last step for the forward elimination
    dp(n - 1) = (d(n - 1) - a(n - 2) * dp(n - 2)) / (b(n - 1) - a(n - 2) * cp(n - 2));

    // Backward substitution
    x(n - 1) = dp(n - 1);
    for (size_t i = n - 2; i != static_cast<size_t>(-1); --i) {
        x(i) = dp(i) - cp(i) * x(i + 1);
    }
}

auto DefaultEulerImplicitMethod::integrate(
    Eigen::SparseVector<Number_t> const& start_v,
    Eigen::SparseMatrix<Number_t> const& A,
    Eigen::SparseVector<Number_t> const& g,
    std::vector<Number_t> const& points
) -> Eigen::SparseMatrix<Number_t>
{
    Eigen::SparseMatrix<Number_t> result(
        A.rows(), points.size());  // Adjust the columns to match intervals

    contract(fun)
    {
        precondition(A.rows() == A.cols(), "A must be a square matrix");
        precondition(A.rows() == start_v.rows(), "A and start_v must have the same number of rows");
        precondition(A.rows() == g.rows(), "A and g must have the same number of rows");
        postcondition(
            result.cols() == points.size(),
            "result must have the same number of columns as intervals"
        );
    };

    result.col(0) = start_v;

    auto E = Eigen::SparseMatrix<Number_t>(A.rows(), A.rows());
    E.setIdentity();

    for (size_t i = 1; i < points.size(); ++i) {  // Iterate over intervals
        auto H = points.at(i) - points.at(i - 1);  // Step size
        Eigen::SparseMatrix<Number_t> M = E - H * A;  // Matrix for the implicit step

        // Prepare the right-hand side of the equation
        Eigen::SparseVector<Number_t> rhs = result.col(i - 1) + H * g;

        // We will now extract the tridiagonal parts of M to solve the system
        size_t n = M.rows();
        
        // Extract diagonals from sparse matrix M
        Eigen::VectorXd a(n - 1);  // Lower diagonal
        Eigen::VectorXd b(n);      // Main diagonal
        Eigen::VectorXd c(n - 1);  // Upper diagonal

        // Fill diagonals by iterating through the sparse matrix
        for (int k = 0; k < n; ++k) {
            b(k) = M.coeff(k, k);  // Main diagonal
            if (k > 0) {
                a(k - 1) = M.coeff(k, k - 1);  // Lower diagonal
            }
            if (k < n - 1) {
                c(k) = M.coeff(k, k + 1);  // Upper diagonal
            }
        }

        // Prepare the right-hand side (rhs) vector
        Eigen::VectorXd rhs_full(n);
        for (size_t k = 0; k < n; ++k) {
            rhs_full(k) = rhs.coeff(k);  // Copy rhs values
        }

        // Now solve the tridiagonal system
        Eigen::VectorXd solution(n);
        solve_tridiagonal(a, b, c, rhs_full, solution);

        // Update the result for this timestep
        result.col(i) = solution.sparseView();
    }

    return result;
}

