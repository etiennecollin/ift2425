use nalgebra::{DMatrix, DVector};

use crate::utils::{FuncMulti, function_vec, gradient_mat, table_formatter};

/// Computes the solution of a linear system using Cramer's rule.
///
/// # Complexity
///
/// - O(n!)
///
/// # Arguments
///
/// - `matrix`: A square matrix of size `n x n`.
/// - `vector`: A vector of size `n`.
///
/// # Returns
///
/// - The solution vector `x`.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn solve_cramer(
    matrix: &DMatrix<f64>,
    vector: &DVector<f64>,
) -> Result<DVector<f64>, &'static str> {
    let mut result = DVector::zeros(vector.len());
    let det = matrix.determinant();
    if det == 0.0 {
        return Err("The matrix is singular");
    }
    println!("╭───────────────────────────────");
    println!("│ Cramer");
    println!("├───────────────────────────────");
    println!("│ det(A) = {det}");

    for i in 0..vector.len() {
        let mut matrix_aux = matrix.clone();
        matrix_aux.set_column(i, vector);
        let det_aux = matrix_aux.determinant();
        if det_aux == 0.0 {
            return Err("The matrix is singular");
        }
        result[i] = det_aux / det;
        println!(
            "│ x_{} = det(A_{})/det(A) = {}/{} = {}",
            i + 1,
            i + 1,
            det_aux,
            det,
            result[i]
        );
    }
    println!("╰───────────────────────────────");
    Ok(result)
}

///putes the solution of a linear system using the triangle method.
///
/// # Complexity
///
/// - O(n)
///
/// # Arguments
///
/// - `matrix`: A square matrix of size `n x n`.
/// - `vector`: A vector of size `n`.
/// - `is_upper`: A boolean indicating if the matrix is upper triangular.
pub fn solve_triangular(
    matrix: &DMatrix<f64>,
    vector: &DVector<f64>,
    is_upper: bool,
) -> DVector<f64> {
    let mut result = DVector::zeros(vector.len());
    let n = vector.len();
    let row_iterator = match is_upper {
        true => (0..n).rev().collect::<Vec<_>>().into_iter(),
        false => (0..n).collect::<Vec<_>>().into_iter(),
    };
    for i in row_iterator {
        let mut sum = 0.0;

        let col_iterator = match is_upper {
            true => (i + 1..n).rev().collect::<Vec<_>>().into_iter(),
            false => (0..i).collect::<Vec<_>>().into_iter(),
        };

        for j in col_iterator {
            sum += matrix[(i, j)] * result[j];
        }
        result[i] = (vector[i] - sum) / matrix[(i, i)];
    }
    result
}

/// Computes the LU decomposition of a matrix without pivoting.
/// The input matrix is modified in-place and contains the packed LU decomposition.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn no_piv_lu_packed_mut(matrix: &mut DMatrix<f64>) -> Result<(), &'static str> {
    let n = matrix.nrows();

    if matrix.ncols() != n {
        return Err("The matrix must be square");
    }

    // Initialize the first row of U
    for j in 1..n {
        matrix[(0, j)] /= matrix[(0, 0)];
    }

    // Compute elements of L and U
    for j in 1..n {
        // Compute elements of L
        for i in j..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + matrix[(i, k)] * matrix[(k, j)]);
            matrix[(i, j)] -= sum;
        }

        // Compute elements of U
        for i in (j + 1)..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + matrix[(j, k)] * matrix[(k, i)]);
            matrix[(j, i)] = (matrix[(j, i)] - sum) / matrix[(j, j)];
        }
    }

    Ok(())
}

/// Computes the LU decomposition of a matrix without pivoting.
///
/// The returned matrix contains the packed LU decomposition.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn no_piv_lu_packed(matrix: &DMatrix<f64>) -> Result<DMatrix<f64>, &'static str> {
    let mut lu = matrix.clone();
    no_piv_lu_packed_mut(&mut lu)?;
    Ok(lu)
}

/// Computes the LU decomposition of a matrix without pivoting.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
///
/// # Returns
///
/// - A tuple containing the lower triangular matrix `L` and the upper triangular matrix `U`.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn no_piv_lu(matrix: &DMatrix<f64>) -> Result<(DMatrix<f64>, DMatrix<f64>), &'static str> {
    let lu = no_piv_lu_packed(matrix)?;
    Ok(lu_unpack(&lu))
}

/// Computes the LU decomposition of a matrix with partial pivoting and returns
/// the solution of the system.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
/// - `b`: A vector of size `n`.
///
/// # Errors
///
/// - If matrix `a` is singular.
pub fn partial_piv_lu(a: &DMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>, &'static str> {
    let decomp = a.clone().lu();

    let x = decomp.solve(b);
    let (p, l, u) = decomp.unpack();
    println!("╭───────────────────────────────");
    println!("│ LU Partial Pivoting");
    println!("├───────────────────────────────");
    println!("│ PA = LU");
    println!("├───────────────────────────────");
    println!("│ P: {:?}", p);
    println!("│ L: {}", l);
    println!("│ U: {}", u);
    println!("├───────────────────────────────");
    match &x {
        Some(x) => {
            println!("│ Ax = b => LUx = b");
            println!("├───────────────────────────────");
            println!("│ x: {}", x);
        }
        None => {
            println!("│ The matrix is singular");
        }
    };
    println!("╰───────────────────────────────");

    match x {
        Some(x) => Ok(x),
        None => Err("The matrix is singular"),
    }
}

/// Unpacks the LU decomposition into the lower and upper triangular matrices.
///
/// # Arguments
///
/// - `lu`: The LU decomposition of a matrix.
pub fn lu_unpack(lu: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let mut u = lu.upper_triangle();
    u.fill_diagonal(1.0);
    let l = lu.lower_triangle();
    (l, u)
}

/// Computes the LU decomposition of a matrix with full pivoting and returns
/// the solution of the system.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
/// - `b`: A vector of size `n`.
///
/// # Errors
///
/// - If matrix `a` is singular.
pub fn full_piv_lu(a: &DMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>, &'static str> {
    let decomp = a.clone().full_piv_lu();

    let x = decomp.solve(b);
    let (p, l, u, q) = decomp.unpack();
    println!("╭───────────────────────────────");
    println!("│ LU Full Pivoting");
    println!("├───────────────────────────────");
    println!("│ PAQ = LU");
    println!("├───────────────────────────────");
    println!("│ P: {:?}", p);
    println!("│ L: {}", l);
    println!("│ U: {}", u);
    println!("│ Q: {:?}", q);
    println!("├───────────────────────────────");
    match &x {
        Some(x) => {
            println!("│ Ax = b => LUx = b");
            println!("├───────────────────────────────");
            println!("│ x: {}", x);
        }
        None => {
            println!("│ The matrix is singular");
        }
    };
    println!("╰───────────────────────────────");

    match x {
        Some(x) => Ok(x),
        None => Err("The matrix is singular"),
    }
}

/// Solves a linear system using an LU decomposition.
///
/// # Arguments
///
/// - `l`: A lower triangular matrix of size `n x n`.
/// - `u`: An upper triangular matrix of size `n x n`.
/// - `b`: A vector of size `n`.
pub fn solve_lu(l: &DMatrix<f64>, u: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    // Step 1: Solve L * y = b
    let y = solve_triangular(l, b, false);
    // Step 2: Solve U * x = y
    let x = solve_triangular(u, &y, true);

    println!("╭───────────────────────────────");
    println!("│ LU");
    println!("│ Ax = b => LUx = b");
    println!("├───────────────────────────────");
    println!("│ L * y = b");
    println!("y: {}", y);
    println!("│ U * x = y");
    println!("x: {}", x);
    println!("╰───────────────────────────────");

    x
}

/// Creates an augmented matrix by combining a matrix and a vector.
///
/// # Arguments
///
/// - `a`: A matrix of size `n x n`.
/// - `b`: A vector of size `n`.
fn create_augmented(a: &DMatrix<f64>, b: &DVector<f64>) -> DMatrix<f64> {
    let mut augmented = DMatrix::<f64>::zeros(a.nrows(), a.ncols() + b.ncols());
    for i in 0..a.nrows() {
        // Copy the original matrix into the augmented matrix
        for j in 0..a.ncols() {
            augmented[(i, j)] = a[(i, j)];
        }
        // Copy the vector b into the last column of the augmented matrix
        augmented[(i, a.ncols())] = b[i];
    }
    augmented
}

/// Computes a triangular matrix with 1's on the diagonal by applying row operations,
/// including handling the augmented matrix (matrix + b), and returns the matrix and vector separately.
///
/// # Arguments
///
/// - `matrix`: A square matrix of size `n x n` (coefficients of the system).
/// - `b`: A vector of size `n` (constants of the system).
/// - `upper`: A boolean indicating if the matrix should be upper triangular.
///
/// # Returns
///
/// A tuple containing:
/// - The triangular matrix with 1's on the diagonal.
/// - The updated vector `b` after the elimination.
///
/// # Errors
///
/// - If a zero pivot is encountered, meaning the matrix is singular.
pub fn make_triangular(
    matrix: &DMatrix<f64>,
    b: &DVector<f64>,
    upper: bool,
) -> Result<(DMatrix<f64>, DVector<f64>), &'static str> {
    let mut augmented_matrix = create_augmented(matrix, b);

    // Perform Gaussian elimination to create a triangular matrix
    let num_rows = matrix.nrows();
    let row_iterator = match upper {
        true => (0..num_rows).collect::<Vec<_>>().into_iter(),
        false => (0..num_rows).rev().collect::<Vec<_>>().into_iter(),
    };

    for pivot in row_iterator {
        let pivot_value = augmented_matrix[(pivot, pivot)];

        // Ensure the pivot element is non-zero
        if pivot_value == 0.0 {
            return Err("Zero pivot encountered, matrix is singular");
        }

        // Normalize the pivot row to ensure the diagonal element is 1
        for col in 0..=num_rows {
            augmented_matrix[(pivot, col)] /= pivot_value;
        }

        // Use the pivot element to eliminate the other elements in the column
        let iterator = if upper { pivot + 1..num_rows } else { 0..pivot };
        for row in iterator {
            // Get number to eliminate
            // This is the element in the pivot column of the current row
            let factor = augmented_matrix[(row, pivot)];

            for col in 0..=num_rows {
                augmented_matrix[(row, col)] -= factor * augmented_matrix[(pivot, col)];
            }
        }
    }

    // Return the lower triangular matrix and the updated vector b
    Ok((
        augmented_matrix.columns(0, num_rows).into(),
        augmented_matrix.column(num_rows).into(),
    ))
}

/// Solves a linear system using the Gauss method.
///
/// # Complexity
///
/// - O((n^2)/2) for backward substitution
/// - O((n^3)/3) for forward substitution
///
/// # Arguments
///
/// - `matrix`: A square matrix of size `n x n`.
/// - `b`: A vector of size `n`.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn solve_gauss(
    matrix: &DMatrix<f64>,
    b: &DVector<f64>,
    backward_substitution: bool,
) -> Result<DVector<f64>, &'static str> {
    let (a, b) = make_triangular(matrix, b, backward_substitution)?;
    let x = solve_triangular(&a, &b, backward_substitution);

    println!("╭───────────────────────────────");
    println!("│ Gauss");
    println!("├───────────────────────────────");
    println!("A: {}", a);
    println!("b: {}", b);
    println!("X: {}", x);
    println!("╰───────────────────────────────");

    Ok(x)
}

/// Computes the solution of a linear system that minimizes the error on the solution.
///
/// An overdetermined system with the form `Ax = b` may not have a solution.
/// We find the `x` that minimizes `r = Ax - b`
///
/// # Arguments
///
/// - `a`: The A matrix.
/// - `b`: The b vector.
///
/// # Errors
///
/// - If matrix `a` is singular.
pub fn least_squares(a: &DMatrix<f64>, b: &DVector<f64>) -> Result<DVector<f64>, &'static str> {
    let a_t = a.transpose();
    let a_t_a = &a_t * a;
    let a_t_b = &a_t * b;

    println!("╭───────────────────────────────");
    println!("│ Least Squares");
    println!("├───────────────────────────────");
    println!("A^T * A: {}", a_t_a);
    println!("A^T * b: {}", a_t_b);
    println!("╰───────────────────────────────");

    let x = partial_piv_lu(&a_t_a, &a_t_b)?;

    Ok(x)
}

/// Solve a system of linear equations using the iterative method.
///
/// It is equivalent to the fixed-point method.
///
/// # Arguments
///
/// - `a`: The A matrix.
/// - `a_inv`: The inverse of the A matrix.
/// - `b`: The b vector.
/// - `iterations`: The number of iterations to perform.
pub fn iterative_solve(
    a: &DMatrix<f64>,
    a_inv: &DMatrix<f64>,
    b: &DVector<f64>,
    iterations: usize,
) -> DVector<f64> {
    // Compute the initial x
    let mut x = a_inv * b;

    // Iterate
    (0..iterations).for_each(|_| {
        let r = b - a * &x;
        let e = a_inv * r;
        x = &x + e;
    });

    x
}

/// Decomposes a matrix A such that A = L + D + U
///
/// # Arguments
/// - `matrix`: The matrix to decompose.
/// - `d_inverse`: A boolean indicating if the inverse of the diagonal matrix D should be computed.
///
/// # Returns
/// - A tuple containing the lower triangular matrix L, the diagonal matrix D or its inverse, and the upper triangular matrix U.
pub fn ldu_decomposition(
    matrix: &DMatrix<f64>,
    d_inverse: bool,
) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let mut l = matrix.lower_triangle();
    l.fill_diagonal(0.0);

    let mut u = matrix.upper_triangle();
    u.fill_diagonal(0.0);

    let d = match d_inverse {
        true => DMatrix::from_diagonal(&matrix.diagonal().map(|x| x.recip())),
        false => DMatrix::from_diagonal(&matrix.diagonal()),
    };

    (l, d, u)
}

/// Perform row permutations to make the diagonal elements the largest in each row.
///
/// # Arguments
///
/// - `a`: The matrix to permute.
///
/// # Returns
///
/// - A tuple containing the permuted matrix and the permutation matrix.
pub fn maximize_diagonal(a: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = a.nrows();
    let mut perm_a = a.clone();
    let mut perm_matrix = DMatrix::identity(n, n); // Identity matrix for tracking row swaps

    // Iterate over each row
    for i in 0..n {
        // Find the index of the largest element in the i-th row
        let max_index = (i..n)
            .max_by(|&x, &y| {
                perm_a[(i, x)]
                    .abs()
                    .partial_cmp(&perm_a[(i, y)].abs())
                    .unwrap()
            })
            .unwrap();

        if max_index != i {
            // Swap the rows in perm_a and perm_b
            perm_a.swap_rows(i, max_index);

            // Swap the corresponding rows in the permutation matrix
            perm_matrix.swap_rows(i, max_index);
        }
    }

    (perm_a, perm_matrix)
}

/// Creates and formats a row for the table.
///
/// # Arguments
///
/// - `i`: The index of the row.
/// - `v`: The vector to format.
/// - `e`: The error to format.
fn create_table_row(i: usize, v: &DVector<f64>, e: f64) -> Vec<String> {
    let n = v.len();
    let mut row = Vec::with_capacity(n + 2);
    row.push(format!("{}", i));
    (0..n).for_each(|i| row.push(format!("{:.6}", v[i])));
    row.push(format!("{:.3e}", e));
    row
}

/// Creates and formats a header for the table.
///
/// # Arguments
///
/// - `i`: The title of the first index column.
/// - `n`: The number of elements in the vector.
/// - `e`: The title of the error column.
fn create_table_header(i: &str, n: usize, e: &str) -> Vec<String> {
    let mut header = Vec::with_capacity(n + 2);
    header.push(i.to_owned());
    (0..n).for_each(|i| header.push(format!("x{}", i + 1)));
    header.push(e.to_owned());
    header
}

/// Checks if the matrix is diagonally dominant.
/// This means that the absolute value of the diagonal element is greater than
/// the sum of the absolute values of the other elements in the row.
///
/// # Arguments
///
/// - `a`: The matrix to check.
pub fn is_diagonally_dominant(a: &DMatrix<f64>) -> bool {
    let n = a.nrows();
    for i in 0..n {
        let sum = (0..n).filter(|&j| j != i).map(|j| a[(i, j)].abs()).sum();
        if a[(i, i)].abs() <= sum {
            return false;
        }
    }
    true
}

/// Performs the Jacobi algorithm to solve the system `Ax = b`.
///
/// # Arguments
///
/// - `a`: The A matrix.
/// - `b`: The b vector.
/// - `x_tolerance`: The tolerance for the x vector between iterations.
/// - `iterations_max`: The maximum number of iterations.
///
/// # Errors
///
/// - If the matrix is not diagonally dominant.
pub fn jacobi_solver(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    x_tolerance: f64,
    iterations_max: usize,
) -> Result<DVector<f64>, &'static str> {
    let n = a.nrows();
    let mut x_prev = DVector::zeros(n); // Initial guess (zero vector)
    let mut x = DVector::zeros(n);

    // Prepare the matrix for the Jacobi method
    let (a, _) = maximize_diagonal(a);
    let (l, d_inv, u) = ldu_decomposition(&a, true);

    if !is_diagonally_dominant(&a) {
        return Err("The matrix is not diagonally dominant");
    }

    // Prepare the table and start the iterations
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    rows.push(create_table_row(0, &x, 0.0));
    for iter in 0..iterations_max {
        // Jacobi update formula
        x = -&d_inv * (&l + &u) * &x_prev + &d_inv * b;
        let x_error = (&x - &x_prev).norm();

        // Store the information for the table
        rows.push(create_table_row(iter + 1, &x, x_error));

        // Check for convergence
        if x_error < x_tolerance {
            break;
        }

        x_prev = x.clone();
    }

    // Prepare the header and information
    let info = vec![
        "Running jacobi iterative solver".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let header = create_table_header("i", n, "x_error");

    // Format and print the table
    let table = table_formatter(info, header, rows).unwrap();
    println!("{}", table);

    Ok(x)
}

/// Performs the Gauss-Seidel algorithm to solve the system `Ax = b`.
///
/// The difference between the Gauss-Seidel and Jacobi methods is that the
/// Gauss-Seidel method uses the updated values of `x` in the same iteration.
///
/// # Arguments
///
/// - `a`: The A matrix.
/// - `b`: The b vector.
/// - `x_tolerance`: The tolerance for the x vector between iterations.
/// - `iterations_max`: The maximum number of iterations.
///
/// # Errors
///
/// - If the matrix is not diagonally dominant.
pub fn gauss_seidel_solver(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    x_tolerance: f64,
    iterations_max: usize,
) -> Result<DVector<f64>, &'static str> {
    let n = a.nrows();
    let mut x_prev = DVector::zeros(n); // Initial guess (zero vector)
    let mut x = DVector::zeros(n);

    // Prepare the matrix for the Jacobi method
    let (a, _) = maximize_diagonal(a);

    if !is_diagonally_dominant(&a) {
        return Err("The matrix is not diagonally dominant");
    }

    // Prepare the table and start the iterations
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    rows.push(create_table_row(0, &x, 0.0));
    for iter in 0..iterations_max {
        // Update each element of x
        for i in 0..n {
            let sum1: f64 = (0..i).map(|j| a[(i, j)] * x[j]).sum();
            let sum2: f64 = (i + 1..n).map(|j| a[(i, j)] * x_prev[j]).sum();

            // Gauss-Seidel update rule
            x[i] = (b[i] - sum1 - sum2) / a[(i, i)];
        }
        let x_error = (&x - &x_prev).norm();

        // Store the information for the table
        rows.push(create_table_row(iter + 1, &x, x_error));

        // Check for convergence
        if x_error < x_tolerance {
            break;
        }

        x_prev = x.clone();
    }

    // Prepare the header and information
    let info = vec![
        "Running Gauss-Seidel iterative solver".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let header = create_table_header("i", n, "x_error");

    // Format and print the table
    let table = table_formatter(info, header, rows).unwrap();
    println!("{}", table);

    Ok(x)
}

/// Performs the relaxation method to solve the system `Ax = b`.
///
/// Convergence usually only happens for `0<w<2`.
///
/// - If `w<1`, the method is called under-relaxation.
/// - If `w=1`, the method is equivalent to the Gauss-Seidel method.
/// - If `w>1`, the method is called over-relaxation.
///
/// # Arguments
///
/// - `a`: The A matrix.
/// - `b`: The b vector.
/// - `w`: The relaxation factor.
/// - `x_tolerance`: The tolerance for the x vector between iterations.
/// - `iterations_max`: The maximum number of iterations.
///
/// # Errors
///
/// - If the matrix is not diagonally dominant.
pub fn relaxation_solver(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    w: f64,
    x_tolerance: f64,
    iterations_max: usize,
) -> Result<DVector<f64>, &'static str> {
    let n = a.nrows();
    let mut x_prev = DVector::zeros(n); // Initial guess (zero vector)
    let mut x = DVector::zeros(n);

    // Prepare the matrix for the Jacobi method
    let (a, _) = maximize_diagonal(a);

    if !is_diagonally_dominant(&a) {
        return Err("The matrix is not diagonally dominant");
    }

    // Prepare the table and start the iterations
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max + 1);
    rows.push(create_table_row(0, &x, 0.0));
    for iter in 0..iterations_max {
        // Update each element of x
        for i in 0..n {
            let sum1: f64 = (0..i).map(|j| a[(i, j)] * x[j]).sum();
            let sum2: f64 = (i + 1..n).map(|j| a[(i, j)] * x_prev[j]).sum();

            // Relaxation (Successive Over-Relaxation) update rule
            let new_value = (b[i] - sum1 - sum2) / a[(i, i)];
            x[i] = (1.0 - w) * x_prev[i] + w * new_value;
        }
        let x_error = (&x - &x_prev).norm();

        // Store the information for the table
        rows.push(create_table_row(iter + 1, &x, x_error));

        // Check for convergence
        if x_error < x_tolerance {
            break;
        }

        x_prev = x.clone();
    }

    // Prepare the header and information
    let info = vec![
        "Running relaxation iterative solver".to_owned(),
        format!("w = {}", w),
        format!("x tolerance = {:e}", x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let header = create_table_header("i", n, "x_error");

    // Format and print the table
    let table = table_formatter(info, header, rows).unwrap();
    println!("{}", table);

    Ok(x)
}

/// Solves non-linear systems using the fixed-point method.
///
/// # Arguments
/// - `system`: A vector of functions representing the system of equations.
/// - `x_initial`: The initial guess for the solution.
/// - `x_tolerance`: The tolerance for the x vector between iterations.
/// - `f_x_tolerance`: The tolerance for the function value.
/// - `iterations_max`: The maximum number of iterations.
///
/// # Errors
///
/// - If the sum of the absolute values of the partial derivatives is greater than 1.
///     This is a necessary condition for the convergence of the fixed-point method.
pub fn fixed_point_system(
    system: &[FuncMulti],
    x_initial: &DVector<f64>,
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<DVector<f64>, &'static str> {
    let ncols = x_initial.len();

    // Initialize the variables
    let mut x_prev = x_initial.clone();
    let mut x = x_initial.clone();
    let mut fx = function_vec(system, x_prev.as_slice());

    if !gradient_mat(system, x_prev.as_slice())
        .abs()
        .row_sum()
        .iter()
        .all(|&x| x < 1.0)
    {
        return Err(
            "For each variable, the sum of the absolute values of the partial derivatives must be less than 1",
        );
    }

    // Prepare the table and start the iterations
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max + 1);
    let mut row = create_table_row(0, &x, 0.0);
    row.push(format!("{:.6e}", fx));
    rows.push(row);

    for iter in 0..iterations_max {
        // Fixed-point update formula
        x = function_vec(system, x_prev.as_slice());

        // Calculate the error
        let x_error = (&x - &x_prev).norm();

        // Update the function value and jacobian
        x_prev = x.clone();
        fx = function_vec(system, x_prev.as_slice());

        // Store the information for the table
        let mut row = create_table_row(iter + 1, &x, x_error);
        row.push(format!("{:.6e}", fx));
        rows.push(row);

        // Check for convergence
        if x_error < x_tolerance || fx.norm() < f_x_tolerance {
            break;
        }
    }

    // Prepare the header and information
    let info = vec![
        "Running fixed point iterative solver".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let mut header = create_table_header("i", ncols, "x_error");
    header.push("f(x)".to_owned());

    // Format and print the table
    let table = table_formatter(info, header, rows).unwrap();
    println!("{}", table);

    Ok(x)
}

/// Solves non-linear systems using the Newton method.
///
/// # Arguments
///
/// - `system`: A vector of functions representing the system of equations.
/// - `x_initial`: The initial guess for the solution.
/// - `x_tolerance`: The tolerance for the x vector between iterations.
/// - `f_x_tolerance`: The tolerance for the function value.
/// - `iterations_max`: The maximum number of iterations.
///
/// # Errors
///
/// - If the Jacobian matrix is singular.
pub fn newton_system(
    system: &[FuncMulti],
    x_initial: &DVector<f64>,
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<DVector<f64>, &'static str> {
    let ncols = x_initial.len();

    // Initialize the variables
    let mut x_prev = x_initial.clone();
    let mut x = x_initial.clone();
    let mut jacobian = gradient_mat(system, x_prev.as_slice());
    let mut fx = function_vec(system, x_prev.as_slice());

    // Prepare the table and start the iterations
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max + 1);
    let mut row = create_table_row(0, &x, 0.0);
    row.push(format!("{:.6e}", fx));
    row.push(format!("{:.6}", jacobian));
    rows.push(row);

    for iter in 0..iterations_max {
        // Newton update formula
        let jacobian_inv = match jacobian.clone().try_inverse() {
            Some(j) => j,
            None => return Err("The jacobi matrix is singular"),
        };
        x = &x_prev - jacobian_inv * &fx;

        // Calculate the error
        let x_error = (&x - &x_prev).norm();

        // Update the function value and jacobian
        x_prev = x.clone();
        fx = function_vec(system, x_prev.as_slice());
        jacobian = gradient_mat(system, x_prev.as_slice());

        // Store the information for the table
        let mut row = create_table_row(iter + 1, &x, x_error);
        row.push(format!("{:.6e}", fx));
        row.push(format!("{:.6}", jacobian));
        rows.push(row);

        // Check for convergence
        if x_error < x_tolerance || fx.norm() < f_x_tolerance {
            break;
        }
    }

    // Prepare the header and information
    let info = vec![
        "Running Newton iterative solver".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let mut header = create_table_header("i", ncols, "x_error");
    header.push("f(x)".to_owned());
    header.push("J(x)".to_owned());

    // Format and print the table
    let table = table_formatter(info, header, rows).unwrap();
    println!("{}", table);

    Ok(x)
}
