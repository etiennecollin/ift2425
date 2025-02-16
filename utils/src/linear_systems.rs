use nalgebra::{DMatrix, DVector};

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
pub fn solve_cramer(matrix: &DMatrix<f64>, vector: &DVector<f64>) -> DVector<f64> {
    let mut result = DVector::zeros(vector.len());
    let det = matrix.determinant();
    println!("╭───────────────────────────────");
    println!("│ Cramer");
    println!("├───────────────────────────────");
    println!("│ det(A) = {det}");

    for i in 0..vector.len() {
        let mut matrix_aux = matrix.clone();
        matrix_aux.set_column(i, vector);
        let det_aux = matrix_aux.determinant();
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
    result
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
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
///
/// # Returns
///
/// - A tuple containing the lower triangular matrix `L` and the upper triangular matrix `U`.
// pub fn no_piv_lu(matrix: &DMatrix<f64>) -> Result<(DMatrix<f64>, DMatrix<f64>), &'static str> {
//     let n = matrix.nrows();
//     let mut l = DMatrix::<f64>::zeros(n, n);
//     let mut u = matrix.clone();
//
//     for i in 0..n {
//         // L diagonal elements are all 1
//         l[(i, i)] = 1.0;
//
//         // Upper Triangular: Make elements below the diagonal zero
//         for j in i..n {
//             u[(i, j)] -= (0..i).map(|k| l[(i, k)] * u[(k, j)]).sum::<f64>();
//         }
//
//         // Lower Triangular: Make elements above the diagonal zero
//         for j in i + 1..n {
//             // Check if the diagonal element is zero
//             // If it is, the matrix is singular
//             if u[(i, i)] == 0.0 {
//                 return Err("Matrix is singular");
//             }
//
//             l[(j, i)] =
//                 (u[(j, i)] - (0..i).map(|k| l[(j, k)] * u[(k, i)]).sum::<f64>()) / u[(i, i)];
//
//             // Zero out the elements in the upper triangle
//             u[(j, i)] = 0.0;
//         }
//     }
//
//     Ok((l, u))
// }

/// Computes the LU decomposition of a matrix without pivoting.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
///
/// # Returns
///
/// - A tuple containing the lower triangular matrix `L` and the upper triangular matrix `U`.
pub fn no_piv_lu(matrix: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = matrix.nrows();
    assert_eq!(n, matrix.ncols(), "The matrix must be square");

    let mut l = DMatrix::zeros(n, n);
    let mut u = DMatrix::identity(n, n);

    // Initialize the first column of L
    for i in 0..n {
        l[(i, 0)] = matrix[(i, 0)];
    }

    // Initialize the first row of U
    for j in 1..n {
        u[(0, j)] = matrix[(0, j)] / l[(0, 0)];
    }

    // Compute elements of L and U
    for j in 1..n {
        // Compute elements of L
        for i in j..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + l[(i, k)] * u[(k, j)]);
            l[(i, j)] = matrix[(i, j)] - sum;
        }

        // Compute elements of U
        for i in (j + 1)..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + l[(j, k)] * u[(k, i)]);
            u[(j, i)] = (matrix[(j, i)] - sum) / l[(j, j)];
        }
    }

    (l, u)
}

/// Computes the LU decomposition of a matrix without pivoting. The input matrix is modified in-place.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
pub fn no_piv_lu_mut(matrix: &mut DMatrix<f64>) {
    let n = matrix.nrows();
    assert_eq!(n, matrix.ncols(), "The matrix must be square");

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
}

/// Computes the LU decomposition of a matrix without pivoting.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
pub fn no_piv_lu_packed(matrix: &DMatrix<f64>) -> DMatrix<f64> {
    let n = matrix.nrows();
    assert_eq!(n, matrix.ncols(), "The matrix must be square");

    let mut lu = DMatrix::zeros(n, n);
    lu.copy_from(matrix);

    // Initialize the first row of U
    for j in 1..n {
        lu[(0, j)] /= lu[(0, 0)];
    }

    // Compute elements of L and U
    for j in 1..n {
        // Compute elements of L
        for i in j..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + lu[(i, k)] * lu[(k, j)]);
            lu[(i, j)] -= sum;
        }

        // Compute elements of U
        for i in (j + 1)..n {
            let sum = (0..j).fold(0.0, |sum, k| sum + lu[(j, k)] * lu[(k, i)]);
            lu[(j, i)] = (lu[(j, i)] - sum) / lu[(j, j)];
        }
    }

    lu
}

/// Computes the LU decomposition of a matrix with partial pivoting.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
/// - `b`: A vector of size `n`.
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

/// Computes the LU decomposition of a matrix with full pivoting.
///
/// # Arguments
///
/// - `a`: A square matrix of size `n x n`.
/// - `b`: A vector of size `n`.
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
pub fn solve_gauss(
    matrix: &DMatrix<f64>,
    b: &DVector<f64>,
    backward_substitution: bool,
) -> DVector<f64> {
    let (a, b) = make_triangular(matrix, b, backward_substitution).unwrap();
    let x = solve_triangular(&a, &b, backward_substitution);

    println!("╭───────────────────────────────");
    println!("│ Gauss");
    println!("├───────────────────────────────");
    println!("A: {}", a);
    println!("b: {}", b);
    println!("X: {}", x);
    println!("╰───────────────────────────────");

    x
}
