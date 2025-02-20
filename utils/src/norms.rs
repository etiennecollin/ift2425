use nalgebra::{DMatrix, DVector};

/// Compute the p-norm of a vector.
///
/// # Arguments
///
/// - `v`: A vector.
/// - `p`: The norm to compute. A value of 0 corresponds to the infinity norm.
pub fn norm_vec(v: &DVector<f64>, p: u32) -> f64 {
    if p == 0 {
        v.abs().max()
    } else if p == 1 {
        v.abs().sum()
    } else if p == 2 {
        v.norm()
    } else {
        v.abs().map(|x| x.powi(p as i32)).sum().powf(1.0 / p as f64)
    }
}

/// Compute the p-norm of a matrix.
///
/// # Arguments
///
/// - `mat`: A matrix.
/// - `p`: The norm to compute. A value of 0 corresponds to the infinity norm.
pub fn norm_mat(mat: &DMatrix<f64>, p: u32) -> f64 {
    if p == 0 {
        mat.abs().column_sum().max()
    } else if p == 1 {
        mat.abs().row_sum().max()
    } else if p == 2 {
        // See definition of the spectral norm of a matrix
        // https://en.wikipedia.org/wiki/Matrix_norm#Spectral_norm_(p_=_2)
        mat.singular_values_unordered().max()
    } else {
        println!("Computing the euclidean norm");
        // mat.abs().map(|x| x.powi(2)).sum().sqrt()
        mat.norm()
    }
}

/// Compute the condition number of a matrix.
///
/// A matrix is well conditioned if its condition number is close to 1.
///
/// # Arguments
///
/// - `mat`: A matrix.
/// - `p`: The norm to compute. A value of 0 corresponds to the infinity norm.
///
/// # Errors
///
/// - If the matrix is singular.
pub fn condition_number(mat: &DMatrix<f64>, p: u32) -> Result<f64, &'static str> {
    let mat_inv = match mat.clone().try_inverse() {
        Some(mat_inv) => mat_inv,
        None => return Err("Matrix is singular"),
    };
    println!("Inverse: {}", mat_inv);
    Ok(norm_mat(mat, p) * norm_mat(&mat_inv, p))
}

/// Computes the range of the error on x given that the error on b is known.
///
/// # Arguments
///
/// - `a`: A matrix.
/// - `b`: The true b vector.
/// - `b_star`: The approximated b vector.
/// - `p`: The norm to use. Use the 1 norm or the 0 (infinite) norm
///
/// # Errors
///
/// - If the matrix is singular.
/// - If the vector norm is not compatible with the matrix norm.
pub fn x_error_range_b(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    b_star: &DVector<f64>,
    p: u32,
) -> Result<(f64, f64), &'static str> {
    if p != 0 && p != 1 {
        return Err("The vector and matrix norms must be compatible. Use the 1 norm or the 0 (infinite) norm.");
    }

    let db = b_star - b;
    let cond = condition_number(a, p)?;
    let frac = norm_vec(&db, p) / norm_vec(b, p);
    let range = (cond.recip() * frac, cond * frac);

    println!("Condition number: {}", cond);
    println!("Fraction: {}", frac);
    println!("Error range: {:?}", range);

    Ok(range)
}

/// Computes the upper bound of the error on x given that the error on A is known.
///
/// # Arguments
///
/// - `a`: The true A matrix.
/// - `a_star`: The approximated A matrix.
/// - `p`: The norm to use. Use the 1 norm or the 0 (infinite) norm
///
/// # Errors
///
/// - If the matrix is singular.
pub fn x_error_bound_a(
    a: &DMatrix<f64>,
    a_star: &DMatrix<f64>,
    p: u32,
) -> Result<f64, &'static str> {
    if p != 0 && p != 1 {
        return Err("The vector and matrix norms must be compatible. Use the 1 norm or the 0 (infinite) norm.");
    }

    let da = a_star - a;
    let cond = condition_number(a, p)?;
    let frac = norm_mat(&da, p) / norm_mat(a, p);
    let bound = cond * frac;

    println!("Condition number: {}", cond);
    println!("Fraction: {}", frac);
    println!("Error upper bound: {}", bound);

    Ok(bound)
}
