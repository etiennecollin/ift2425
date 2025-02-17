use nalgebra::{DMatrix, DVector};

/// Compute the p-norm of a vector.
///
/// # Arguments
///
/// - `v`: A vector.
/// - `p`: The norm to compute.
pub fn norm_vec(v: &DVector<f64>, p: u32) -> f64 {
    if p == 0 {
        println!("Computing the norm with p=infinity");
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
/// - `p`: The norm to compute.
pub fn norm_mat(mat: &DMatrix<f64>, p: u32) -> f64 {
    if p == 0 {
        println!("Computing the norm with p=infinity");
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
/// - `p`: The norm to compute.
pub fn condition_number(mat: &DMatrix<f64>, p: u32) -> Result<f64, &'static str> {
    let mat_inv = match mat.clone().try_inverse() {
        Some(mat_inv) => mat_inv,
        None => return Err("Matrix is singular"),
    };
    println!("Inverse: {}", mat_inv);
    Ok(norm_mat(mat, p) * norm_mat(&mat_inv, p))
}
