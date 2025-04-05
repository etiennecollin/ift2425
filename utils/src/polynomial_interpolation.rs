use core::f64;

use crate::utils::{FuncSingle, choose_float, factorial};

/// Computes the Lagrange interpolation polynomial at a given point `x`.
///
/// The degree of the polynomial is `n`, where `n` is the number of data points minus one.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `x`: The point at which to evaluate the polynomial.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The y-coordinates of the data points.
pub fn lagrange(degree: usize, x: f64, xs: &[f64], fs: &[f64]) -> Result<f64, &'static str> {
    if xs.len() < 2 {
        return Err("The length of xs must be at least 2.");
    }
    if xs.len() != fs.len() {
        return Err("The length of xs and fs must be equal.");
    }
    if !xs.iter().is_sorted() {
        return Err("The x-coordinates must be sorted.");
    }
    if degree < 1 {
        return Err("The degree must be at least 1.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let mut sum = 0.0;
    for i in 0..=degree {
        sum += polynomial(degree, i, x, xs)? / polynomial(degree, i, xs[i], xs)? * fs[i]
    }

    println!("╭───────────────");
    println!("│ Lagrange Polynomial Interpolation");
    println!("├─");
    println!("│ P_{}({}) = {}", degree, x, sum);
    println!("╰───────────────");
    Ok(sum)
}

/// Computes the Lagrange basis polynomial for a given index `i`.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `i`: The index of the basis polynomial.
/// - `x`: The point at which to evaluate the polynomial.
/// - `xs`: The x-coordinates of the data points.
fn polynomial(degree: usize, i: usize, x: f64, xs: &[f64]) -> Result<f64, &'static str> {
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let mut product = 1.0;
    for (j, xj) in xs.iter().take(degree + 1).enumerate() {
        if i == j {
            continue;
        }
        product *= x - xj;
    }

    Ok(product)
}

/// Computes the error range of a Lagrange polynomial interpolation at a given point `x`.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `derivative`: The `n+1`th derivative of the function approximated by the polynomial of
///   degree `n`
/// - `x`: The point at which to evaluate the polynomial.
/// - `xs`: The x-coordinates of the data points.
///
/// # Errors
///
/// - If the length of `xs` is less than 2, an error is returned.
/// - If the length of `xs` is greater than 170, an error is returned.
pub fn lagrange_polynomial_error_range(
    degree: usize,
    derivative: FuncSingle,
    x: f64,
    xs: &[f64],
) -> Result<(f64, f64), &'static str> {
    if xs.len() < 2 {
        return Err("The length of xs must be at least 2.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;

    // Get min and max of xs
    xs.iter().take(degree + 1).for_each(|&xi| {
        if xi < x_min {
            x_min = xi;
        }
        if xi > x_max {
            x_max = xi;
        }
    });

    // Pre-compute the factorial
    let factorial = factorial((degree + 1) as u8)?;

    // Pre-compute the product
    let mut product = 1.0;
    for xi in xs.iter().take(degree + 1) {
        product *= x - xi;
    }

    // Compute the error bounds
    let factor = product / factorial;
    let mut bound_min = derivative(x_min) * factor;
    let mut bound_max = derivative(x_max) * factor;

    // Get the min and max bounds
    if bound_min > bound_max {
        std::mem::swap(&mut bound_min, &mut bound_max);
    }

    println!("╭───────────────");
    println!("│ Polynomial Error Range");
    println!("├─");
    println!("│ {:.6e} <= E({}) <= {:.6e}", bound_min, x, bound_max);
    println!("╰───────────────");

    Ok((bound_min, bound_max))
}

/// Computes the Forward Newton-Gregory interpolation polynomial at a given point `x`.
///
/// This method assumes all x-coordinates are equally spaced.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `x`: The point at which to evaluate the polynomial.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The y-coordinates of the data points.
pub fn newton_gregory_forward(
    degree: usize,
    x: f64,
    xs: &[f64],
    fs: &[f64],
) -> Result<f64, &'static str> {
    if xs.len() < 2 {
        return Err("The length of xs must be at least 2.");
    }
    if xs.len() != fs.len() {
        return Err("The length of xs and fs must be equal.");
    }
    if !xs.iter().is_sorted() {
        return Err("The x-coordinates must be sorted.");
    }
    if degree < 1 {
        return Err("The degree must be at least 1.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    // Pre-compute the finite differences
    let h = xs[1] - xs[0];
    let s = (x - xs[0]) / h;

    // Newton-Gregory Forward Polynomial Interpolation
    let mut result = fs[0];
    for k in 1..=degree {
        let factor = delta_f_i(k, 0, fs) / factorial(k as u8)?;
        let mut product = 1.0;
        for j in 0..k {
            product *= s - j as f64
        }
        result += factor * product
    }

    println!("╭───────────────");
    println!("│ Newton-Gregory Forward Polynomial Interpolation");
    println!("├─");
    println!("│ Assuming equally spaced x-coordinates");
    println!("├─");
    println!("│ P_{}({}) = {}", degree, x, result);
    println!("╰───────────────");

    Ok(result)
}

pub fn scratch(degree: usize, x: f64, xs: &[f64], fs: &[f64]) -> Result<f64, &'static str> {
    if xs.len() < 2 {
        return Err("The length of xs must be at least 2.");
    }
    if xs.len() != fs.len() {
        return Err("The length of xs and fs must be equal.");
    }
    if !xs.iter().is_sorted() {
        return Err("The x-coordinates must be sorted.");
    }
    if degree < 1 {
        return Err("The degree must be at least 1.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    // Pre-compute the finite differences
    let h = xs[1] - xs[0];
    let s = (x - xs[0]) / h;

    // Newton-Gregory Forward Polynomial Interpolation
    let mut result = fs[0];
    for n in 1..=degree {
        result += choose_float(s, n as u64)? * delta_f_i(n, 0, fs);
    }

    // FIXME: This does not give the same result as newton_gregory_forward
    // and it shouvlve be the same
    println!("╭───────────────");
    println!("│ Newton-Gregory Forward Polynomial Interpolation");
    println!("├─");
    println!("│ Assuming equally spaced x-coordinates");
    println!("│ P_(n+1)(x) = P_n(x) + `s choose n+1` * delta^(n+1)f_0");
    println!("├─");
    println!("│ P_{}({}) = {}", degree, x, result);
    println!("╰───────────────");

    Ok(result)
}

pub fn newton_gregory_forward_error() -> Result<f64, &'static str> {
    todo!("Ch.4 p.21")
}

/// Computes the `delta`th finite difference of the `i`th element in the array `fs`.
///
/// # Arguments
///
/// - `delta`: The order of the finite difference.
/// - `i`: The index of the element in the array.
/// - `fs`: The array of function values.
pub fn delta_f_i(delta: usize, i: usize, fs: &[f64]) -> f64 {
    let mut dp = vec![vec![None; fs.len()]; delta + 1];
    _delta_f_i(delta, i, fs, &mut dp)
}

/// Recursively computes the `delta`th finite difference of the `i`th element in the array `fs`.
/// This function uses memoization to store previously computed results in the `dp` table.
///
/// # Arguments
///
/// - `delta`: The order of the finite difference.
/// - `i`: The index of the element in the array.
/// - `fs`: The array of function values.
/// - `dp`: The memoization table.
fn _delta_f_i(delta: usize, i: usize, fs: &[f64], dp: &mut [Vec<Option<f64>>]) -> f64 {
    // Check if the value is already computed
    if let Some(result) = dp[delta][i] {
        return result;
    }

    // Base cases
    let result = match delta {
        0 => fs[i],
        1 => fs[i + 1] - fs[i],
        _ => _delta_f_i(delta - 1, i + 1, fs, dp) - _delta_f_i(delta - 1, i, fs, dp),
    };

    // Store the result in dp table for future reuse
    dp[delta][i] = Some(result);

    result
}
