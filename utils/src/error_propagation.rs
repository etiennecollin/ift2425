use crate::utils::{FuncMulti, compute_gradient};

/// Computes the position of the last significant digit in the error.
///
/// Gets the position `r` of the digit in `error` such that `error < 0.5 * 10^r`.
///
/// # Arguments
///
/// - `error`: The error on the function.
pub fn compute_sig_position(error: f64) -> i32 {
    (error * 2.0).log10().ceil() as i32
}

/// Round a number to the given number of digits.
///
/// # Arguments
///
/// - `x`: The number to round.
/// - `digits`: The number of digits to round to.
fn round_to_digits(x: f64, digits: i32) -> f64 {
    (x * 10f64.powi(digits)).round() / 10f64.powi(digits)
}

/// Correctly rounds value of a function at the given point with proper error.
///
/// # Arguments
///
/// - `f`: A function that takes a slice of f64 and returns a f64. Each element of the slice represents a variable in the function.
/// - `point`: A slice of f64 representing the point where the function is evaluated.
/// - `error`: A slice of f64 representing the error on each variable of the function.
pub fn cse_taylor(f: FuncMulti, point: &[f64], error: &[f64]) -> f64 {
    // Compute the function on the point
    let res = f(point);

    // Compute the gradient of the function at the point
    let gradient = compute_gradient(f, point);

    // Get the error on f
    let df: f64 = gradient
        .iter()
        .enumerate()
        .fold(0.0, |acc, (i, x)| acc + x.abs() * error[i]);

    // Get the position r of the digit in df such that df < 0.5 * 10^r
    let r = compute_sig_position(df);

    // Round to r decimal places
    let res_rounded = round_to_digits(res, -r);
    let df_rounded = round_to_digits(df, -r + 2);

    // Print results
    println!("╭───────────────");
    println!("│ f*(point) = {}", res);
    println!("│ ∇f* = {:?}", gradient);
    println!("├─");
    println!("│ Using the first terms of Taylor's formula:");
    println!("│ Δf = {}", df);
    println!("│ Δf/f = {:.3}%", df / res * 100.0);
    println!("├─");
    println!("│ Δf < 0.5 * 10^r -> r = {}", r);
    println!(
        "│ f(point) = {} ± {} or {:.3}%",
        res_rounded,
        df_rounded,
        df / res * 100.0
    );
    println!("╰───────────────");

    res_rounded
}

/// Correctly rounds value of a function at the given point with proper error.
///
/// This uses the fork method to compute the error on the function.
/// This method is only valid if the function is continuous and monotonic.
///
/// # Arguments
///
/// - `f`: A function that takes a slice of f64 and returns a f64. Each element of the slice represents a variable in the function.
/// - `point`: A slice of f64 representing the point where the function is evaluated.
/// - `error`: A slice of f64 representing the error on each variable of the function.
pub fn cse_fork(f: FuncMulti, point: &[f64], error: &[f64]) -> f64 {
    // Number of variables
    let n = point.len();

    // Iterate through all combinations of +/- error for each variable
    // For each variable, we can either add or subtract the error
    // That means there are 2^n possible combinations
    let mut all_points = Vec::new();
    for i in 0..(1 << n) {
        let mut new_point = Vec::with_capacity(n);

        for j in 0..n {
            let sign = if (i >> j) & 1 == 0 { -1.0 } else { 1.0 };
            new_point.push(point[j] + sign * error[j]);
        }

        all_points.push(new_point);
    }

    // Evaluate the function at each of these points
    let mut function_values = Vec::new();
    for p in all_points {
        function_values.push(f(&p));
    }

    // Get the minimum and maximum values of the function
    let min = function_values
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let max = function_values
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);

    // Get the error on f
    let df = (max - min) / 2.0;

    // Get the position r of the digit in df such that df < 0.5 * 10^r
    let r = compute_sig_position(df);

    // Round to r decimal places
    let res = f(point);
    let res_rounded = round_to_digits(res, -r);
    let df_rounded = round_to_digits(df, -r + 2);

    println!("╭───────────────");
    println!("│ f*(point) = {}", res);
    println!("├─");
    println!("│ Using the fork method, assuming f is continuous and monotonic:");
    println!("│ f_min = {}", min);
    println!("│ f_max = {}", max);
    println!("│ Δf = {}", df);
    println!("│ Δf/f = {:.3}%", df / res * 100.0);
    println!("├─");
    println!("│ Δf < 0.5 * 10^r -> r = {}", r);
    println!(
        "│ f(point) = {} ± {} or {:.3}%",
        res_rounded,
        df_rounded,
        df / res * 100.0
    );
    println!("╰───────────────");

    res_rounded
}
