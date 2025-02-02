use crate::utils::{compute_gradient, FuncMulti};

/// This function the correctly rounded value of the function at the given point with proper error.
pub fn cse(f: FuncMulti, point: &[f64], error: &[f64], epsilon: f64) -> f64 {
    // Compute the function on the point
    let res = f(point);

    // Compute the gradient of the function at the point
    let gradient = compute_gradient(f, point, epsilon);

    // Get the error on f
    let df: f64 = gradient
        .iter()
        .enumerate()
        .fold(0.0, |acc, (i, x)| acc + x.abs() * error[i]);

    // Get the position r of the digit in df such that df < 0.5 * 10^r
    let r: i32 = -(df * 2.0).log10().ceil() as i32;

    // Round to e decimal places
    let res_rounded = (res * 10f64.powi(r)).round() / 10f64.powi(r);

    // Print results
    println!("╭───────────────");
    println!("│ f*(point) = {}", res);
    println!("│ ∇f* = {:?}", gradient);
    println!("├─");
    println!("│ Using the first terms of Taylor's formula:");
    println!("│ Δf = {}", df);
    println!("│ Δf/f = {:.3}%", df / res * 100.0);
    println!("├─");
    println!("│ r = {}", -r);
    println!("│ f(point) = {}", res_rounded);
    println!("╰───────────────");

    res_rounded
}
