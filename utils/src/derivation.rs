use crate::{
    interpolation::finite_diff_table,
    utils::{FuncSingle, choose_float},
};

/// Computes the derivative of the Forward Newton-Gregory interpolation polynomial at a given point `x`.
///
/// This method assumes all x-coordinates are equally spaced.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `x_index`: The index of the point in the x-coordinates array.
/// - `h`: The spacing between the x-coordinates.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The y-coordinates of the data points.
pub fn newton_gregory_forward_derivative(
    degree: usize,
    x_index: usize,
    h: f64,
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
    if x_index > xs.len() - 1 {
        return Err("The x_index must be less than the number of data points.");
    }

    // Pre-compute the finite differences
    let finite_diff = finite_diff_table(fs.len() - 1, fs)?;

    println!("╭───────────────");
    println!("│ Newton-Gregory Forward Polynomial Interpolation Derivative");
    println!("├─");
    println!("│ Assuming equally spaced x-coordinates");
    println!("│ Assuming s = 0");
    println!("│ h = {}", h);
    println!("├─");
    println!("│ sum = 0.0");

    // Newton-Gregory Forward Polynomial Interpolation Derivative
    let mut sum = 0.0;
    for k in 1..=degree {
        let sign = if k % 2 == 0 { -1.0 } else { 1.0 };
        sum += sign * finite_diff[x_index][k] / k as f64;
        println!(
            "│ sum += {sign}/{} * {:.4e} = {:.6e}",
            k, finite_diff[x_index][k], sum
        );
    }

    println!("│ sum = {:.6e}", sum);
    let result = sum / h;

    println!("├─");
    println!("│ P'_{}({}) = sum/h = {:.6e}", degree, xs[x_index], result);
    println!("╰───────────────");

    Ok(result)
}

/// FIXME: Does not work
///
/// Computes the error estimate for the Newton-Gregory derivative at the base point `x0`
/// using the finite difference table.
///
/// # Arguments
///
/// - `degree`: Degree `n` of the Newton-Gregory polynomial.
/// - `x`: The point at which to estimate the derivative.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The function values.
/// - `h`: Step size.
/// - `interval_indices`: The indices in `xs` of the interval in which to estimate the error.
///
pub fn newton_gregory_derivative_error_estimate(
    degree: usize,
    x_index: usize,
    xs: &[f64],
    fs: &[f64],
    h: f64,
) -> Result<(f64, f64), &'static str> {
    if xs.len() < 2 {
        return Err("The length of xs must be at least 2.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let mut x_min: (usize, f64) = (0, f64::INFINITY);
    let mut x_max: (usize, f64) = (0, f64::NEG_INFINITY);

    // Get min and max of xs
    xs.iter().enumerate().for_each(|(i, &xi)| {
        if xi < x_min.1 {
            x_min = (i, xi);
        }
        if xi > x_max.1 {
            x_max = (i, xi);
        }
    });

    let finite_differences = finite_diff_table(fs.len() - 1, fs)?;
    let sign = if degree % 2 == 0 { 1.0 } else { -1.0 };
    let factor = sign / ((degree + 1) as f64) * h.powi(degree as i32);
    let order = degree + 1;

    let mut bound_min = finite_differences[x_min.0][order] / h.powi(order as i32) * factor;
    let mut bound_max = finite_differences[x_max.0][order] / h.powi(order as i32) * factor;

    // Get the min and max bounds
    if bound_min > bound_max {
        std::mem::swap(&mut bound_min, &mut bound_max);
    }

    println!("╭───────────────");
    println!("│ Newton-Gregory Derivative Error Estimate");
    println!("├─");
    println!("│ x={} <= E({}) <= x={}", x_min.1, xs[x_index], x_max.1);
    println!(
        "│ {:.6e} <= E({}) <= {:.6e}",
        bound_min, xs[x_index], bound_max
    );
    println!("╰───────────────");

    Ok((bound_min, bound_max))
}

/// Computes the nth derivative of a function using finite differences.
///
/// # Arguments
///
/// - `f`: Function to find derivative of.
/// - `x`: Value of x to find derivative at.
/// - `h`: Stepsize.
/// - `order`: Order of the derivative to compute.
fn finite_difference_nth(f: FuncSingle, x: f64, h: f64, order: u8) -> Result<f64, &'static str> {
    match order {
        0 => {
            let res = f(x);
            println!("f(x) = {res}");
            Ok(res)
        }
        1 => {
            let res = (f(x + h) - f(x - h)) / (2.0 * h);
            println!("f'(x) ≈ (f(x+h) - f(x-h)) / 2h  = {res:.4e}");
            Ok(res)
        }
        2 => {
            let res = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h.powi(2));
            println!("f''(x) ≈ (f(x + h) - 2f(x) + f(x - h)) / h^2 = {res:.4e}");
            Ok(res)
        }
        3 => {
            let res = (f(x + 2.0 * h) - 2.0 * f(x + h) + 2.0 * f(x - h) - f(x - 2.0 * h))
                / (2.0 * h.powi(3));
            println!(
                "f'''(x) ≈ (f(x + 2h) - 2f(x + h) + 2f(x - h) - f(x - 2h)) / (2h^3) = {res:.4e}"
            );
            Ok(res)
        }
        n => {
            let mut sum = 0.0;
            for i in 0..=n {
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                sum += sign * choose_float(n as f64, i)? * f(x + (n as f64 - i as f64) * h)
                    / h.powi(n as i32);
            }

            println!(
                "f^{n}(x) ≈ sum_i=0^n (-1)^i * C(n, i) * f(x + (n/2 - i)*h)) / h^{n} = {sum:.4e}"
            );
            Ok(sum)
        }
    }
}

/// Richardson's Extrapolation to approximate f'(x) at a particular x.
///
///
/// # Arguments
///
/// - `f`: Function to find derivative of.
/// - `x`: Value of x to find derivative at.
/// - `h_init`: Initial stepsize.
/// - `levels`: Number of levels of extrapolation.
/// - `order`: Order of the derivative to compute.
/// - `accurate`: If true, use the accurate version of Richardson's
///   extrapolation, else use a more numerically stable version.
///   See [Mode Comparison]
///
/// # Returns
///
/// - The Richardson table. The best estimate is at the bottom-right of the
///   Richardson table: `array[levels-1][levels-1]`.
///
/// # Mode Comparison
///
/// ## Stable
///
/// - This version increases `h` for better stability with `h*2^n`. The tradeoff
///   is a worse approximation of the derivative as `h` is larger.
/// - **Pros**
///   - Avoids very small `h` -> less roundoff errors and subtractive cancellation errors.
///   - Better stability in presence of noisy or tabulated data.
/// - **Cons**
///   - Requires more extrapolation steps to reach the same precision.
///   - Might not achieve the same level of accuracy if `h` is too large initially.
///
/// ## Accurate
///
/// - This version decreases `h` for better accuracy with `h/2^n`. The tradeoff
///   is a worse stability of the derivative as `h` gets smaller and smaller.
/// - **Pros**
///   - Increased accuracy as `h -> 0`.
///   - Enables higher-order convergence with fewer steps.
/// - **Cons**
///   - Floating-point roundoff error increases when `h` becomes too small.
///   - Subtractive cancellation: `f(x + h) - f(x - h)` can suffer when `h` is tiny.
pub fn richardson(
    f: FuncSingle,
    x: f64,
    h_init: f64,
    levels: usize,
    order: u8,
    accurate: bool,
) -> Result<Vec<Vec<f64>>, &'static str> {
    let levels = levels + 1;
    let mut r = vec![vec![0.0; levels]; levels];

    println!("╭───────────────");
    println!(
        "│ Richardson's Extrapolation - {}",
        if accurate { "Accurate" } else { "Stable" }
    );
    println!("├─");
    println!("│ n = {order}");
    println!("│ x = {x}");

    // Fill first column: central difference approximations at increasing h
    for (i, row) in r.iter_mut().enumerate() {
        let h = if accurate {
            print!("│ h = {h_init} / 2^{i}, ");
            h_init / 2f64.powi(i as i32)
        } else {
            print!("│ h = {h_init} * 2^{i}, ");
            h_init * 2f64.powi(i as i32)
        };
        row[0] = finite_difference_nth(f, x, h, order)?;
    }

    println!("├─");

    // Apply Richardson extrapolation
    for j in 1..levels {
        let mut power_of_4 = 1.0;
        for i in j..levels {
            power_of_4 *= 4.0;

            if accurate {
                r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (power_of_4 - 1.0);
                println!(
                    "│ f^{order}(x) = {:.4e} + ({:.4e} - {:.4e}) / ({power_of_4} - 1) = {:.6e}",
                    r[i][j - 1],
                    r[i][j - 1],
                    r[i - 1][j - 1],
                    r[i][j]
                );
            } else {
                r[i][j] = r[i - 1][j - 1] + (r[i - 1][j - 1] - r[i][j - 1]) / (power_of_4 - 1.0);
                println!(
                    "│ f^{order}(x) = {:.4e} + ({:.4e} - {:.4e}) / {} = {:.6e} + O(h^{})",
                    r[i - 1][j - 1],
                    r[i - 1][j - 1],
                    r[i][j - 1],
                    power_of_4 - 1.0,
                    r[i][j],
                    2 * (j + 1)
                );
            }
        }
    }

    println!("├─");
    println!(
        "│ Best estimate: {:.6e} + O(h^{})",
        r[levels - 1][levels - 1],
        2 * levels
    );
    println!("╰───────────────");

    Ok(r)
}
