use crate::{derivation::finite_difference_nth, utils::FuncSingle};

/// Computes the integral of a Newton-Gregory Forward Polynomial Interpolation
/// using the Newton-Cotes formula.
///
/// This method assumes all x-coordinates are equally spaced.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The y-coordinates of the data points.
/// - `range`: A tuple (a, b) such that the integral is computed from x_a to x_b.
pub fn newton_cotes(
    degree: usize,
    xs: &[f64],
    fs: &[f64],
    range: (usize, usize),
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
    if degree > 4 {
        return Err("Degree not supported");
    }
    if range.0 > range.1 {
        return Err("The range is invalid.");
    }
    if range.1 > xs.len() - 1 || range.0 > xs.len() - 1 {
        return Err("The range is out of bounds.");
    }

    let a = xs[range.0];
    let b = xs[range.1];
    let a_index = range.0;
    let b_index = range.1;
    let length = b_index - a_index;

    let h = (b - a) / degree as f64;

    println!("╭───────────────");
    println!("│ Newton-Gregory Forward Polynomial Interpolation Integral");
    println!("├─");
    println!("│ Degree: {}", degree);
    println!("│ h: {}", h);
    println!("├─");

    let integral = match degree {
        1 => {
            let result = (fs[a_index] + fs[b_index]) * h / 2.0;

            // Print the result
            println!("│ P_1(x) = (f(x_{}) + f(x_{})) * (h/2)", a_index, b_index);
            println!("│ P_1(x) = {:.6e}", result);
            result
        }
        2 => {
            let middle = length / 2;
            let result = (fs[a_index] + 4.0 * fs[middle] + fs[b_index]) * h / 3.0;

            // Print the result
            println!(
                "│ P_2(x) = (f(x_{}) + 4f(x_{}) + f(x_{})) * (h/3)",
                a_index, middle, b_index
            );
            println!("│ P_2(x) = {:.6e}", result);
            result
        }
        3 => {
            let first_third = length / 3;
            let second_third = length * 2 / 3;
            let result =
                ((fs[a_index] + 3.0 * fs[first_third] + 3.0 * fs[second_third] + fs[b_index])
                    * (3.0 * h))
                    / 8.0;

            // Print the result
            println!(
                "│ P_3(x) = (f(x_{}) + 3f(x_{}) + 3f(x_{}) + f(x_{})) * (3h/8)",
                a_index, first_third, second_third, b_index
            );
            println!("│ P_3(x) = {:.6e}", result);
            result
        }
        4 => {
            let first_quarter = length / 4;
            let third_quarter = length * 3 / 4;
            let middle = length / 2;
            let result = ((7.0 * fs[a_index]
                + 32.0 * fs[first_quarter]
                + 12.0 * fs[middle]
                + 32.0 * fs[third_quarter]
                + 7.0 * fs[b_index])
                * (2.0 * h))
                / 45.0;

            // Print the result
            println!(
                "│ P_4(x) = (7f(x_{}) + 32f(x_{}) + 12f(x_{}) + 32f(x_{}) + 7f(x_{})) * (2h/45)",
                a_index, first_quarter, middle, third_quarter, b_index
            );
            println!("│ P_4(x) = {:.6e}", result);

            result
        }
        _ => unreachable!(),
    };

    println!("╰───────────────");
    Ok(integral)
}

pub fn newton_cotes_error() {
    todo!("Ch.5 p.19-22");
}

pub enum QuadratureMethod {
    Trapezoidal,
    Simpson13,
    Simpson38,
}

impl core::fmt::Display for QuadratureMethod {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            QuadratureMethod::Trapezoidal => write!(f, "Trapezoidal"),
            QuadratureMethod::Simpson13 => write!(f, "Simpson's 1/3"),
            QuadratureMethod::Simpson38 => write!(f, "Simpson's 3/8"),
        }
    }
}

/// Computes the integral of a function using the Composite Quadrature Rule.
///
/// # Arguments
///
/// - `f`: The function to integrate.
/// - `h`: The step size.
/// - `range`: A tuple (a, b) such that the integral is computed from x_a to x_b.
/// - `method`: The quadrature method to use.
pub fn composite_quadrature(
    f: FuncSingle,
    h: f64,
    range: (f64, f64),
    method: QuadratureMethod,
) -> Result<f64, &'static str> {
    if range.0 > range.1 {
        return Err("The range is invalid.");
    }

    let num_intervals = ((range.1 - range.0) / h).ceil() as usize;
    let num_points = num_intervals - 1;
    let mut points: Vec<f64> = Vec::with_capacity(num_points);

    // Fill points from range.0 to range.1 with step h
    // We skip the first and last points
    for i in 1..=num_points {
        points.push(f(range.0 + i as f64 * h));
    }

    println!("╭───────────────");
    println!("│ Composite Quadrature");
    println!("│ Method: {}", method);
    println!("├─");
    println!("│ h: {}", h);
    println!("│ num_intervals: {}", num_intervals);
    println!("│ num_points: {} + 2", num_points);
    println!("│ range: ({}, {})", range.0, range.1);
    println!("├─");

    let integral = match method {
        QuadratureMethod::Trapezoidal => {
            println!("│ The error on the integral is O(h^2).");
            println!("├─");
            println!("│ Integral: h * (0.5f_a + f_1 + ... + f_(n-1) + 0.5f_b)");
            (0.5 * f(range.0) + points.iter().sum::<f64>() + 0.5 * f(range.1)) * h
        }
        QuadratureMethod::Simpson13 => {
            println!("│ If f(x) is a polynomial of degree n <= 3, then the integral is exact.");
            println!("│ The error on the integral is O(h^4).");
            println!("├─");
            println!("│ Integral: h/3 * (f_a + 4f_1 + 2f_2 + 4f_3 + ... + 4f_(n-1) + f_b)");

            (f(range.0)
                + points
                    .iter()
                    .enumerate()
                    // `i` is indexed from 0, but we want to start from 1
                    // so we inverse the `if` and `else` statements
                    .map(|(i, x)| if i % 2 == 0 { 4.0 * x } else { 2.0 * x })
                    .sum::<f64>()
                + f(range.1))
                * h
                / 3.0
        }
        QuadratureMethod::Simpson38 => {
            println!("│ The error on the integral is O(h^4).");
            println!("├─");
            println!(
                "│ Integral: 3h/8 * (f_a + 3f_1 + 3f_2 + + 2f_3 ... + 2f_(n-3) + 3f_(n-2) + 3f_(n-1) + f_b)"
            );

            (f(range.0)
                + points
                    .iter()
                    .enumerate()
                    // `i` is indexed from 0, but we want to start from 1
                    // so we add 1 to `i` to start from 1
                    .map(|(i, x)| if (i + 1) % 3 == 0 { 2.0 * x } else { 3.0 * x })
                    .sum::<f64>()
                + f(range.1))
                * (3.0 * h)
                / 8.0
        }
    };

    println!("│ Integral: {:.6e}", integral);
    println!("╰───────────────");

    Ok(integral)
}

/// Computes the error estimate for the Composite Trapezoidal Rule.
///
/// # Arguments
///
/// - `second_derivative`: Approximation of the maximum value of the second
///   derivative of the function in the integration `range`.
/// - `h`: The step size.
/// - `range`: A tuple (a, b) such that the integral is computed from x_a to x_b.
/// - `target_error`: The target error for the integration.
///
/// # Returns
///
/// - The maximum error estimate for the Composite Trapezoidal Rule.
pub fn composite_quadrature_trapezoidal_error(
    f: FuncSingle,
    h: f64,
    range: (f64, f64),
    target_error: f64,
) -> Result<f64, &'static str> {
    if range.0 > range.1 {
        return Err("The range is invalid.");
    }
    if target_error <= 0.0 {
        return Err("The target error must be positive.");
    }

    let num_intervals = ((range.1 - range.0) / h).ceil() as usize;
    let num_points = num_intervals - 1;
    let mut points: Vec<f64> = Vec::with_capacity(num_points);

    // Fill points from range.0 to range.1 with step h
    // We skip the first and last points
    for i in 1..=num_points {
        points.push(range.0 + i as f64 * h);
    }

    println!("╭───────────────");
    println!("│ Composite Trapezoidal Rule Error Estimate");
    println!("├─");
    println!("│ h: {}", h);
    println!("│ range: ({}, {})", range.0, range.1);
    println!("├─");
    println!("│ Testing derivatives:");

    // Find the maximum value of the second derivative
    let mut max_derivative = 0.0;
    let mut max_derivative_point = 0.0;

    // Test first point
    let derivative = (f(range.0 + 2.0 * h) - 2.0 * f(range.0 + h) + f(range.0)) / (h.powi(2));
    if derivative.abs() > max_derivative {
        max_derivative = derivative.abs();
        max_derivative_point = range.0;
    }
    println!(
        "│ x = {}, f''(x) ≈ (f(x + 2h) - 2f(x + h) + f(x)) / h^2 = {derivative:.4e}",
        range.0
    );

    // Test middle points
    for val in points.iter() {
        print!("│ x = {}, ", *val);
        let derivative = finite_difference_nth(f, *val, h, 2)?.abs();
        if derivative > max_derivative {
            max_derivative = derivative;
            max_derivative_point = *val;
        }
    }

    // Test last point
    let derivative = (f(range.1 - 2.0 * h) - 2.0 * f(range.1 - h) + f(range.1)) / (h.powi(2));
    if derivative.abs() > max_derivative {
        max_derivative = derivative.abs();
        max_derivative_point = range.1;
    }
    println!(
        "│ x = {}, f''(x) ≈ (f(x - 2h) - 2f(x - h) + f(x)) / h^2 = {derivative:.4e}",
        range.1
    );

    // Print max derivative
    println!("├─");
    println!(
        "│ Max derivative: {} at x = {}",
        max_derivative, max_derivative_point
    );

    // Compute the maximum error
    let range_length = range.1 - range.0;
    let e_max = (-h.powi(2) * range_length * max_derivative.abs()).abs() / 12.0;

    // Find required h to reach target error
    let target_h = ((12.0 * target_error) / (range_length * max_derivative.abs())).sqrt();
    // Find required number of intervals to reach target error
    let target_n_intervals = (range_length / target_h).ceil();

    println!("├─");
    println!("│ |E| <={:.6e}", e_max);
    println!("├─");
    println!("│ Target error: {:.6e}", target_error);
    println!("│ To reach target error:");
    println!("│   h <= {:.6e}", target_h);
    println!("│   number of intervals >= {}", target_n_intervals);
    println!("╰───────────────");

    Ok(e_max)
}

pub fn composite_quadrature_simpson13_error() {
    todo!("Ch.5 p.30");
}

pub fn composite_quadrature_simpson38_error() {
    todo!("Ch.5 p.32");
}

/// Computes the integral of a function using the Romberg integration technique.
///
/// # Arguments
///
/// - `f`: The function to integrate.
/// - `range`: A tuple `(a, b)` representing the integration limits.
/// - `max_level`: Maximum depth `k` of extrapolation (the number of rows in the Romberg table).
///
/// # Returns
///
/// - Final integral estimate using Romberg extrapolation.
pub fn romberg(
    f: FuncSingle,
    h_init: f64,
    range: (f64, f64),
    levels: usize,
) -> Result<Vec<Vec<f64>>, &'static str> {
    if range.0 >= range.1 {
        return Err("Invalid integration range.");
    }

    let levels = levels + 1;
    let mut r = vec![vec![0.0; levels]; levels];

    println!("╭─────────────────────────────");
    println!("│ Romberg Integration Table");
    println!("├─");
    println!("│ Integral: h * (0.5f_a + f_1 + ... + f_(n-1) + 0.5f_b)");
    println!("├─");

    // Fill first column: trapezoidal rule
    for (i, row) in r.iter_mut().enumerate() {
        // Compute new h
        let h = h_init * 2f64.powi(i as i32);

        // Get number of steps
        let num_steps = ((range.1 - range.0) / h).ceil() as usize - 1;

        // Fill points from range.0 to range.1 with step h
        // We skip the first and last points
        let mut steps: Vec<f64> = Vec::with_capacity(num_steps);
        for i in 1..=num_steps {
            let tmp = range.0 + i as f64 * h;
            // TODO: Remove this rounding
            steps.push(f((tmp * 100.0).round() / 100.0));
        }

        // Compute the integral
        row[0] = (0.5 * f(range.0) + steps.iter().sum::<f64>() + 0.5 * f(range.1)) * h;
        println!("│ h = {h_init} * 2^{i}, Integral = {:.4e} + O(h^2)", row[0]);
    }

    println!("├─");

    // Apply Richardson extrapolation
    let mut power_of_4 = 1;
    for j in 1..levels {
        power_of_4 *= 4;
        for i in j..levels {
            r[i][j] = r[i - 1][j - 1] + (r[i - 1][j - 1] - r[i][j - 1]) / (power_of_4 - 1) as f64;
            println!(
                "│ T_({},(n/{})) = {:.4e} + ({:.4e} - {:.4e}) / {} = {:.6e} + O(h^{})",
                j,
                i,
                r[i - 1][j - 1],
                r[i - 1][j - 1],
                r[i][j - 1],
                power_of_4 - 1,
                r[i][j],
                2 * (j + 1)
            );
        }
    }

    println!("├─");
    println!(
        "│ Best estimate: {:.6e} + O(h^{})",
        r[levels - 1][levels - 1],
        2 * levels
    );
    println!("╰─────────────────────────────");

    Ok(r)
}

/// Gauss–Legendre quadrature nodes and weights for n = 1 to 5
///
/// Each entry is a tuple: (nodes: &[f64], weights: &[f64])
pub fn gauss_legendre_table(n: usize) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    match n {
        1 => Ok((vec![0.0], vec![2.0])),
        2 => {
            let x = (1.0_f64 / 3.0_f64).sqrt();
            Ok((vec![-x, x], vec![1.0, 1.0]))
        }
        3 => {
            let x = (3.0_f64 / 5.0_f64).sqrt();
            let w1 = 8.0 / 9.0;
            let w2 = 5.0 / 9.0;
            Ok((vec![-x, 0.0, x], vec![w2, w1, w2]))
        }
        4 => {
            let a = 3.0 / 7.0;
            let b = (2.0 / 7.0) * (6.0_f64 / 5.0_f64).sqrt();
            let x1 = (a - b).sqrt();
            let x2 = (a + b).sqrt();

            let sqrt_30 = 30.0_f64.sqrt();
            let w1 = (18.0 + sqrt_30) / 36.0;
            let w2 = (18.0 - sqrt_30) / 36.0;

            Ok((vec![-x2, -x1, x1, x2], vec![w2, w1, w1, w2]))
        }
        5 => {
            let a = 1.0 / 3.0;
            let b = 2.0_f64 * (10.0_f64 / 7.0_f64).sqrt();
            let x1 = a * (5.0_f64 - b).sqrt();
            let x2 = a * (5.0_f64 + b).sqrt();

            let sqrt_70 = 70.0_f64.sqrt();
            let w0 = 128.0 / 225.0;
            let w1 = (322.0 + 13.0 * sqrt_70) / 900.0;
            let w2 = (322.0 - 13.0 * sqrt_70) / 900.0;

            Ok((vec![-x2, -x1, 0.0, x1, x2], vec![w2, w1, w0, w1, w2]))
        }
        _ => Err("Supported n: 1 to 5"),
    }
}

/// Computes the integral of a function using Gaussian quadrature.
///
/// # Arguments
///
/// - `f`: The function to integrate.
/// - `degree`: The degree of the polynomial to use.
/// - `range`: A tuple (a, b) such that the integral is computed from x_a to x_b.
pub fn gaussian_quadrature(
    f: FuncSingle,
    degree: usize,
    range: (f64, f64),
) -> Result<f64, &'static str> {
    if degree < 1 {
        return Err("The degree must be at least 1.");
    }
    if range.0 > range.1 {
        return Err("The range is invalid.");
    }

    let legendre = gauss_legendre_table(degree)?;

    // Rescale between -1 and 1
    let factor = (range.1 - range.0) / 2.0;
    // Change of variable to rescale
    let xt = |t: f64, a: f64, b: f64| ((b - a) * t + (a + b)) * 0.5;

    println!("╭───────────────");
    println!("│ Gaussian Quadrature");
    println!("├─");
    println!("│ Degree: {}", degree);
    println!("│ range: ({}, {})", range.0, range.1);
    println!("├─");
    println!("│ Legendre nodes (t_i): {:?}", legendre.0);
    println!("│ Legendre weights (w_i): {:?}", legendre.1);
    println!("├─");
    println!("│ Variable change: x(t) = 0.5 * ((b - a)t + (a + b))");

    // Compute the integral
    let integral: f64 = legendre
        .0
        .iter()
        .zip(legendre.1.iter())
        .enumerate()
        .map(|(i, (t, w))| {
            let x = xt(*t, range.0, range.1);
            let ft = f(x);
            let res = ft * w;
            println!(
                "│ f(x(t_{})) * w_{} = f({:.4e}) * w_{} = {:.4e} * w_{} = {:.4e}",
                i, i, x, i, ft, i, res
            );
            res
        })
        .sum();

    let result = factor * integral;

    println!("├─");
    println!("│ Integral: w_0 * f(x_0) + ... + w_n * f(x_n)");
    println!("│ Factor: dx = 0.5 * (b - a) dt = {} dt", factor);
    println!("│ Integral: Factor * ({:.6e}) = {:.6e}", integral, result);
    println!("╰───────────────");

    Ok(result)
}
