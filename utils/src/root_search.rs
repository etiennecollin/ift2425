use crate::utils::{compute_derivative, table_formatter, FuncSingle};

/// Compute the number of iterations needed to reach the error threshold using the bissection method.
///
/// # Arguments
///
/// - `interval`: A tuple representing the interval where the root is located.
/// - `max_error`: The maximum error allowed for the root.
///
/// # Errors
///
/// - If the interval is such that a >= b.
pub fn get_bissection_iterations(
    interval: (f64, f64),
    max_error: f64,
) -> Result<usize, &'static str> {
    // Check that a < b
    if interval.0 >= interval.1 {
        return Err("The interval must be such that a < b");
    }

    let iterations = ((interval.1 - interval.0).abs() / max_error).log2().ceil() as usize;

    println!("╭───────────────");
    println!("│ Bissection method");
    println!("├─");
    println!("│ |b - a| / 2^n = error");
    println!(
        "│ |{:.3} - {:.3}| / 2^n = {:.3e}",
        interval.1, interval.0, max_error
    );
    println!("├─");
    println!("│ Convergence after {} iterations", iterations);
    println!("╰───────────────");

    Ok(iterations)
}

/// This function computes the root of a function using the bissection method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `interval`: A tuple representing the interval where the root is located.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the bissection method.
///
/// # Errors
///
/// - If the interval is such that a >= b.
/// - If the function does not cross the x-axis (f(a)f(b) < 0).
pub fn bissection(
    f: FuncSingle,
    interval: (f64, f64),
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    // Initialize the interval
    let mut a = interval.0;
    let mut b = interval.1;

    // Check that a < b
    if a >= b {
        return Err("The interval must be such that a < b");
    }

    // If the number of iterations is 0, return the middle of the interval
    let iterations = iterations_max.min(get_bissection_iterations(interval, x_tolerance)?);
    if iterations == 0 {
        return Ok((interval.0 + interval.1) / 2.0);
    }

    let f_a = f(a);
    let f_b = f(b);

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("f(a)f(b) must be less than zero");
    }

    // Compute the bissection method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations);
    let mut x;
    let mut f_x;
    for i in 0..iterations {
        // Compute the middle of the interval and the function value
        x = (a + b) / 2.0;
        f_x = f(x);

        // Store the information for the table
        rows.push(
            format!(
                "{}, {}, {}, {}, {:e}, {:e}",
                i + 1,
                x,
                a,
                b,
                (b - a).abs() / 2.0,
                f_x,
            )
            .split(",")
            .map(|s| s.trim().to_owned())
            .collect(),
        );

        // Check if the function changes sign in the interval
        if f(a) * f_x <= 0.0 {
            b = x;
        } else {
            a = x;
        }

        // Check if the function value is close enough to zero
        if f_x.abs() < f_x_tolerance || (b - a).abs() < x_tolerance {
            break;
        }
    }
    x = (a + b) / 2.0;

    // Prepare the header and information
    let info = vec![
        format!("Running bissection with {} iterations", iterations),
        format!(
            "Assuming f is continuous on the interval [{}, {}]",
            interval.0, interval.1
        ),
        format!("f(a) = {:e}, f(b) = {:e}", f_a, f_b),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("x = {}", x),
    ];
    let header = "i, x, a, b, |b-a|/2^i, f(x)"
        .split(",")
        .map(|s| s.trim().to_owned())
        .collect();

    // Format and print the table
    let table = table_formatter(info, header, rows)?;
    println!("{}", table);

    Ok(x)
}

/// This function computes the root of a function using the linear interpolation method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `interval`: A tuple representing the interval where the root is located.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the linear interpolation method.
///
/// # Errors
///
/// - If the interval is such that a >= b.
/// - If the function does not cross the x-axis (f(a)f(b) < 0).
pub fn linear_interpolation(
    f: FuncSingle,
    interval: (f64, f64),
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    // Initialize the interval
    let mut a = interval.0;
    let mut b = interval.1;
    let mut f_a = f(a);
    let mut f_b = f(b);

    // Check that a < b
    if a >= b {
        return Err("The interval must be such that a < b");
    }

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("f(a)f(b) must be less than zero");
    }

    // Initialize the values
    let mut x = b - (f_b * (b - a)) / (f_b - f_a);
    let mut f_x = f(x);

    // Compute the linear interpolation method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    for i in 0..iterations_max {
        // Store the information for the table
        rows.push(
            format!(
                "{}, {}, {:e}, {:e}, {:e}, {:e}",
                i + 1,
                x,
                f_x,
                (b - a).abs(),
                (x - a).abs(),
                (x - b).abs(),
            )
            .split(",")
            .map(|s| s.trim().to_owned())
            .collect(),
        );

        // Check if the error is below the tolerance
        if f_x.abs() < f_x_tolerance || (x - a).abs() < x_tolerance {
            break;
        }

        // Check if the function changes sign in the interval
        if f(a) * f_x <= 0.0 {
            b = x;
        } else {
            a = x;
        }

        // Update the values
        f_a = f(a);
        f_b = f(b);
        x = b - (f_b * (b - a)) / (f_b - f_a);
        f_x = f(x);
    }

    // Prepare the header and information
    let info = vec![
        "Running linear interpolation method".to_owned(),
        format!(
            "Assuming f is continuous on the interval [{}, {}]",
            interval.0, interval.1
        ),
        format!("f(a) = {:e}, f(b) = {:e}", f_a, f_b),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let header = "i, x, f(a), |b-a|, |x-a|, |x-b|"
        .split(",")
        .map(|s| s.trim().to_owned())
        .collect();

    // Format and print the table
    let table = table_formatter(info, header, rows)?;
    println!("{}", table);

    Ok(x)
}

/// This function computes the root of a function using the secant method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `interval`: A tuple representing the interval where the root is located.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the linear interpolation method.
///
/// # Errors
///
/// - If the interval is such that a >= b.
/// - If the function does not cross the x-axis (f(a)f(b) < 0).
pub fn secant(
    f: FuncSingle,
    interval: (f64, f64),
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    // Initialize the interval
    let mut a = interval.0;
    let mut b = interval.1;
    let mut f_a = f(a);
    let mut f_b = f(b);

    // Check that a < b
    if a >= b {
        return Err("The interval must be such that a < b");
    }

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("f(a)f(b) must be less than zero");
    }

    // Initialize the values
    let mut x = b - (f_b * (b - a)) / (f_b - f_a);
    let mut f_x = f(x);

    // Compute the secant method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    for i in 0..iterations_max {
        // Store the information for the table
        rows.push(
            format!(
                "{}, {}, {}, {}, {:e}, {:e}",
                i + 1,
                a,
                b,
                x,
                f_x,
                (b - a).abs(),
            )
            .split(",")
            .map(|s| s.trim().to_owned())
            .collect(),
        );

        // Check if the error is below the tolerance
        if f_x.abs() < f_x_tolerance || (b - a).abs() < x_tolerance {
            break;
        }

        // Check if the function changes sign in the interval

        a = b;
        b = x;

        // Update the values
        f_a = f(a);
        f_b = f(b);
        x = b - (f_b * (b - a)) / (f_b - f_a);
        f_x = f(x);
    }

    // Prepare the header and information
    let info = vec![
        "Running secant method".to_owned(),
        format!(
            "Assuming f is continuous on the interval [{}, {}]",
            interval.0, interval.1
        ),
        format!("f(a) = {:e}, f(b) = {:e}", f_a, f_b),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x),
    ];
    let header = "i, a, b, x, f(x), |b-a|"
        .split(",")
        .map(|s| s.trim().to_owned())
        .collect();

    // Format and print the table
    let table = table_formatter(info, header, rows)?;
    println!("{}", table);

    Ok(x)
}

/// This function computes the root of a function using Newton's method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `x_initial`: The initial guess for the root.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the Newton's method.
///
/// # Errors
///
/// - If the derivative of the function at the initial point is equal to zero.
/// - If the derivative of the function is equal to zero around the root.
pub fn newton(
    f: FuncSingle,
    x_initial: f64,
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    // Initialize the variables
    let mut x1;
    let mut x2 = x_initial;
    let mut f_x = f(x_initial);
    let mut f_x_derivative = compute_derivative(f, x_initial);

    if f_x_derivative == 0.0 {
        return Err(
            "The derivative of the function at the initial point must be different from zero",
        );
    }

    // Compute Newton's method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    for i in 0..iterations_max {
        // Update the values
        x1 = x2;
        x2 = x1 - f_x / f_x_derivative;
        f_x = f(x2);
        f_x_derivative = compute_derivative(f, x2);

        if f_x_derivative == 0.0 {
            return Err(
                "The derivative of the function must be different from zero around the root",
            );
        }

        // Store the information for the table
        rows.push(
            format!(
                "{}, {}, {}, {:e}, {:e}",
                i + 1,
                x1,
                x2,
                (x2 - x1).abs(),
                f_x,
            )
            .split(",")
            .map(|s| s.trim().to_owned())
            .collect(),
        );

        // Check if the error is below the tolerance
        if f_x.abs() < f_x_tolerance || (x2 - x1).abs() < x_tolerance || f_x_derivative == 0.0 {
            break;
        }
    }

    // Prepare the header and information
    let info = vec![
        "Running Newton's method".to_owned(),
        "Assuming f' exists".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x2),
    ];
    let header = "i, x1, x2, |x2-x1|, f(x2)"
        .split(",")
        .map(|s| s.trim().to_owned())
        .collect();

    // Format and print the table
    let table = table_formatter(info, header, rows)?;

    println!("{}", table);

    Ok(x2)
}

/// Computes the tolerance to use on x for the fixed point method to reach the error threshold.
///
/// This uses the finite increment theorem.
///
/// # Arguments
///
/// - `f`: The initial function
/// - `g`: The formulation of f as g(x) = x
/// - `interval`: A tuple representing the interval where the root is located.
/// - `error_threshold`: The maximum error allowed for the root.
///
/// # Errors
///
/// - If the interval is such that a >= b.
/// - If the function does not cross the x-axis (f(a)f(b) < 0).
/// - If the derivative of the function is greater or equal to 1.
pub fn get_fixed_point_x_tolerance_fit(
    f: FuncSingle,
    g: FuncSingle,
    interval: (f64, f64),
    error_threshold: f64,
) -> Result<f64, &'static str> {
    // Check that a < b
    if interval.0 >= interval.1 {
        return Err("The interval must be such that a < b");
    }

    // Check that the function crosses the x-axis
    if f(interval.0) * f(interval.1) >= 0.0 {
        return Err("f(a)f(b) must be less than zero");
    }

    // Get max value of g'(x) on the interval
    let derivative_a = compute_derivative(g, interval.0).abs();
    let derivative_b = compute_derivative(g, interval.1).abs();
    let max_derivative = derivative_a.max(derivative_b);

    // Check that the derivative is less than 1 for convergence
    if max_derivative >= 1.0 {
        return Err("The derivative of the function must be less than 1 for convergence");
    }

    // Compute the x tolerance
    let max_derivative_frac = 1.0 / (1.0 - max_derivative).abs();
    let x_tolerance = error_threshold / max_derivative_frac;

    // Print results
    println!("╭───────────────");
    println!("│ Finite increment theorem");
    println!("├─");
    println!("│ Assuming:");
    println!("│ That f is continuous and monotonic on the interval where the root is located");
    println!("│ That f crosses the x-axis on the interval: f(a)f(b) < 0");
    println!("│ That there is a unique root in the interval: f'(x) > 0 on the interval");
    println!("├─");
    println!("│ f(a) = {}", f(interval.0));
    println!("│ f(b) = {}", f(interval.1));
    println!("├─");
    println!("│ Using the finite increment theorem.");
    println!("│ r_(n+1) - r_n = (r_(n+1) - r) - (r_n - r)");
    println!("│               = (g(r_n) - g(r)) - (r_n - r)");
    println!("│               = ((g(r_n) - g(r)) - (r_n - r)) / (r_n - r) * (r_n - r)");
    println!("│               = (g'(ζ) - 1) * (r_n - r) with ζ between r_n and r");
    println!("│ Therefore, |r_n - r| = 1/|1 - g'(ζ)| * |r_(n+1) - r_n| for all x in the interval");
    println!(
        "│ Knowing that |g'(ζ)| ≤ {:.3}, we have that 1/|1 - g'(ζ)| ≤ {:.3}",
        max_derivative, max_derivative_frac
    );
    println!(
        "│ Hence, |r_n - r| ≤ {:.3} * |r_(n+1) - r_n|",
        max_derivative_frac,
    );
    println!("├─");
    println!("│ Therefore, if we want an error |r_n - r| < {:e}, we can choose a n (stop the algo) when:", error_threshold);
    println!(
        "│ |r_(n+1) - r_n| ≤ {:e} / {:.3} ≈ {:.3e}",
        error_threshold, max_derivative_frac, x_tolerance
    );
    println!("╰───────────────");

    Ok(x_tolerance)
}

/// Computes the number of iterations needed to reach the error threshold using the fixed point method.
///
/// This uses the mean value theorem.
///
/// # Arguments
///
/// - `f`: The initial function
/// - `g`: The formulation of f as g(x) = x
/// - `interval`: A tuple representing the interval where the root is located.
/// - `error_threshold`: The maximum error allowed for the root.
///
/// # Errors
///
/// - If the interval is such that a >= b.
/// - If the function does not cross the x-axis (f(a)f(b) < 0).
/// - If the derivative of the function is greater or equal to 1.
pub fn get_fixed_point_iterations_mvt(
    f: FuncSingle,
    g: FuncSingle,
    interval: (f64, f64),
    error_threshold: f64,
) -> Result<usize, &'static str> {
    // Check that a < b

    if interval.0 >= interval.1 {
        return Err("The interval must be such that a < b");
    }

    // Check that the function crosses the x-axis
    if f(interval.0) * f(interval.1) >= 0.0 {
        return Err("f(a)f(b) must be less than zero");
    }

    // Get max value of g'(x) on the interval
    let derivative_a = compute_derivative(g, interval.0).abs();
    let derivative_b = compute_derivative(g, interval.1).abs();
    let max_derivative = derivative_a.max(derivative_b);

    // Check that the derivative is less than 1 for convergence
    if max_derivative >= 1.0 {
        return Err("The derivative of the function must be less than 1 for convergence");
    }

    // Compute the number of iterations
    let max_error = interval.1 - interval.0;
    let iterations = (error_threshold / max_error).log(max_derivative).ceil() as usize;

    // Print results
    println!("╭───────────────");
    println!("│ Mean value theorem");
    println!("├─");
    println!("│ Assuming:");
    println!("│ That f is continuous and monotonic on the interval where the root is located");
    println!("│ That f crosses the x-axis on the interval: f(a)f(b) < 0");
    println!("│ That there is a unique root in the interval: f'(x) > 0 on the interval");
    println!("├─");
    println!("│ f(a) = {}", f(interval.0));
    println!("│ f(b) = {}", f(interval.1));
    println!("├─");
    println!("│ Using the mean value theorem.");
    println!("│ r_n - r = (g(r_(n-1)) - g(r)) / (r_(n-1) - r) * (r_(n-1) - r)");
    println!("│         = g'(ζ) * (r_(n-1) - r) with ζ in interval");
    println!("│ Knowing that |g'(ζ)| ≤ {:.3}", max_derivative);
    println!(
        "│ We have |r_n - r| ≤ {:.3} * |r_(n-1) - r| ≤ ... ≤ {:.3}^n * |r_0 - r| ≤ {:.3}^n * (b-a) = {:.3}^n * {:.3e}",
        max_derivative, max_derivative, max_derivative, max_error, max_error
    );
    println!(
        "│ Therefore, |r_n - r| ≤ {:.3}^n * {:.3e} ≤ {:.3e}",
        max_derivative, max_error, error_threshold
    );
    println!("├─");
    println!("│ Convergence after a maximum of {} iterations", iterations);
    println!("╰───────────────");

    Ok(iterations)
}

/// This function computes the root of a function using the fixed point method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `x_initial`: The initial guess for the root.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the Newton's method.
///
/// # Errors
///
/// - If the derivative of the function at the initial point is greater or equal to 1.
pub fn fixed_point(
    f: FuncSingle,
    x_initial: f64,
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    if compute_derivative(f, x_initial).abs() >= 1.0 {
        return Err("The derivative of the function at the initial point must be less than 1 for convergence");
    }

    // Initialize the variables
    let mut x1 = x_initial;
    let mut x2 = f(x1);

    // Compute the fixed point method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    for i in 0..iterations_max {
        // Store the information for the table
        rows.push(
            format!("{}, {}, {}, {:e}", i + 1, x1, x2, (x2 - x1).abs(),)
                .split(",")
                .map(|s| s.trim().to_owned())
                .collect(),
        );

        // Check if the error is below the tolerance
        if f(x2).abs() < f_x_tolerance || (x2 - x1).abs() < x_tolerance {
            break;
        }

        // Update the values
        x1 = x2;
        x2 = f(x1);
    }

    // Prepare the header and information
    let info = vec![
        "Running fixed point method".to_owned(),
        "Assuming |g'(x_initial)| < 1 and x = g(x)".to_owned(),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("max iterations = {}", iterations_max),
        format!("x = {}", x2),
    ];
    let header = "i, x1, x2, |x2-x1|"
        .split(",")
        .map(|s| s.trim().to_owned())
        .collect();

    // Format and print the table
    let table = table_formatter(info, header, rows)?;
    println!("{}", table);

    Ok(x2)
}
