use crate::utils::{compute_derivative, table_formatter, FuncSingle};

/// Compute the number of iterations needed to reach the error threshold using the bissection method.
pub fn get_bissection_iterations(interval: (f64, f64), max_error: f64) -> usize {
    ((interval.1 - interval.0).abs() / max_error).log2().ceil() as usize
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
pub fn bissection(
    f: FuncSingle,
    interval: (f64, f64),
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    // Get the number of iterations
    let iterations = iterations_max.min(get_bissection_iterations(interval, x_tolerance));

    // If the number of iterations is 0, return the middle of the interval
    if iterations == 0 {
        return Ok((interval.0 + interval.1) / 2.0);
    }

    // Initialize the interval
    let mut a = interval.0;
    let mut b = interval.1;
    let f_a = f(a);
    let f_b = f(b);

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("The function does not change sign in the interval");
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

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("The function does not change sign in the interval");
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
        if f_x.abs() < f_x_tolerance || (b - a).abs() < x_tolerance {
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

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        return Err("The function does not change sign in the interval");
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

    // Compute Newton's method and prepare the rows of the table
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(iterations_max);
    for i in 0..iterations_max {
        // Update the values
        x1 = x2;
        x2 = x1 - f_x / f_x_derivative;
        f_x = f(x2);
        f_x_derivative = compute_derivative(f, x2);

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

/// This function computes the root of a function using the fixed point method.
///
/// # Arguments
///
/// - `f`: A function that takes a single f64 argument and returns a f64.
/// - `x_initial`: The initial guess for the root.
/// - `x_tolerance`: The maximum error allowed for the root.
/// - `f_x_tolerance`: The maximum error allowed for the function value at the root.
/// - `iterations_max`: The maximum number of iterations to run the Newton's method.
pub fn fixed_point(
    f: FuncSingle,
    x_initial: f64,
    x_tolerance: f64,
    f_x_tolerance: f64,
    iterations_max: usize,
) -> Result<f64, &'static str> {
    if compute_derivative(f, x_initial) >= 1.0 {
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
