use crate::utils::{table_formatter, FuncSingle};

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
/// - `min_iterations`: The minimum number of iterations to run the bissection method.
///     - If [None], the number of iterations is computed using the interval and the max error to reach the error threshold.
///     - If [Some], the algorithm will run for at least `min_iterations` and until the error threshold is reached.
pub fn bissection(
    f: FuncSingle,
    interval: (f64, f64),
    x_tolerance: f64,
    f_x_tolerance: f64,
    min_iterations: Option<usize>,
) -> f64 {
    // Get the number of iterations
    // The minimum number of iterations is computed using the interval and the max error
    let iterations = match min_iterations {
        Some(min) => min.max(get_bissection_iterations(interval, x_tolerance)),
        None => get_bissection_iterations(interval, x_tolerance),
    };

    // If the number of iterations is 0, return the middle of the interval
    if iterations == 0 {
        return (interval.0 + interval.1) / 2.0;
    }

    // Initialize the interval
    let mut a = interval.0;
    let mut b = interval.1;
    let f_a = f(a);
    let f_b = f(b);

    // Check if the function changes sign in the interval
    if f_a * f_b >= 0.0 {
        panic!("The function does not change sign in the interval");
    }

    // Compute the bissection method and prepare the rows of the table
    let mut rows: Vec<String> = Vec::with_capacity(iterations);
    let mut x;
    for i in 0..iterations {
        // Compute the middle of the interval and the function value
        x = (a + b) / 2.0;
        let f_x = f(x);

        // Store the information for the table
        rows.push(format!(
            "{}, {}, {}, {}, {}, {:e}",
            i + 1,
            x,
            f_x,
            a,
            b,
            (b - a).abs() / 2.0
        ));

        // Check if the function changes sign in the interval
        if f(a) * f_x <= 0.0 {
            b = x;
        } else {
            a = x;
        }

        // Check if the function value is close enough to zero
        if f_x.abs() < f_x_tolerance {
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
        format!("f(a) = {}, f(b) = {}", f_a, f_b),
        format!("x tolerance = {:e}", x_tolerance),
        format!("f(x) tolerance = {:e}", f_x_tolerance),
        format!("x = {}", x),
    ];
    let header = "iteration, x, f(x), a, b, |b-a|/2^i".to_string();

    // Format and print the table
    let table = table_formatter(info, header, rows);
    println!("{}", table);

    x
}
