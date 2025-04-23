use core::f64;

use nalgebra::{Matrix2, Matrix4, Vector2};

use crate::utils::{FuncSingle, choose_float, factorial, table_formatter};

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

    println!("╭───────────────");
    println!("│ Lagrange Polynomial Interpolation");
    println!("├─");
    println!("│ x = {}, degree = {}", x, degree);
    println!(
        "│ P_{}({}) = sum_(i=0)^{} π_i(x)/π_i(x_i) * y_i",
        degree, x, degree
    );
    println!("│ π_i(x) = product_(j=0, j≠i)^{} (x - x_j)", degree);
    println!("├─");
    println!("│ sum = 0.0");

    let mut sum = 0.0;
    for i in 0..=degree {
        let a = lagrange_basis_polynomial(degree, i, x, xs)?;
        let b = lagrange_basis_polynomial(degree, i, xs[i], xs)?;
        sum += a / b * fs[i];
        println!(
            "│ sum += π_{}({}) / π_{}({}) * y_{} = {:.4e} / {:.4e} * {:.4e}",
            i, x, i, xs[i], i, a, b, fs[i]
        );
    }

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
fn lagrange_basis_polynomial(
    degree: usize,
    i: usize,
    x: f64,
    xs: &[f64],
) -> Result<f64, &'static str> {
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

/// Computes the finite forward difference table for a given degree and function
/// values.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `fs`: The function values.
///
/// # Returns
///
/// - An array that is `degree+1 x degree+1` such that `array[i][j]` contains `delta^j f_i`.
pub fn finite_forward_diff_table(degree: usize, fs: &[f64]) -> Result<Vec<Vec<f64>>, &'static str> {
    let n = degree + 1;

    if degree > fs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let mut array = vec![vec![0.0; degree + 1]; n];

    // Initialize the first column of the array with the function values
    for i in 0..n {
        array[i][0] = fs[i];
    }

    // Fill the finite difference table
    for j in 1..n {
        for i in 0..n - j {
            array[i][j] = array[i + 1][j - 1] - array[i][j - 1];
        }
    }

    Ok(array)
}

/// Computes the finite centered (or divided) difference table for a given
/// degree and function values.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The function values.
///
/// # Returns
///
/// - An array that is `degree+1 x degree+1` such that `array[i][j]` contains `delta^j f_i`.
pub fn finite_centered_diff_table(
    degree: usize,
    xs: &[f64],
    fs: &[f64],
) -> Result<Vec<Vec<f64>>, &'static str> {
    if xs.len() != fs.len() {
        return Err("The length of xs and fs must be equal.");
    }
    if degree > xs.len() - 1 {
        return Err("The degree must be less than the number of data points minus 1.");
    }

    let n_rows = degree + 1;
    let n_cols = degree + 2;
    let mut array = vec![vec![0.0; n_cols]; n_rows];

    // Initialize the first columns of the array with xs and fs
    for i in 0..n_rows {
        array[i][0] = xs[i];
        array[i][1] = fs[i];
    }

    // Fill the finite difference table
    for j in 2..n_cols {
        for i in 0..(n_rows - j + 1) {
            let numerator = array[i + 1][j - 1] - array[i][j - 1];
            let denominator = array[i + j - 1][0] - array[i][0];
            if denominator == 0.0 {
                return Err("Duplicate x values detected.");
            }
            array[i][j] = numerator / denominator;
        }
    }

    Ok(array)
}

/// Formats the finite centered (or divided) difference table into a string.
/// This string, once printed, will look like a table.
///
/// # Arguments
///
/// - `table`: The finite centered difference table.
pub fn finite_centered_diff_table_string(table: &[Vec<f64>]) -> Result<String, &'static str> {
    let n_rows = table.len();
    let n_cols = table[0].len();
    let degree = n_cols - 2;

    let info = vec![
        "Finite Centered Difference Table".to_owned(),
        format!("Degree: {}", degree),
    ];

    // Format the header
    let mut header = vec!["x".to_owned(), "y".to_owned()];
    for i in 1..=degree {
        let mut content = Vec::new();
        for j in 0..=i {
            content.push(format!("x{}", j));
        }
        header.push(format!("∆^{}[{}]", i, content.join(", ")));
    }

    // Format the table rows
    let mut rows = Vec::with_capacity(n_rows);
    for (i, table_row) in table.iter().enumerate() {
        let mut row = Vec::with_capacity(n_cols);
        for (j, &value) in table_row.iter().enumerate() {
            // Only print the values that are not empty
            // (That's the upper triangular part of the table)
            if j < n_cols - i {
                row.push(format!("{:>2}", value));
            } else {
                row.push("".to_owned());
            }
        }
        rows.push(row);
    }

    // Format and print the table
    let table_string = table_formatter(info, header, rows)?;

    Ok(table_string)
}

/// Formats the finite forward difference table into a string.
/// This string, once printed, will look like a table.
///
/// # Arguments
///
/// - `xs`: The x-coordinates of the data points.
/// - `table`: The finite forward difference table.
pub fn finite_forward_diff_table_string(
    xs: &[f64],
    table: &[Vec<f64>],
) -> Result<String, &'static str> {
    let n_rows = table.len();
    let n_cols = table[0].len();
    let degree = n_cols - 1;

    if xs.len() != n_rows {
        return Err("The length of xs and table must be equal.");
    }

    let info = vec![
        "Finite Forward Difference Table".to_owned(),
        format!("Degree: {}", degree),
    ];

    // Format the header
    let mut header = vec!["x".to_owned(), "y".to_owned()];
    for i in 1..=degree {
        header.push(format!("∆^{}y", i));
    }

    // Format the table rows
    let mut rows = Vec::with_capacity(n_rows);
    for (i, table_row) in table.iter().enumerate() {
        let mut row = vec![xs[i].to_string()];
        for (j, &value) in table_row.iter().enumerate() {
            // Only print the values that are not empty
            // (That's the upper triangular part of the table)
            if j < n_cols - i {
                row.push(format!("{:>2}", value));
            } else {
                row.push("".to_owned());
            }
        }
        rows.push(row);
    }

    // Format and print the table
    let table_string = table_formatter(info, header, rows)?;

    Ok(table_string)
}

/// Computes the Forward Newton-Gregory interpolation polynomial at a given point `x`.
///
/// This method assumes all x-coordinates are equally spaced.
///
/// # Arguments
///
/// - `degree`: The degree of the polynomial.
/// - `x`: The point at which to evaluate the polynomial.
/// - `h`: The spacing between the x-coordinates.
/// - `xs`: The x-coordinates of the data points.
/// - `fs`: The y-coordinates of the data points.
/// - `choose`: If true, use the choose function to compute the finite differences.
pub fn newton_gregory_forward(
    degree: usize,
    x: f64,
    h: f64,
    xs: &[f64],
    fs: &[f64],
    choose: bool,
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
    let s = (x - xs[0]) / h;

    // Newton-Gregory Forward Polynomial Interpolation
    let finite_diff = finite_forward_diff_table(degree, fs)?;
    let mut result = fs[0];
    for k in 1..=degree {
        if choose {
            result += choose_float(s, k as u8)? * finite_diff[0][k];
        } else {
            let factor = finite_diff[0][k] / factorial(k as u8)?;
            let mut product = 1.0;
            for j in 0..k {
                product *= s - j as f64
            }
            result += factor * product
        }
    }

    println!("╭───────────────");
    println!("│ Newton-Gregory Forward Polynomial Interpolation");
    println!("├─");
    println!("│ Assuming equally spaced x-coordinates");
    if choose {
        println!("│ P_(n+1)(x) = P_n(x) + `s choose n+1` * delta^(n+1)f_0");
    }
    println!("├─");
    println!("│ P_{}({}) = {}", degree, x, result);
    println!("╰───────────────");

    Ok(result)
}

pub fn newton_gregory_forward_error() -> Result<f64, &'static str> {
    todo!("Ch.4 p.21")
}

pub fn newton_divided_difference() {
    todo!("Ch.4 p.26")
}

pub enum SplineType {
    Type1,
    Type2,
    Type3,
    Type4,
}

pub fn spline(xs: &[f64], fs: &[f64], spline_type: SplineType) {
    // Compute the h (distance between each points)
    // Compute the S values
    // Return the S values and Q equations
    todo!("Ch.4 p.35,45,46")
}

pub fn tridiagonal_solve() {
    todo!("Ch.4 p.55")
}

pub fn bezier() {
    todo!("Ch.4 p.56")
}

pub fn bezier_mat() {
    todo!("Ch.4 p.63")
}

pub fn bezier_3d() {
    todo!("Ch.4 p.65")
}

pub fn bspline() {
    todo!("Ch.4 p.68")
}

pub fn bspline_mat() {
    todo!("Ch.4 p.69")
}

pub fn bspline_3d() {
    todo!("Ch.4 p.74")
}

pub fn collocation_polynomial_surface() {
    todo!("Ch.4 p.75")
}

/// Least Squares Interpolation
///
/// # Arguments
///
/// - `xs`: The x-coordinates of the data points.
/// - `ys`: The y-coordinates of the data points.
///
/// # Returns
///
/// - A vector containing the coefficients `[a, b]` of the linear regression `y = ax+b`.
pub fn linear_regression(xs: &[f64], fs: &[f64]) -> Result<Vector2<f64>, &'static str> {
    let sumx = xs.iter().sum::<f64>();
    let sumx_squared = xs.iter().map(|&x| x.powi(2)).sum::<f64>();
    let sumy = fs.iter().sum::<f64>();
    let sumxy = xs.iter().zip(fs).map(|(&x, &y)| x * y).sum::<f64>();
    let a = Matrix2::new(sumx_squared, sumx, sumx, xs.len() as f64);
    let b = Vector2::new(sumxy, sumy);

    let x = match a.lu().solve(&b) {
        Some(x) => x,
        None => {
            return Err("Matrix is singular");
        }
    };

    println!("╭───────────────");
    println!("│ Linear Regression");
    println!("├─");
    println!("│ AX = B");
    println!("│ A_00 = sum(x_i^2) = {:.4e}", sumx_squared);
    println!("│ A_01 = sum(x_i) = {:.4e}", sumx);
    println!("│ A_10 = sum(x_i) = {:.4e}", sumx);
    println!("│ A_11 = n = {}", xs.len());
    println!("│ B_0 = sum(x_i * y_i) = {:.4e}", sumxy);
    println!("│ B_1 = sum(y_i) = {:.4e}", sumy);
    println!("├─");
    println!("│ y = ax + b");
    println!("│ a = {:.4e}", x[0]);
    println!("│ b = {:.4e}", x[1]);
    println!("╰───────────────");

    Ok(x)
}

pub enum NonLinearRegressionType {
    Polynomial,
    Exponential,
    Logarithmic,
}

/// Least Squares Interpolation
pub fn nonlinear_regression(regression_type: NonLinearRegressionType) {
    todo!("Ch.4 p.84")
}
