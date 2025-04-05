use nalgebra::{DMatrix, DVector};
use tabled::{
    builder::Builder,
    settings::{Panel, Style, Theme, style::HorizontalLine},
};

/// A f64 representing the perturbation used to compute the gradient.
pub const EPSILON: f64 = 1e-6;

/// This type defines a function that takes a slice of f64 and returns an f64.
/// Each element of the slice represents a variable in the function.
pub type FuncMulti = fn(&[f64]) -> f64;

/// This type defines a function that takes a single f64 and returns an f64.
pub type FuncSingle = fn(f64) -> f64;

/// Computes the gradient of an arbitrary function at a given point.
pub fn compute_gradient(func: FuncMulti, point: &[f64]) -> Vec<f64> {
    let mut gradient = Vec::with_capacity(point.len());

    for (i, &x) in point.iter().enumerate() {
        // Copy the point to avoid mutating the original
        let mut point_plus = point.to_vec();
        let mut point_minus = point.to_vec();

        // Perturb the i-th dimension
        point_plus[i] = x + EPSILON;
        point_minus[i] = x - EPSILON;

        // Compute the finite difference approximation of the derivative
        let derivative = (func(&point_plus) - func(&point_minus)) / (2.0 * EPSILON);

        // Store the derivative for this variable in the gradient vector
        gradient.push(derivative);
    }

    gradient
}

/// Computes the derivative of an arbitrary function at a given point.
pub fn compute_derivative(func: FuncSingle, x: f64) -> f64 {
    (func(x + EPSILON) - func(x - EPSILON)) / (2.0 * EPSILON)
}

// /// Computes the n-th derivative of an arbitrary function at a given point.
// ///
// /// WARNING: There are precision issues with n > 2
// pub fn compute_nth_derivative(func: FuncSingle, x: f64, n: u64) -> f64 {
//     let sum = (0..=n).fold(0.0, |acc, k| {
//         let sign = if (k + n) % 2 == 0 { 1.0 } else { -1.0 };
//         let binom = binomial(n, k) as f64;
//         let fx = func(x + (k as f64) * EPSILON);
//         acc + sign * binom * fx
//     });
//     sum / EPSILON.powi(n as i32)
// }

/// Computes at a given point the value of every function in a system of equations.
pub fn function_vec(system: &[FuncMulti], point: &[f64]) -> DVector<f64> {
    let data: Vec<_> = system.iter().map(|f| f(point)).collect();
    DVector::from_vec(data)
}

/// Computes at a given point the gradient of every function in a system of equations.
pub fn gradient_mat(system: &[FuncMulti], point: &[f64]) -> DMatrix<f64> {
    let data: Vec<_> = system
        .iter()
        .flat_map(|f| compute_gradient(*f, point))
        .collect();
    DMatrix::from_vec(system.len(), point.len(), data).transpose()
}

/// Defines the possible errors that can occur when formatting a table.
#[derive(Debug)]
pub enum TableError {
    EmptyTable,
    HeaderRowsLengthMismatch,
}

impl From<TableError> for &'static str {
    fn from(val: TableError) -> Self {
        match val {
            TableError::EmptyTable => "The table is empty.",
            TableError::HeaderRowsLengthMismatch => {
                "The header and rows have different number of columns."
            }
        }
    }
}

/// This function formats a table with a header and contents.
///
///
/// # Arguments
///
/// - `info`: A vector of strings representing the information to put at the top of the table.
/// - `header`: A vector of strings representing the header of the table.
/// - `rows`: A vector of strings representing the contents of the table.
///     The columns of the rows should be comma-separated.
///
/// # Returns
///
/// - A string representing the formatted table.
///
/// # Errors
///
/// - If the length of the header and rows does not match.
/// - If the table is empty.
pub fn table_formatter<T>(
    info: Vec<String>,
    header: Vec<String>,
    rows: Vec<Vec<T>>,
) -> Result<String, TableError>
where
    T: Into<String>,
{
    // Check if the header and rows have the same number of columns
    if (header.is_empty() && !rows.is_empty())
        || (!header.is_empty() && rows.is_empty())
        || rows.iter().any(|row| row.len() != header.len())
    {
        return Err(TableError::HeaderRowsLengthMismatch);
    }

    // Check if the table is empty
    if info.is_empty() && (header.is_empty() && rows.is_empty()) {
        return Err(TableError::EmptyTable);
    }

    // Create the table and prepare rows
    let mut builder = Builder::default();
    builder.push_record(header);
    rows.into_iter().for_each(|row| {
        builder.push_record(row);
    });
    let mut table = builder.build();

    // Prepare the style of the table
    let mut theme = Theme::from_style(Style::modern_rounded());
    let hline = HorizontalLine::inherit(Style::modern_rounded());
    theme.remove_horizontal_lines();

    // Add a horizontal line after the header and the info
    theme.insert_horizontal_line(1, hline);
    if !info.is_empty() {
        theme.insert_horizontal_line(2, hline);
    }

    // Build the table
    let info_str = info.join("\n");
    table.with(theme).with(Panel::header(info_str));
    Ok(table.to_string())
}

/// Computes the factorial of a number.
///
/// # Errors
///
/// - If the number is greater than 170, an error is returned as this is the largest
///   factorial that can be represented as a f64 without overflow.
pub fn factorial(num: u8) -> Result<f64, &'static str> {
    if num > 170 {
        return Err("Cannot compute the factorial of a number greater than 170.");
    }

    let mut result = 1.0;
    for x in 1..=num {
        result *= x as f64;
    }
    Ok(result)
}

/// Computes the binomial coefficient "n choose k".
pub fn choose_float(n: f64, k: u64) -> Result<f64, &'static str> {
    if k as f64 > n {
        return Ok(0.0); // n choose k is 0 when k > n
    }
    if k == 0 || k as f64 == n {
        return Ok(1.0); // n choose 0 and n choose n are always 1
    }

    let product = (0..k).fold(1.0, |acc, i| acc * (n - i as f64));
    let result = product / factorial(k as u8)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test factorial of small numbers
    #[test]
    fn test_factorial_small_numbers() {
        assert_eq!(factorial(0), Ok(1.0)); // 0! = 1
        assert_eq!(factorial(1), Ok(1.0)); // 1! = 1
        assert_eq!(factorial(2), Ok(2.0)); // 2! = 2
        assert_eq!(factorial(3), Ok(6.0)); // 3! = 6
        assert_eq!(factorial(4), Ok(24.0)); // 4! = 24
        assert_eq!(factorial(5), Ok(120.0)); // 5! = 120
    }

    // Test factorial for the largest allowed number (170)
    #[test]
    fn test_factorial_170() {
        assert!(factorial(170).is_ok()); // factorial(170) should not return an error
    }

    // Test that factorial of numbers greater than 170 returns an error
    #[test]
    fn test_factorial_greater_than_170() {
        assert!(factorial(171).is_err(),);
        assert!(factorial(200).is_err());
    }

    // Test that factorial returns correct large values
    #[test]
    fn test_factorial_large_numbers() {
        // Here you can validate the factorial for some smaller large numbers if needed
        // For instance, let's test factorial of 10, 15, 20
        assert_eq!(factorial(10), Ok(3_628_800.0)); // 10! = 3,628,800
        assert_eq!(factorial(15), Ok(1_307_674_368_000.0)); // 15! = 1,307,674,368,000
        assert_eq!(factorial(20), Ok(2_432_902_008_176_640_000.0)); // 20! = 2,432,902,008,176,640,000
    }

    // Test that n choose 0 and n choose n are always 1
    #[test]
    fn test_n_choose_0_and_n_choose_n() {
        assert_eq!(choose_float(5.0, 0), Ok(1.0)); // 5 choose 0 = 1
        assert_eq!(choose_float(5.0, 5), Ok(1.0)); // 5 choose 5 = 1
        assert_eq!(choose_float(10.0, 0), Ok(1.0)); // 10 choose 0 = 1
        assert_eq!(choose_float(10.0, 10), Ok(1.0)); // 10 choose 10 = 1
    }

    // Test that n choose k returns 0 when k > n
    #[test]
    fn test_choose_greater_than_n() {
        assert_eq!(choose_float(5.0, 6), Ok(0.0)); // 5 choose 6 = 0
        assert_eq!(choose_float(10.0, 11), Ok(0.0)); // 10 choose 11 = 0
    }

    // Test symmetry: n choose k == n choose (n - k)
    #[test]
    fn test_choose_symmetry() {
        assert_eq!(choose_float(5.0, 2), choose_float(5.0, 3)); // 5 choose 2 == 5 choose 3
        assert_eq!(choose_float(10.0, 4), choose_float(10.0, 6)); // 10 choose 4 == 10 choose 6
    }

    // Test some common known values
    #[test]
    fn test_choose_common_values() {
        assert_eq!(choose_float(5.0, 2), Ok(10.0)); // 5 choose 2 = 10
        assert_eq!(choose_float(6.0, 3), Ok(20.0)); // 6 choose 3 = 20
        assert_eq!(choose_float(7.0, 4), Ok(35.0)); // 7 choose 4 = 35
    }

    // Test large values
    #[test]
    fn test_choose_large_values() {
        assert_eq!(choose_float(30.0, 15), Ok(155117520.0)); // 30 choose 15
        assert_eq!(choose_float(20.0, 10), Ok(184756.0)); // 20 choose 10
    }
}
