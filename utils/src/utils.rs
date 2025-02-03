use tabled::{
    builder::Builder,
    settings::{style::HorizontalLine, Panel, Style, Theme},
};

/// A f64 representing the perturbation used to compute the gradient.
const EPSILON: f64 = 1e-6;

/// This type defines a function that takes a slice of f64 and returns an f64.
/// Each element of the slice represents a variable in the function.
pub type FuncMulti = fn(&[f64]) -> f64;

/// This type defines a function that takes a single f64 and returns an f64.
pub type FuncSingle = fn(f64) -> f64;

// Function to compute the gradient of an arbitrary function at a given point.
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

// Function to compute the gradient of an arbitrary function at a given point.
pub fn compute_derivative(func: FuncSingle, x: f64) -> f64 {
    (func(x + EPSILON) - func(x - EPSILON)) / (2.0 * EPSILON)
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
