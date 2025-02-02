use std::fmt::Write;

/// This type defines a function that takes a slice of f64 and returns an f64.
/// Each element of the slice represents a variable in the function.
pub type FuncMulti = fn(&[f64]) -> f64;

/// This type defines a function that takes a single f64 and returns an f64.
pub type FuncSingle = fn(f64) -> f64;

// Function to compute the gradient of an arbitrary function at a given point.
pub fn compute_gradient(func: FuncMulti, point: &[f64], epsilon: f64) -> Vec<f64> {
    let mut gradient = Vec::with_capacity(point.len());

    for (i, &x) in point.iter().enumerate() {
        // Copy the point to avoid mutating the original
        let mut point_plus = point.to_vec();
        let mut point_minus = point.to_vec();

        // Perturb the i-th dimension
        point_plus[i] = x + epsilon;
        point_minus[i] = x - epsilon;

        // Compute the finite difference approximation of the derivative
        let derivative = (func(&point_plus) - func(&point_minus)) / (2.0 * epsilon);

        // Store the derivative for this variable in the gradient vector
        gradient.push(derivative);
    }

    gradient
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
/// # Panics
///
/// - If the number of columns in the header is different from the number of columns in the rows.
///
/// # Quirks
///
/// - This works well enough for me.
/// - The table width only depends on the `header` and `rows`, not the `info`.
/// - The behaviour is weird and might [panic] if the `info`, `header` or `rows` is empty.
pub fn table_formatter(info: Vec<String>, header: String, rows: Vec<String>) -> String {
    // These variables will store the maximum width of each column
    let column_number = header.split(",").count();

    // Check that the number of columns in the header is the same as the number of columns in the rows
    rows.iter().for_each(|row| {
        assert_eq!(row.split(",").count(), column_number);
    });

    // Compute the width of each column
    let mut widths: Vec<usize> = vec![0; column_number];
    rows.iter().for_each(|row| {
        row.split(",")
            .map(|col| col.trim().len())
            .enumerate()
            .for_each(|(i, width)| widths[i] = widths[i].max(width))
    });
    header
        .split(",")
        .map(|col| col.trim().len())
        .enumerate()
        .for_each(|(i, width)| widths[i] = widths[i].max(width));

    // Compute the width of the widest formatted row
    // The width of the row is:
    // - Content: The width of each column's content
    // - Separators: columns + 1
    // - Spaces: columns * 2
    let row_max_width = widths.iter().sum::<usize>() + 3 * column_number + 1;

    // Correctly format the header
    let header_string = header.split(",").map(|col| col.trim()).enumerate().fold(
        String::with_capacity(row_max_width),
        |mut output, (i, col)| {
            let width = widths[i];
            if i == 0 {
                let _ = write!(output, "│ {:<width$} │", col);
            } else {
                let _ = write!(output, " {:<width$} │", col);
            }
            output
        },
    );

    // Correctly format each row
    let rows_string = rows.iter().enumerate().fold(
        String::with_capacity(row_max_width * rows.len()),
        |mut output, (i, row)| {
            // Split the row into columns
            // Trim each column
            // Write the formatted column to the output
            row.split(",")
                .map(|col| col.trim())
                .enumerate()
                .for_each(|(j, col)| {
                    // Get the max width of the column
                    let width = widths[j];

                    // Format the column
                    if j == 0 {
                        // If it is the first column, add a separator
                        let _ = write!(output, "│ {:<width$} │", col);
                    } else if j == column_number - 1 && i != rows.len() - 1 {
                        // If it is the last column, add a newline except for the last row
                        let _ = writeln!(output, " {:<width$} │", col);
                    } else {
                        // If it is a middle column or the last column of the last row, add a separator
                        let _ = write!(output, " {:<width$} │", col);
                    }
                });

            output
        },
    );

    // Get the maximum width of the table
    let table_width: usize = header_string.len() - 2 * column_number - 4;

    // Get the width of the header (there are 2 delimiters)
    let info_width = table_width - 2;

    // Format the info
    let info_string = info
        .iter()
        .enumerate()
        .fold(String::new(), |mut output, (i, inf)| {
            if i == info.len() - 1 {
                let _ = write!(output, "│ {:<info_width$} │", inf,);
            } else {
                let _ = writeln!(output, "│ {:<info_width$} │", inf,);
            }
            output
        });

    // Create the table
    let table: Vec<String> = vec![
        format!("╭{}╮", "─".repeat(table_width)),
        info_string,
        format!("├{}┤", "─".repeat(table_width)),
        header_string,
        format!("├{}┤", "─".repeat(table_width)),
        rows_string,
        format!("╰{}╯", "─".repeat(table_width)),
    ];

    // Concatenate the table into a single string
    table.join("\n")
}
