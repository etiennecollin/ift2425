use utils::{
    bissection::bissection,
    cse::cse,
    utils::{table_formatter, FuncMulti, FuncSingle},
};

fn main() {
    // // Error calculations
    // let epsilon = 1e-6;
    // let point = vec![549.12, 1327.5, 10250.65];
    // let error = vec![0.5e-2, 0.5e-1, 0.5e-2];
    // let func: FuncMulti = |point| {
    //     let g = |x: f64| x.powi(2) - 2.0 * x.sqrt();
    //     g((point[0] + point[2]) / point[1])
    // };
    // let _ = cse(func, &point, &error, epsilon);

    // Bissection method
    let func: FuncSingle = |x: f64| -3.0 - 3.0 * x + x.powi(2) + x.powi(3);
    let interval = (1.0, 2.0);
    let x_tolerance = 5e-4;
    let f_x_tolerance = 1e-8;
    let _ = bissection(func, interval, x_tolerance, f_x_tolerance, None);
}
