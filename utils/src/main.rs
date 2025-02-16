use nalgebra::{dmatrix, dvector, DMatrix};
use utils::{cse::*, linear_systems::*, root_search::*, utils::*};

// fn main() {
//     const ITERATIONS_MAX: usize = 100;
//
//     // -------------------------------------------------------------------------
//     // Error calculations
//     // -------------------------------------------------------------------------
//     let func: FuncMulti = |point| {
//         let g = |x: f64| x.powi(2) - 2.0 * x.sqrt();
//         g((point[0] + point[2]) / point[1])
//     };
//     let point = vec![549.12, 1327.5, 10250.65];
//     let error = vec![0.5e-2, 0.5e-1, 0.5e-2];
//     let _ = cse(func, &point, &error);
//
//     // -------------------------------------------------------------------------
//     // Root search
//     // -------------------------------------------------------------------------
//     let func: FuncSingle = |x: f64| -3.0 - 3.0 * x + x.powi(2) + x.powi(3);
//     let x_initial = 2.0; // Only for Newton's method
//     let interval = (1.0, 2.0); // Only for bissection, linear interpolation and secant
//     let x_tolerance = 5e-4;
//     let f_x_tolerance = 1e-5;
//
//     // Bissection method
//     let _ = bissection(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
//
//     // Linear interpolation method
//     let _ =
//         linear_interpolation(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
//
//     // Secant method
//     let _ = secant(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
//
//     // Newton's method
//     let _ = newton(func, x_initial, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
//
//     // Fixed-point method
//     // To converge, |g'(x_initial)| < 1
//     let g: FuncSingle = |x: f64| (-x.powi(3) - x.powi(2) + 3.0) / -3.0;
//     let x_initial = 0.3;
//     let _ = fixed_point(g, x_initial, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
// }

fn main() {
    let a = dmatrix![
        2.0, -1.0, -1.0;
        0.0,  -4.0, 2.0;
        6.0, -3.0, 0.0;
    ];
    let b = dvector![0.0, -2.0, 3.0];
    // let (l, u) = no_piv_lu(&a);
    // println!("L = {}", l);
    // println!("U = {}", u);

    let lu = no_piv_lu_packed(&a);
    let (l, u) = lu_unpack(&lu);
    println!("LU = {}", a);
    println!("L = {}", l);
    println!("U = {}", u);

    // solve_lu(&l, &u, &b);

    // // Solve with lu
    // let (l, u) = no_piv_lu(&a).unwrap();
    // let _x = solve_lu(&l, &u, &b);
    //
    // // LU with partial pivoting
    // let x = partial_piv_lu(&a, &b);
    //
    // LU with full pivoting
    // let x = full_piv_lu(&a, &b);

    // Other methods
    // let x = solve_gauss(&a, &b, true);
    // let x = solve_cramer(&a, &b);
}
