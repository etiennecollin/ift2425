#![allow(unused_variables)]

use nalgebra::{dmatrix, dvector};
#[allow(unused_imports)]
use utils::{cse::*, norms::*, root_search::*, systems::*, utils::*};

#[allow(dead_code)]
const ITERATIONS_MAX: usize = 100;

fn main() {
    error_analysis();
    root_search();
    linear_systems();
    norms();
    non_linear_systems();
}

fn error_analysis() {
    let func: FuncMulti = |p| {
        // let g = |x: f64| x.powi(2) - 2.0 * x.sqrt();
        // g((p[0] + p[2]) / p[1])
        p[0].powf(p[1]) / p[1].powf(p[0])
    };

    // let point = [549.12, 1327.5, 10250.65];
    // let error = [0.5e-2, 0.5e-1, 0.5e-2];
    let point = [2, 3].map(|x| x as f64);
    let error = [0.2, 0.3];

    // let _ = cse_taylor(func, &point, &error);
    // let _ = cse_fork(func, &point, &error);
}

fn root_search() {
    let func: FuncSingle = |x: f64| -3.0 - 3.0 * x + x.powi(2) + x.powi(3);
    let x_tolerance = 5e-4;
    let f_x_tolerance = 1e-5;
    let x_initial = 2.0; // Only for Newton's method
    let interval = (1.0, 2.0); // Only for bissection, linear interpolation and secant

    // Bissection method
    // let _ = bissection(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();

    // Linear interpolation method
    // let _ =
    //     linear_interpolation(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();

    // Secant method
    // let _ = secant(func, interval, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();

    // Newton's method
    // let _ = newton(func, x_initial, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();

    // Fixed-point method
    // To converge, |g'(x_initial)| < 1
    // let g: FuncSingle = |x: f64| (-x.powi(3) - x.powi(2) + 3.0) / -3.0;
    // let x_initial = 0.3;
    // let _ = fixed_point(g, x_initial, x_tolerance, f_x_tolerance, ITERATIONS_MAX).unwrap();
}

fn linear_systems() {
    let x_tolerance = 1e-5;
    let a = dmatrix![
        2.0, -1.0, -1.0;
        0.0,  -4.0, 2.0;
        6.0, -3.0, 0.0;
    ];
    let b = dvector![0.0, -2.0, 3.0];

    // Other methods
    // let x = solve_gauss(&a, &b, true);
    // let x = solve_cramer(&a, &b);

    // Solve with lu unpacked
    // let (l, u) = no_piv_lu(&a);
    // let x = solve_lu(&l, &u, &b);

    // Solve with lu packed
    // let lu = no_piv_lu_packed(&a);
    // let (l, u) = lu_unpack(&lu);
    // let x = solve_lu(&l, &u, &b);

    // Solve with lu packed mutable
    // no_piv_lu_mut(&mut a);
    // let (l, u) = lu_unpack(&a);
    // let x = solve_lu(&l, &u, &b);

    // LU with partial pivoting
    // let x = partial_piv_lu(&a, &b);

    // LU with full pivoting
    // let x = full_piv_lu(&a, &b);

    // Least squares
    // let a = dmatrix![
    //     2.0, 1.0;
    //     1.0,  -2.0;
    //     -3.0, 2.0;
    // ];
    // let b = dvector![4.0, -3.0, 3.0];
    // let x = least_squares(&a, &b).unwrap();

    // Iterative solve
    // let a = dmatrix![
    //     4.0, 2.0;
    //     2.0, 3.0;
    // ];
    // let a_inv = a.clone().try_inverse().unwrap();
    // let b = dvector![6.0, 5.0];
    // let x = iterative_solve(&a, &a_inv, &b, 1000).unwrap();
    // println!("Solution: {}", x);

    // Iterative solvers
    // let x = jacobi_solver(&a, &b, x_tolerance, ITERATIONS_MAX);
    // let x = gauss_seidel_solver(&a, &b, x_tolerance, ITERATIONS_MAX);
    // let x = relaxation_solver(&a, &b, 0.95, x_tolerance, ITERATIONS_MAX);
}

fn norms() {
    // Vector norm
    // let norm = norm_vec(&b, 0);
    // println!("Norm: {}", norm);

    // Matrix norm
    // let norm = norm_mat(&mat, 3);
    // println!("Norm: {}", norm);

    // Condition number
    // let cond = condition_number(&a, 0).unwrap();
    // println!("Condition number: {}", cond);

    // Error propagation
    // let a = dmatrix![
    //     10, 7, 8, 7;
    //     7, 5, 6, 5;
    //     8, 6, 10, 9;
    //     7, 5, 9, 10;
    // ]
    // .map(|x| x as f64);
    // let a_star = dmatrix![
    //     10.0, 7.0, 8.1, 7.2;
    //     7.0, 5.01, 6.0, 5.02;
    //     8.0, 6.0, 10.0, 9.0;
    //     6.99, 4.99, 9.1, 10.0;
    // ];
    // let b = dvector![32, 23, 33, 31].map(|x| x as f64);
    // let b_star = dvector![32.1, 22.9, 33.1, 30.9];
    // let range = x_error_range_b(&a, &b, &b_star, 0).unwrap();
    // let upper = x_error_bound_a(&a, &a_star, 0).unwrap();
}

fn non_linear_systems() {
    let x_initial = dvector![-2.0, 1.0];
    let x_tolerance = 1e-5;
    let f_x_tolerance = 1e-5;

    // #[rustfmt::skip]
    // let system: Vec<FuncMulti> = vec![
    //     |p| p[0].powi(2) + p[1].powi(2) - 4.0,
    //     |p| p[0].exp() + p[1] - 1.0,
    // ];
    //
    // Newton's method
    // let x = newton_system(
    //     &system,
    //     &x_initial,
    //     x_tolerance,
    //     f_x_tolerance,
    //     ITERATIONS_MAX,
    // );

    // #[rustfmt::skip]
    // let system: Vec<FuncMulti> = vec![
    //     |p| -(4.0 - p[1].powi(2)).sqrt(),
    //     |p| 1.0 - p[0].exp(),
    // ];
    //
    // Fixed-point method
    // let x = fixed_point_system(
    //     &system,
    //     &x_initial,
    //     x_tolerance,
    //     f_x_tolerance,
    //     ITERATIONS_MAX,
    // );
}
