#![allow(unused_variables, dead_code)]

use nalgebra::{dmatrix, dvector};
#[allow(unused_imports)]
use utils::{cse::*, norms::*, root_search::*, systems::*, utils::*};

const ITERATIONS_MAX: usize = 100;
const X_TOLERANCE: f64 = 1e-5;
const F_X_TOLERANCE: f64 = 1e-5;

fn main() {
    error_analysis();
    root_search();
    linear_systems();
    norms();
    non_linear_systems();
}

fn error_analysis() {
    let func: FuncMulti = |p| {
        let g = |x: f64| x.powi(2) - 2.0 * x.sqrt();
        g((p[0] + p[2]) / p[1])
    };

    let point = [549.12, 1327.5, 10250.65];
    let error = [0.5e-2, 0.5e-1, 0.5e-2];

    // let _ = cse_taylor(func, &point, &error);
    // let _ = cse_fork(func, &point, &error);
}

fn root_search() {
    // Bissection method
    // let f: FuncSingle = |x: f64| x.ln() - (-x).exp();
    // let interval = (1.0, 2.0);
    // let _ = get_bissection_iterations(interval, X_TOLERANCE);
    // let _ = bissection(f, interval, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();

    // Linear interpolation method
    // let _ =
    //     linear_interpolation(f, interval, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();

    // Secant method
    // let _ = secant(f, interval, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();

    // Newton's method
    // To converge, |f'(x_initial)| > 0
    // let x_initial = 1.5;
    // let f: FuncSingle = |x: f64| x.ln() - (-x).exp();
    // let _ = newton(f, x_initial, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();

    // Fixed-point method
    // To converge, |g'(x_initial)| < 1
    // let f: FuncSingle = |x: f64| x.ln() - (-x).exp();
    // let g: FuncSingle = |x: f64| (-x).exp().exp();
    // let g: FuncSingle = |x: f64| -x.ln().ln();
    // let x_initial = 1.5;
    // let interval = (1.0, 2.0);
    // let max_error = 1e-5;
    // let _ = get_fixed_point_iterations_mvt(f, g, interval, max_error);
    // let _ = get_fixed_point_x_tolerance_fit(f, g, interval, max_error);
    // let _ = fixed_point(g, x_initial, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();
}

fn linear_systems() {
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
    // let x = jacobi_solver(&a, &b, X_TOLERANCE, ITERATIONS_MAX);
    // let x = gauss_seidel_solver(&a, &b, X_TOLERANCE, ITERATIONS_MAX);
    // let x = relaxation_solver(&a, &b, 0.95, X_TOLERANCE, ITERATIONS_MAX);
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
    //     X_TOLERANCE,
    //     F_X_TOLERANCE,
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
    //     X_TOLERANCE,
    //     F_X_TOLERANCE,
    //     ITERATIONS_MAX,
    // );
}
