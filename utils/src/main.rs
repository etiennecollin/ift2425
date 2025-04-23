#![allow(unused_variables, dead_code)]

use std::f64::consts::{FRAC_PI_4, PI};

use nalgebra::{dmatrix, dvector};
#[allow(unused_imports)]
use utils::{
    derivation::*, error_propagation::*, integration::*, interpolation::*, linear_systems::*,
    norms::*, root_search::*, utils::*,
};

const ITERATIONS_MAX: usize = 100;
const X_TOLERANCE: f64 = 1e-5;
const F_X_TOLERANCE: f64 = 1e-5;

fn main() {
    // error_analysis();
    // root_search();
    // linear_systems();
    // norms();
    // non_linear_systems();
    // interpolation();
    // derivative();
    integration();
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

    let error = 3.505e-3;
    let r = compute_sig_position(error);
    println!("Sig fig position: {}", r);
}

fn root_search() {
    // Bissection method
    // let f: FuncSingle = |x: f64| x.ln() - (-x).exp();
    // let interval = (1.0, 2.0);
    // let _ = get_bissection_iterations(interval, X_TOLERANCE).unwrap();
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
    // let _ = get_fixed_point_iterations_mvt(f, g, interval, max_error).unwrap();
    // let _ = get_fixed_point_x_tolerance_fit(f, g, interval, max_error).unwrap();
    // let _ = fixed_point(g, x_initial, X_TOLERANCE, F_X_TOLERANCE, ITERATIONS_MAX).unwrap();
}

fn linear_systems() {
    let a = dmatrix![
        1.424, 2.083;
        2.083, 4.000;
    ];
    let b = dvector![-1.402, -1.389];

    // Other methods
    // let _ = solve_gauss(&a, &b, true).unwrap();
    // let _ = solve_cramer(&a, &b).unwrap();

    // Solve with lu unpacked
    // let (l, u) = no_piv_lu(&a).unwrap();
    // let _ = solve_lu(&l, &u, &b);

    // Solve with lu packed
    // let lu = no_piv_lu_packed(&a).unwrap();
    // let (l, u) = lu_unpack(&lu);
    // let _ = solve_lu(&l, &u, &b);

    // Solve with lu packed mutable
    // let _ = no_piv_lu_packed_mut(&mut a).unwrap();
    // let (l, u) = lu_unpack(&a);
    // let _ = solve_lu(&l, &u, &b);

    // Solve LU with partial pivoting
    // let _ = partial_piv_lu(&a, &b).unwrap();

    // Solve LU with full pivoting
    // let _ = full_piv_lu(&a, &b).unwrap();

    // Least squares
    // let a = dmatrix![
    //     1.0, 1.0, 1.0, 1.0;
    //     8.0,  4.0, 2.0, 1.0;
    //     27.0, 9.0, 3.0, 1.0;
    //     64.0, 16.0, 4.0, 1.0;
    //     125.0, 25.0, 5.0, 1.0;
    //     216.0, 36.0, 6.0, 1.0;
    // ];
    // let b = dvector![5.04, 8.12, 10.64, 13.18, 16.20, 20.04];
    // let _ = least_squares(&a, &b).unwrap();

    // Iterative solve
    // let a = dmatrix![
    //     4.0, 2.0;
    //     2.0, 3.0;
    // ];
    // let a_inv = a.clone().try_inverse().unwrap();
    // let b = dvector![6.0, 5.0];
    // let _ = iterative_solve(&a, &a_inv, &b, 1000);
    // println!("Solution: {}", x);

    // Iterative solvers
    // let _ = jacobi_solver(&a, &b, X_TOLERANCE, ITERATIONS_MAX).unwrap();
    // let _ = gauss_seidel_solver(&a, &b, X_TOLERANCE, ITERATIONS_MAX).unwrap();
    // let _ = relaxation_solver(&a, &b, 0.95, X_TOLERANCE, ITERATIONS_MAX).unwrap();
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
    let x_initial = dvector![0.0, 0.0];

    // #[rustfmt::skip]
    // let system: Vec<FuncMulti> = vec![
    //     |p| p[1].cos()/3.0 + 1.0/6.0,
    //     |p| (0.5*p[0].powi(2) + 0.5).sqrt(),
    // ];

    // Newton's method
    // let _ = newton_system(
    //     &system,
    //     &x_initial,
    //     X_TOLERANCE,
    //     F_X_TOLERANCE,
    //     ITERATIONS_MAX,
    // )
    // .unwrap();

    #[rustfmt::skip]
    let system: Vec<FuncMulti> = vec![
        |p| -(4.0 - p[1].powi(2)).sqrt(),
        |p| 1.0 - p[0].exp(),
    ];

    // Fixed-point method
    // let _ = fixed_point_system(
    //     &system,
    //     &x_initial,
    //     X_TOLERANCE,
    //     F_X_TOLERANCE,
    //     ITERATIONS_MAX,
    // )
    // .unwrap();
}

fn interpolation() {
    let xs = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30];
    let fs = [1.0, 1.02470, 1.04881, 1.07238, 1.09544, 1.11803, 1.14017];

    let x = 1.28;
    let h = 0.05;
    let degree = 3;

    // let _ = lagrange(degree, x, &xs, &fs);

    // With f = sin(x) -> f^(degree+1) = -cos(x)
    // let derivative: FuncSingle = |x| -x.cos();
    // let _ = lagrange_polynomial_error_range(degree, derivative, x, &xs).unwrap();

    // Compute delta^k f_i
    // for d in 0..fs.len() {
    //     for i in 0..(fs.len() - d) {
    //         println!("delta: {}, i: {}, res: {}", d, i, delta_f_i(d, i, &fs));
    //     }
    // }

    // let _ = newton_gregory_forward(degree, x, h, &xs, &fs, true).unwrap();
    // let table = finite_forward_diff_table(degree, &fs).unwrap();
    // println!(
    //     "{}",
    //     finite_forward_diff_table_string(&xs[..=degree], &table).unwrap()
    // );

    // let xs = [0.0, 1.0, 2.0, 3.0];
    // let fs = [7.39, 4.66, 1.79, 1.0];
    // let degree = 3;
    // let table = finite_centered_diff_table(degree, &xs, &fs).unwrap();
    // println!("{}", finite_centered_diff_table_string(&table).unwrap());

    // let xs = [1.0, 0.54, -0.416, -0.99];
    // let fs = [2.0, 1.539, 0.582, 0.0];
    // let vector = linear_regression(&xs, &fs).unwrap();
    // let x: f64 = 1.0;
    // let x_true = 4.66;
    // let lin_f = (vector[0] * x.cos() + vector[1]).exp();
    // let error = (lin_f - x_true).abs();
    // println!("f({x}) = {lin_f}");
    // println!("Error: {error}");
    // println!("Sig fig Position: {}", compute_sig_position(error));
}

fn derivative() {
    // let f: FuncSingle = |x: f64| match x {
    //     2.3 => 0.34718,
    //     2.4 => 0.31729,
    //     2.5 => 0.28587,
    //     2.6 => 0.25337,
    //     2.7 => 0.22008,
    //     _ => panic!("Invalid x value, x = {}", x),
    // };
    // let x = 2.5;
    // let h_init = 0.1;
    // let level = 1;
    // let order = 1;
    // let _ = richardson(f, x, h_init, level, order, false).unwrap();

    // let xs = [1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5];
    // let fs = [3.669, 4.482, 5.474, 6.686, 8.166, 9.974, 12.182];
    // let derivative: FuncSingle = |x| x.exp();
    // let x_index = 2;
    // let h = 0.2;
    // let degree = 1;
    // let _ = newton_gregory_forward_derivative(degree, x_index, h, &xs, &fs).unwrap();
    // let _ = newton_gregory_forward_derivative_error_estimate(x_index, &xs, h, degree, derivative);
}

fn integration() {
    // let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    // let fs = [0.0, 1.0, 8.0, 27.0, 64.0, 125.0, 216.0];
    // let h = 1.0;
    // let degree = 4;
    // let range = (0, 6);
    // let _ = newton_cotes(degree, &xs, &fs, range).unwrap();

    let f: FuncSingle = |x| x.powi(2).exp();
    let h = 0.208139;
    let range = (0.0, 2f64.ln().sqrt());
    let method = QuadratureMethod::Simpson13;
    let _ = composite_quadrature(f, h, range, method).unwrap();

    // Max value of the second derivative in the `range`.
    // let target_error = 5e-4;
    // let error = composite_quadrature_trapezoidal_error(f, h, range, target_error).unwrap();
    // let r = compute_sig_position(error);
    // println!("Sig fig position: {}", r);

    // let f: FuncSingle = |x| match x {
    //     0.0 => 0.0,
    //     0.2 => 0.199,
    //     0.4 => 0.389,
    //     0.6 => 0.565,
    //     0.8 => 0.717,
    //     _ => panic!("Invalid x value, x = {}", x),
    // };
    // let h_init = 0.2;
    // let range = (0.0, 0.8);
    // let levels = 2;
    // let _ = romberg(f, h_init, range, levels).unwrap();

    // let f = |x| x.sin();
    // let degree = 2;
    // let range = (0.0, PI);
    // let _ = gaussian_quadrature(f, degree, range).unwrap();
}
