fn f(x: f64, y: &[f64]) -> f64 {
    // sum1 : y_i^(x) ln(y_i)
    let sum1 = y.iter().map(|&y_i| y_i.powf(x) * y_i.ln()).sum::<f64>();
    // sum2 : y_i^(x)
    let sum2 = y.iter().map(|&y_i| y_i.powf(x)).sum::<f64>();
    // sum3 : ln(y_i)
    let sum3 = y.iter().map(|&y_i| y_i.ln()).sum::<f64>();
    (sum1 / sum2) - (1.0 / x) - (sum3 / y.len() as f64)
}

fn derivative_f(x: f64, y: &[f64]) -> f64 {
    let h = 1e-5;
    (f(x + h, y) - f(x, y)) / h
}

pub fn newton(x: f64, y: &[f64]) -> f64 {
    let delta = 1e-5;
    let epsilon = 1e-6;

    let mut x1 = x;
    let mut x2 = 0.0;

    while (x2 - x1).abs() >= delta && f(x1, y).abs() >= epsilon && derivative_f(x1, y) != 0.0 {
        x2 = x1;
        x1 = x1 - f(x1, y) / derivative_f(x1, y);
    }

    x1
}
