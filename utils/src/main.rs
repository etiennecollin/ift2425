use utils::newton::newton;

fn main() {
    let x = 0.25;
    let y = &[0.11, 0.24, 0.27, 0.52, 1.13, 1.54, 1.71, 1.84, 1.92, 2.01];

    let res = newton(x, y);
    println!("f({}) = {}", x, res);
}
