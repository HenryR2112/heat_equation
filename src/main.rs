extern crate nalgebra as na;

fn main() {
    let m = 10;
    let h = 1.0 / (m + 1) as f64;;
    let gamma = 4.0*h*h;
    let N = 20;
    let alpha = 1.5;
    //heat_equation(m, h as f64, gamma as f64, N, alpha)
    let (Adiaglower, Adiagmain, Adiagupper) = create_matrix_A_sparse(alpha, gamma, h, m);
    println!("adiaglower: {:?}", Adiaglower);
    println!("Adiagmain: {:?}", Adiagmain);
    println!("Adiagupper: {:?}", Adiagupper);
}

fn heat_equation(m: usize, h: f64, gamma: f64, N: usize, alpha: f64){
    let mut G = vec![0.0; m + 2];
    let mut T = vec![0.0; N + 1];

    for i in 0..=m {
        G[i] = i as f64 * h; // Calculate the value of G[i]
    }

// Populate T vector
    for j in 0..=N {
        T[j] = j as f64 * gamma; // Calculate the value of T[j]
    }

    //println!("G: {:?}", G);
    //println!("T: {:?}", T);
}

fn f1(x: f64, t: f64) -> f64 {
    let condition = t == 0.0;
    let value = if condition {
        2.5 * (f64::abs((f64::sin(f64::sin(x * std::f64::consts::PI * 3.0) * 5.0 * std::f64::consts::PI))) * x.powi(2))
    } else {
        0.0
    };
    value
}

fn f2(x: f64, t: f64) -> f64 {
    let condition = x >= 0.4 && x <= 0.8;
    let value = f64::abs(f64::sin(f64::sin(x * std::f64::consts::PI * 3.0) * 5.0 * std::f64::consts::PI)) * (10.0 * t + 1.0)
        + if condition { f64::abs(x) } else { 0.0 };
    value
}

fn create_matrix_A_sparse(alpha: f64, gamma: f64, h: f64, m: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let c = alpha * gamma / (h * h);
    let mut Ddiaglower = vec![1.0; m + 1];
    let mut DdiagMain = vec![-2.0; m + 2];
    DdiagMain[0] = 0.0;
    DdiagMain[m + 1] = 0.0;
    let mut Ddiagupper = vec![1.0; m + 1];

    let Adiaglower: Vec<f64> = Ddiaglower.iter().map(|&x| -c * x).collect();
    let Adiagmain: Vec<f64> = DdiagMain.iter().map(|&x| 1.0 - c * x).collect();
    let Adiagupper: Vec<f64> = Ddiagupper.iter().map(|&x| -c * x).collect();

    (Adiaglower, Adiagmain, Adiagupper)
}

fn tridiag_LU(lowerdiag: Vec<f64>, maindiag: Vec<f64>, upperdiag: Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut Ldiag = lowerdiag;
    let mut Udiag1 = maindiag;
    let mut Udiag2 = upperdiag;
    let n = Udiag1.len();
    for i in 1..=n{
       let mut factor = Ldiag[i] / Udiag1[i];
        Udiag1[i+1] = Udiag1[i+1] - factor*Udiag2[i];
        Ldiag[i] = factor;
    }
    (Ldiag,Udiag1, Udiag2)
}

