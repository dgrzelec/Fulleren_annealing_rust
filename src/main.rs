use std::{io::{Write}, collections::{VecDeque, HashSet}, ops::Index};
use ndarray::{prelude::*, IndexLonger};
use rand::prelude::*;
use utilities::save_gnuplot2D;

use crate::utilities::get_file_buffer;

mod utilities;

//################# params ###################
const R0: f64 = 1.315;
const R1: f64 = 1.7;
const R2: f64 = 2.0;
const De: f64 = 6.325;
const S: f64 = 1.29;
const lambda: f64 = 1.5;
const del: f64 = 0.80469;
const a0: f64 = 0.011304;
const c0: f64 = 19.;
const d0: f64 = 2.5;
// ##############################
const N: usize = 10+2;
type MatrixInt = Array2<i32>;
type VectorInt = Array1<i32>;
type VectorFloat = Array1<f64>;

// ############# structs and implementations
#[derive( Debug, Clone)]
struct Point6 {
    x: f64,
    y: f64,
    z: f64,
    r: f64,
    phi: f64,
    theta: f64
}

impl Point6 {
    fn new() -> Point6 {
        // Point6 {}
        todo!()
    }

    fn from_cartesian<T: Index<usize, Output = f64>>(data: &T) -> Point6 {
        let xt:f64 = data[0];
        let yt = data[1];
        let zt = data[2];
        let rt = (xt.powi(2) + yt.powi(2) + zt.powi(2)).sqrt();
        Point6 { x: xt, 
                 y: yt, 
                 z: zt, 
                 r: rt, 
                 phi: (yt/xt).atan(), 
                 theta: (zt/rt).acos() }
    }
}
#[derive( Debug)]
struct Fuleren {
    positions: Array1<Point6>,

}

impl Fuleren {
    fn new() -> Fuleren {
        todo!()
    }
    
}

// ####################################
// ########### functions #############



// ##################################

fn main() {
    let mut rng = rand::thread_rng();
    // let distr = rand::distributions::Uniform::new_inclusive(0, 3);

    
    let mut Board = MatrixInt::from_elem((N,N), 5);

    let point1 = Point6::from_cartesian(&[1.,2.,3.]);

    println!("{:?}", point1);
}

