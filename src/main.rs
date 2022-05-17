use std::{io::{Write, self, BufRead, BufReader}, collections::{VecDeque, HashSet}, ops::Index, f64::consts::PI, fs::File, path::Path, iter::Map};
use ndarray::{prelude::*, IndexLonger, AssignElem};
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
        Point6 {x: 0., y: 0., z: 0., r: 0., phi: 0., theta: 0.}
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

    fn from_spherical<T: Index<usize, Output = f64>>(data: &T) -> Point6 {
        let r = data[0];
        let phi = data[1];
        let theta = data[2];

        Point6 { x: r*theta.sin()*phi.cos(), 
                 y: r*theta.sin()*phi.sin(), 
                 z: r*theta.cos(), 
                 r, 
                 phi, 
                 theta }
    }
    // methods

    fn assert_angles(&mut self) {
        //phi [0, 2*PI]
        if self.phi < 0. { self.phi += 2.*PI}
        else if self.phi >2.*PI { self.phi -= 2.*PI  }

        //theta [0, PI]
        if self.theta < 0. { self.theta += PI}
        else if self.theta > PI { self.theta -= PI  }

    }
}

impl std::fmt::Display for Point6 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}",
                 self.x, self.y, self.z, self.r, self.phi, self.theta)
    }
}

type Point6Array = Array1<Point6>;

#[derive( Debug)]
struct Fuleren {
    positions: Point6Array,
    size: usize,
    E: f64,
}

impl Fuleren {
    // constructors
    fn new(size: usize) -> Fuleren {
        Fuleren { positions: Point6Array::from_elem(size, Point6::new()),
                  size,
                  E: 0. }
    }
    
    fn from_file(path: &str) -> Result<Fuleren, String>  {
        
        if let Ok(lines) = read_lines(path) {
            let iter = lines
                                                    .map(|line| line
                                                        .expect("wrong line")
                                                        .split_ascii_whitespace()
                                                        .map(|num_str| num_str.parse::<f64>().expect("error duting parsing"))
                                                        .collect::<Array1<f64>>())
                                                    .map(|data| Point6::from_cartesian(&data));
            let pos_array: Point6Array = iter.collect();
        Ok(Fuleren {size: pos_array.len(), E: 0.,
                positions: pos_array} )
        }
        else {
            Err("Error during reading from file".to_string())
        }

    }

    // methods
    fn randomize_on_sphere(&mut self, r: f64) {
        let phi_distr = rand::distributions::Uniform::new_inclusive(0., 2.*PI);
        let theta_distr = rand::distributions::Uniform::new_inclusive(0., PI);
        let mut rng = rand::thread_rng();

        self.positions.iter_mut()
                      .for_each(|point| 
                                point.assign_elem(Point6::from_spherical(&[r, 
                                                                            rng.sample(phi_distr), 
                                                                            rng.sample(theta_distr)]) ));
    }
}

impl std::fmt::Display for Fuleren {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut res = write!(f, "Fuleren with {} atoms, Energy: {:8.3}\n", self.size, self.E);
        res = write!(f, "{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}\t{:<10.5}\n", "x", "y", "z", "r", "phi", "theta");
        for point in self.positions.iter(){
            res = write!(f, "{}\n", *point);
        }
        res
    }
}



// ####################################
// ########### functions #############


fn check_angles(mut phi: f64, mut theta: f64) -> (f64, f64) {
    //phi [0, 2*PI]
    if phi < 0. { phi += 2.*PI}
    else if phi >2.*PI { phi -= 2.*PI  }

    //theta [0, PI]
    if theta < 0. { theta += PI}
    else if theta > PI { theta -= PI  }

    (phi, theta)
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename).expect("cannot read the file");
    Ok(io::BufReader::new(file).lines())
}

// ##################################

fn main() {
    let mut rng = rand::thread_rng();
    // let distr = rand::distributions::Uniform::new_inclusive(0, 3);

    let N = 10;

    // let point1 = Point6::from_cartesian(&[1.,2.,3.]);

    // let mut F = Fuleren::new(10);

    // println!("{}", F);

    // F.randomize_on_sphere(1.);
    // println!("{}", F);

    let F = Fuleren::from_file("data/atoms_test.dat").unwrap();
    println!("{}", F);

    
    
}

