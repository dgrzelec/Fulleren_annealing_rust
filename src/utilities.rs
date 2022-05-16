use ndarray::prelude::*;


use ndarray::{Array1, Array2};
use std::fmt::{Display, Debug};
/// files
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn get_file_buffer(path: &str) -> BufWriter<File>{
    let f = File::create(path).expect("unable to create file");
    BufWriter::new(f)
}

/// saves given 1D ndarray to file named in path argument; Produces Gnuplot ready files
pub fn save_gnuplot1D<T: Display>(data: &Array1<T>, path: &str){
    
    let mut f = get_file_buffer(path);    

    let i_width = std::cmp::max(5,data.len().to_string().len()+2);
    let data_width = std::cmp::max(8, data[0].to_string().len());

    for i in 0..data.len(){
        write!(f, "{:<i_width$} {:<data_width$}\n", i, data[i]).expect("nie udało sie zapisac");
    }
    write!(f, "\n").expect("nie udało sie zapisac");
}


/// saves given 2D ndarray to file named in path argument; Produces Gnuplot ready files
pub fn save_gnuplot2D<T: Display>(data: &Array2<T>, path: &str){
    
    let mut f = get_file_buffer(path);    


    // calculates width of given variable in string to save;
    let i_width = std::cmp::max(5,data.shape()[0].to_string().len()+2);
    let j_width = std::cmp::max(5,data.shape()[1].to_string().len()+2);
    // notice assumption that single data element is enough to determine the width
    let data_width = std::cmp::max(8, data[[0,0]].to_string().len());

    for i in 0..data.shape()[0]{
        for j in 0..data.shape()[1]{
            
            write!(f, "{:<i_width$} {:<j_width$} {:<data_width$}\n", i, j, data[[i,j]]).expect("nie udało sie zapisac");
        }
        write!(f, "\n").expect("nie udało sie zapisac");
    }
    write!(f, "\n").expect("nie udało sie zapisac");
}


fn main() {
    
    
    // test save of 1D array
    let a = array![1./8.,2.55555,3./7.,4.,5.];

    save_gnuplot1D(&a, "test_1D_array.txt");

    // test save of 2D array
    let b = arr2(&[[1,2,3,4,5],
                    [2,3,4,5,6],
                    [3,4,5,6,7],
                    [4,5,6,7,8],
                    [5,6,7,8,9]]);
    save_gnuplot2D(&b, "test_2D_array.txt");
    
}