use std::cmp;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{Write, Error, BufReader};
use rand::Rng;
use std::io::{self, BufRead};
use std::fs::read_to_string;
use std::path::Path;
use ::std::cmp::Reverse;
use ordered_float::NotNan;

const N_i: usize = 2048;

fn generate_file(num: i64, name: &str) -> Result<(), Error> {
    let mut rng = rand::thread_rng();
    let mut output = File::create(name.to_string())?;
    for _i in 1..num {
        let f: f64 = rng.gen_range(-100_000.0..100_000.0);
        write!(output, "{f}\n");
    }

    Ok(())
}

struct Stats {
    filename: String,
    sample_mean: f64,
    sample_dispersion: f64,
    min: f64,
    max: f64,
    median: f64,
    left_tail: BinaryHeap<NotNan<f64>>,
    right_tail: BinaryHeap<Reverse<NotNan<f64>>>,
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}


impl Stats {
    fn get_len(filename: &str) -> Result<usize, Error> {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines().count())
    }


    fn stats(filename: &str) -> Self {
        let len = Stats::get_len(filename).unwrap();
        let heap_size = (len as f64 * 0.1) as usize;
        let mut amount_of_buckets = len / N_i; //each bucket stores average of N_i elements
        if len % N_i != 0 {
            amount_of_buckets += 1;
        }
        let mut average_buckets: Vec<f64> = vec![0.0];
        let mut median_buckets: Vec<f64> = vec![0.0];
        let mut dispersion_buckets: Vec<f64> = vec![0.0];
        let mut buffer: Vec<f64> = vec![];
        let mut lines_counter = 0;
        let mut buck_ind: usize = 0;
        let mut min = f64::MAX;
        let mut max = f64::MIN;
        let mut left_tail: BinaryHeap<NotNan<f64>> = BinaryHeap::new();
        let mut right_tail: BinaryHeap<Reverse<NotNan<f64>>> = BinaryHeap::new();

        if let Ok(lines) = read_lines(filename) {
            for line in lines {
                if let Ok(ip) = line {
                    let cur = ip.parse::<f64>().unwrap();
                    min = f64::min(min, cur);
                    max = f64::max(max, cur);
                    let not_nan_cur = NotNan::new(cur);
                    if left_tail.is_empty() || &not_nan_cur.unwrap() < left_tail.peek().unwrap() {
                        if left_tail.len() == heap_size {
                            left_tail.pop();
                        }
                        left_tail.push(not_nan_cur.unwrap());
                    }

                    if right_tail.is_empty() || &Reverse(not_nan_cur.unwrap()) < right_tail.peek().unwrap() {
                        if right_tail.len() == heap_size {
                            right_tail.pop();
                        }
                        right_tail.push(Reverse(NotNan::new(cur).unwrap()));
                    }

                    buffer.push(cur);
                    lines_counter += 1;
                    average_buckets[buck_ind] += cur;
                    if lines_counter == N_i {
                        average_buckets[buck_ind] /= (N_i as f64);
                        buffer.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        median_buckets[buck_ind] = buffer[N_i / 2];
                        for i in &buffer {
                            dispersion_buckets[buck_ind] += (i - average_buckets[buck_ind]).powi(2);
                        }
                        dispersion_buckets[buck_ind] /= (N_i as f64);
                        lines_counter = 0;
                        buck_ind += 1;
                        average_buckets.push(0.0);
                        dispersion_buckets.push(0.0);
                        median_buckets.push(0.0);
                        buffer.clear();
                    }
                }
            }
        }

        if lines_counter != 0 {
            average_buckets[buck_ind] /= lines_counter as f64;
            buffer.sort_by(|a, b| a.partial_cmp(b).unwrap());
            median_buckets[buck_ind] = buffer[lines_counter / 2];

            for i in buffer {
                dispersion_buckets[buck_ind] += (i - average_buckets[buck_ind]).powi(2);
            }
            dispersion_buckets[buck_ind] /= (lines_counter as f64);
        }

        let mut sample_mean: f64 = 0.0;
        let mut sample_dispersion: f64 = 0.0;
        let mut median: f64 = 0.0;
        if amount_of_buckets == 1 {
            sample_mean = average_buckets[0];
            sample_dispersion = dispersion_buckets[0];
            median = median_buckets[0];
        } else {
            for ind in 0..buck_ind {
                sample_mean += average_buckets[ind] / ((buck_ind) as f64);
            }
            for ind in 0..buck_ind {
                sample_dispersion += (dispersion_buckets[ind] + (average_buckets[ind] - sample_mean).powi(2)) / (buck_ind as f64);
            }
            let mut prev_mean = sample_mean;
            sample_mean = (sample_mean * ((len - lines_counter) as f64) + (average_buckets[buck_ind] * (lines_counter as f64))) / (len as f64);

            sample_dispersion = (((len - lines_counter) as f64) * (sample_dispersion + (prev_mean - sample_mean).powi(2)) + (lines_counter as f64) * (dispersion_buckets[buck_ind] + (sample_mean - average_buckets[buck_ind]).powi(2))) / (len as f64);

            median_buckets.sort_by(|a, b| a.partial_cmp(b).unwrap());
            median = median_buckets[median_buckets.len() / 2];
        }

        Self {
            filename: filename.to_string(),
            sample_mean,
            sample_dispersion,
            min,
            max,
            median,
            left_tail,
            right_tail,
        }
    }
}

//TODO: пофиксить медиану
fn main() -> io::Result<()> {
    let filename = "small.txt";
    //generate_file(133_875, filename);
    let st = Stats::stats(filename);
    let len = Stats::get_len(filename).unwrap() as f64;
    let mut av2: f64 = 0.0;
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(ip) = line {
                av2 += ip.parse::<f64>().unwrap();
            }
        }
    }
    av2 /= len;
    println!("sample mean: {}, {}", st.sample_mean, av2);
    let mut disp2: f64 = 0.0;
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(ip) = line {
                disp2 += (ip.parse::<f64>().unwrap() - av2).powi(2);
            }
        }
    }
    disp2 /= len;
    println!("sample dispersion: {}, {}", st.sample_dispersion, disp2);

    // println!("left_tail: ");
    // for i in st.left_tail.into_sorted_vec() {
    //     println!("{i}");
    // }
    // println!();
    // println!("right_tail: ");
    // for i in st.right_tail.into_sorted_vec() {
    //     println!("{}", i.0);
    // }
    // println!();
    Ok(())
}