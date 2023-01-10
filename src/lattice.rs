use crate::phasevector::PhaseVector;
use anyhow::Ok;
use num_complex::{Complex, ComplexFloat};
use fastrand::Rng;
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;

const UNIT_VECTORS: [[usize; 4]; 4] = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
const ACCEPTANCE_CONSTANT: f64 = 0.2105137;

#[derive(Clone, Debug)]
pub struct Lattice {
    /* the actual lattice holding the configuration */
    lattice: Vec<Vec<Vec<Vec<PhaseVector>>>>,
    width: usize,
}

impl Lattice {
    pub fn new_uniform(width: usize) -> Self {
        Self {
            lattice: vec![vec![vec![vec![PhaseVector::new_uniform(); width]; width]; width]; width],
            width: width,
        }
    }

    pub fn new_random(width: usize, rng: &mut Rng) -> Self {
        let mut new_lattice = Lattice::new_uniform(width);

        for axis_1 in new_lattice.lattice.iter_mut() {
            for axis_2 in axis_1.iter_mut() {
                for axis_3 in axis_2.iter_mut() {
                    for phase_vector in axis_3.iter_mut() {
                        *phase_vector = PhaseVector::new_random(rng);
                    }
                }
            }
        }

        return new_lattice;
    }

    /* compute the average action per plaquette */
    pub fn average_action(&self) -> f64 {
        let mut sum = 0f64;
        /* in 4d there are 6 plaquettes per vertex */
        let num_plaquettes = (6 * self.width.pow(4)) as f64;

        /* Sum over all vertices and plaquettes at those vertices */
        for i in 0..self.width {
            for j in 0..self.width {
                for k in 0..self.width {
                    for l in 0..self.width {
                        for m in 0..3 {
                            for n in m + 1..4 {
                                let phase1 = self.lattice[i][j][k][l].phases[m]; /* U_\mu(n) */
                                let phase2 = self.lattice[(i + UNIT_VECTORS[m][0]) % self.width]
                                    [(j + UNIT_VECTORS[m][1]) % self.width]
                                    [(k + UNIT_VECTORS[m][2]) % self.width]
                                    [(l + UNIT_VECTORS[m][3]) % self.width]
                                    .phases[n]; /* U_\nu(n+ \hat{\mu}) */
                                let phase3 = self.lattice[(i + UNIT_VECTORS[n][0]) % self.width]
                                    [(j + UNIT_VECTORS[n][1]) % self.width]
                                    [(k + UNIT_VECTORS[n][2]) % self.width]
                                    [(l + UNIT_VECTORS[n][3]) % self.width]
                                    .phases[m]; /* U_\mu(n+ \hat{\nu}) */
                                let phase4 = self.lattice[i][j][k][l].phases[n];
                                /* U_\nu(n) */

                                sum += 1.0 - (phase1 + phase2 - phase3 - phase4).cos();
                                /* take complex conjugate of last two */
                            }
                        }
                    }
                }
            }
        }

        return sum / num_plaquettes;
    }

    fn plaquettes_without_link(
        &self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        m: usize,
    ) -> Complex<f64> {
        let mut lambda_sum = Complex::from_polar(0.0, 0.0);

        for n in 0..4 {
            if m != n {
                let phase1 = self.lattice[(i + UNIT_VECTORS[m][0]) % self.width]
                    [(j + UNIT_VECTORS[m][1]) % self.width][(k + UNIT_VECTORS[m][2]) % self.width]
                    [(l + UNIT_VECTORS[m][3]) % self.width]
                    .phases[n]; /* U_\nu(n+ \hat{\mu}) */
                let phase2 = self.lattice[(i + UNIT_VECTORS[n][0]) % self.width]
                    [(j + UNIT_VECTORS[n][1]) % self.width][(k + UNIT_VECTORS[n][2]) % self.width]
                    [(l + UNIT_VECTORS[n][3]) % self.width]
                    .phases[m]; /* U_\mu(n+ \hat{\nu}) */
                let phase3 = self.lattice[i][j][k][l].phases[n]; /* U_\mu(n) */

                let lambda1 = Complex::from_polar(1.0, phase1 - phase2 - phase3);
                lambda_sum += lambda1;

                let phase4 = self.lattice[(i + self.width - UNIT_VECTORS[n][0]) % self.width]
                    [(j + self.width - UNIT_VECTORS[n][1]) % self.width]
                    [(k + self.width - UNIT_VECTORS[n][2]) % self.width]
                    [(l + self.width - UNIT_VECTORS[n][3]) % self.width]
                    .phases[m]; /* U_\mu(n+ \hat{\nu}) */
                let phase5 = self.lattice
                    [(i + self.width - UNIT_VECTORS[n][0] + UNIT_VECTORS[m][0]) % self.width]
                    [(j + self.width - UNIT_VECTORS[n][1] + UNIT_VECTORS[m][1]) % self.width]
                    [(k + self.width - UNIT_VECTORS[n][2] + UNIT_VECTORS[m][2]) % self.width]
                    [(l + self.width - UNIT_VECTORS[n][3] + UNIT_VECTORS[m][3]) % self.width]
                    .phases[n]; /* U_\nu(n - \hat{\nu} + \hat{\mu}) */
                let phase6 = self.lattice[(i + self.width - UNIT_VECTORS[n][0]) % self.width]
                    [(j + self.width - UNIT_VECTORS[n][1]) % self.width]
                    [(k + self.width - UNIT_VECTORS[n][2]) % self.width]
                    [(l + self.width - UNIT_VECTORS[n][3]) % self.width]
                    .phases[n]; /* U_\nu(n - \hat{\nu}) */

                let lambda2 = Complex::from_polar(1.0, -phase4 - phase5 + phase6);
                lambda_sum += lambda2;
            }
        }
        return lambda_sum;
    }

    pub fn heatbath_sweep(&mut self, beta: f64, rng: &mut Rng) {
        for i in 0..self.width {
            for j in 0..self.width {
                for k in 0..self.width {
                    for l in 0..(self.width) {
                        for m in 0..4 {
                            let other_plaquettes = self.plaquettes_without_link(i, j, k, l, m);
                            let alpha = other_plaquettes.abs();
                            let theta_0 = -other_plaquettes.arg();

                            let new_theta = sample_theta(alpha, beta, rng);

                            (*self).lattice[i][j][k][l].phases[m] = new_theta + theta_0;
                        }
                    }
                }
            }
        }
    }

    pub fn visualize_3d_lattice(&self, file: &mut File) -> anyhow::Result<()>  {
        writeln!(file, "\\tdplotsetmaincoords{{22}}{{22}}")?;
        writeln!(file, "\\begin{{tikzpicture}}[tdplot_main_coords]")?;
        let plane_index = self.width / 2; // take a plane somewhere in the middle

        for i in 0..self.width {
            for j in 0..self.width {
                for k in 0..self.width{
                writeln!(file, "\\filldraw[black] ({},{},{}) circle (2pt) ;", i,j,k)?;
                let color_x1 = phase_to_rgb(self.lattice[i][j][k][plane_index].phases[0]);
                let color_x2 = phase_to_rgb(self.lattice[i][j][k][plane_index].phases[1]);
                let color_x3 = phase_to_rgb(self.lattice[i][j][k][plane_index].phases[2]);

                writeln!(file, "\\definecolor{{color{}{}{}1}}{{RGB}}{{{},{},{}}} ;",i,j,k, color_x1.0, color_x1.1, color_x1.2)?;
                writeln!(file, "\\draw[color{0}{1}{2}1, thick] ({0},{1},{2}) -- ({3},{1},{2}) ;",i,j,k,i+1)?;
                writeln!(file, "\\definecolor{{color{}{}{}2}}{{RGB}}{{{},{},{}}} ;",i,j,k, color_x2.0, color_x2.1, color_x2.2)?;
                writeln!(file, "\\draw[color{0}{1}{2}2, thick] ({0},{1},{2}) -- ({0},{3},{2}) ;",i,j,k,j+1)?;

                writeln!(file, "\\definecolor{{color{}{}{}3}}{{RGB}}{{{},{},{}}} ;",i,j,k, color_x3.0, color_x3.1, color_x3.2)?;
                writeln!(file, "\\draw[color{0}{1}{2}3, thick] ({0},{1},{2}) -- ({0},{1},{3}) ;",i,j,k,k+1)?;
            }
        }
    }
        writeln!(file, "\\end{{tikzpicture}}")?;
        Ok(())
    }

    pub fn visualize_plaquettes_plane(&self, file: &mut File) -> anyhow::Result<()> {
        writeln!(file, "\\begin{{tikzpicture}}")?;
        let plane_index = self.width/2;

        for i in 0..self.width {
            for j in 0..self.width {

                let mut plaquette = self.lattice[i][j][plane_index][plane_index].phases[0]
                                        +self.lattice[(i + 1) % self.width][j][plane_index][plane_index].phases[1]
                                        -self.lattice[i][(j + 1) % self.width][plane_index][plane_index].phases[0]
                                        -self.lattice[i][j][plane_index][plane_index].phases[1];

                while plaquette < 0.0 {
                    plaquette += 2.0*PI;
                }

                while plaquette > 2.0 *PI {
                    plaquette -= 2.0 * PI;
                }

                let color = phase_to_rgb(plaquette);
                writeln!(file, "\\definecolor{{color{}{}}}{{RGB}}{{{},{},{}}} ;",i,j, color.0, color.1, color.2)?;
                writeln!(file, "\\fill[color{0}{1}] ({0},{1}) rectangle ({2},{3}) ;",i,j,i+1,j+1)?;
            }
        }

        for i in 0..self.width {
            for j in 0..self.width {
                writeln!(file, "\\filldraw[black] ({},{}) circle (2pt) ;", i ,j)?;
            }
        }
        writeln!(file,"\\end{{tikzpicture}}")?;
        Ok(())
    }

    pub fn visualize_plaquettes_plane_svg(&self, file: &mut File) -> anyhow::Result<()> {
        writeln!(file, "<svg width=\"{0}\" height=\"{0}\">", 50*self.width+20)?;
        let plane_index = self.width/2;

        for i in 0..self.width {
            for j in 0..self.width {

                let mut plaquette = self.lattice[i][j][plane_index][plane_index].phases[0]
                +self.lattice[(i + 1) % self.width][j][plane_index][plane_index].phases[1]
                -self.lattice[i][(j + 1) % self.width][plane_index][plane_index].phases[0]
                -self.lattice[i][j][plane_index][plane_index].phases[1];

                while plaquette < 0.0 {
                    plaquette += 2.0*PI;
                }

                while plaquette > 2.0 *PI {
                    plaquette -= 2.0 * PI;
                }

                let (r,g,b) = phase_to_rgb(plaquette);
                writeln!(file, "<rect x=\"{0}\" y=\"{1}\" width=\"50\" height=\"50\" fill=\"#{r:02X?}{g:02X?}{b:02X?}\"/>",i*50+10, j*50+10)?;
            }
        }

        for i in 0..self.width {
            for j in 0..self.width {
                writeln!(file, "<circle cx=\"{}\" cy=\"{}\" r=\"5\" fill=\"#FFFFF\"/>", i*50+10 ,j*50+10)?;
            }
        }
        writeln!(file,"</svg>")?;
        Ok(())
    }
}

fn acceptance_probability(x: f64, prefactor: f64) -> f64 {
    return ((((PI/2.0)*(1.0-x)).cos() - x) * prefactor).exp() / (ACCEPTANCE_CONSTANT * prefactor).exp();
}

pub fn sample_theta(alpha: f64, beta: f64, rng: &mut Rng) -> f64 {
    let prefactor = alpha * beta;

    loop {
        let sample_x = -1.0
            + (1.0 / prefactor) * (1.0 + ((2.0 * prefactor).exp() - 1.0) * rng.f64()).ln();

        if rng.f64() < acceptance_probability(sample_x, prefactor) {
            let mut theta = (PI / 2.0) * (1.0 - sample_x);
            if rng.bool() {
                theta = -theta;
            }

            return theta;
        }
    }
}

fn phase_to_rgb(phi: f64) -> (u8,u8,u8)  {
    let division = PI / 3.0;
     if phi >= 0.0 && phi <= division {
        return (255, (phi * 255.0 / division) as u8, 0);
    } else if phi > division && phi <= 2.0 * division {
        return ( 255 - ((phi - division) * 255.0 / division) as u8 ,255,0);
    } else if phi > 2.0 * division && phi <= 3.0 * division {
        return (0,255, ((phi - 2.0 * division) * 255.0 / division) as u8);    
    } else if phi > 3.0 * division && phi <= 4.0 * division {
        return (0,255 - ((phi - 3.0 * division) * 255.0 / division) as u8 ,255);
    } else if phi > 4.0 * division && phi <= 5.0 * division {
        return (((phi - 4.0 * division) * 255.0 / division) as u8,0,255);
    } else if phi > 5.0 * division && phi <= 6.0 * division {
        return (255,0,255 - ((phi - 5.0 * division)* 255.0 /division) as u8);
    }
    
    else {
        return (0,0,0);
    }
}