use crate::phasevector::PhaseVector;
use num_complex::{Complex, ComplexFloat};
use rand::distributions::Uniform;
use rand::prelude::*;
use rand::rngs::ThreadRng;
use std::f64::consts::PI;

const UNIT_VECTORS: [[usize; 4]; 4] = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
const ACCEPTANCE_CONSTANT: f64 = 0.2105137;

#[derive(Clone, Debug)]
pub struct Lattice {
    /* the actual lattice holding the configuration */
    lattice: Vec<Vec<Vec<Vec<PhaseVector>>>>,
    width: usize,
    index_distribution: Uniform<usize>,
    link_distribution: Uniform<usize>,
}

impl Lattice {
    pub fn new_uniform(width: usize) -> Self {
        Self {
            lattice: vec![vec![vec![vec![PhaseVector::new_uniform(); width]; width]; width]; width],
            width: width,
            index_distribution: Uniform::from(0..width),
            link_distribution: Uniform::from(0..4),
        }
    }

    pub fn new_random(width: usize, rng: &mut ThreadRng) -> Self {
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

    pub fn heatbath_update(&mut self, beta: f64, rng: &mut ThreadRng) {
        let random_i = self.index_distribution.sample(rng);
        let random_j = self.index_distribution.sample(rng);
        let random_k = self.index_distribution.sample(rng);
        let random_l = self.index_distribution.sample(rng);
        let random_m = self.link_distribution.sample(rng);

        let other_plaquettes =
            self.plaquettes_without_link(random_i, random_j, random_k, random_l, random_m);
        let alpha = other_plaquettes.abs();
        let theta_0 = -other_plaquettes.arg();

        let new_theta = sample_theta(alpha, beta, rng);

        (*self).lattice[random_i][random_j][random_k][random_l].phases[random_m] =
            new_theta + theta_0;
    }
}

fn acceptance_probability(x: f64, prefactor: f64) -> f64 {
    return ((((1.0 - x) * PI / 2.0).cos() - ACCEPTANCE_CONSTANT) * prefactor - x).exp();
}

pub fn sample_theta(alpha: f64, beta: f64, rng: &mut ThreadRng) -> f64 {
    let prefactor = alpha * beta;

    loop {
        let sample_x = -1.0
            + (1.0 / prefactor) * (1.0 + ((2.0 * prefactor).exp() - 1.0) * rng.gen::<f64>()).ln();

        if rng.gen::<f64>() < acceptance_probability(sample_x, prefactor) {
            let mut theta = (PI / 2.0) * (1.0 - sample_x);
            if rng.gen_bool(0.5) {
                theta = -theta;
            }

            return theta;
        }
    }
}
