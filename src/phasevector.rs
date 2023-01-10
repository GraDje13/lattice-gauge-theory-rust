use fastrand::Rng;
use std::f64::consts::PI;

#[derive(Copy, Clone, Debug)]
pub struct PhaseVector {
    pub phases: [f64; 4], /* simple struct to hold phases */
}

impl PhaseVector {
    pub fn new_uniform() -> Self {
        Self { phases: [0.0; 4] }
    }

    pub fn new_random(rng: &mut Rng) -> Self {
        let mut new_phase_vector = Self::new_uniform();

        for phase in new_phase_vector.phases.iter_mut() {
            *phase = rng.f64() * 2.0 * PI;
        }

        return new_phase_vector;
    }
}
