pub mod lattice;
pub mod phasevector;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use hdf5::File;
use lattice::Lattice;

#[derive(Parser)]
#[command(author, version, about, long_about=None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// resume from old config
    Resume(Resume),
    /// create new config
    New(New),
}

#[derive(Args)]
struct Resume {}

#[derive(Args)]
struct New {
    /// name for new save file
    #[arg(short, long)]
    name: String,

    /// specify value of beta
    #[arg(short, long)]
    beta: f64,

    /// specify lattice width
    #[arg(short, long)]
    lattice_width: usize,

    /// specify if state should start in ordered config
    #[arg(short, long)]
    ordered: bool,

    /// specify number of measurements
    #[arg(short, long)]
    measurements: usize,

    /// specify number of equilibration sweeps
    #[arg(short, long)]
    equilibration_sweeps: usize,

    /// specify number of sweeps between measurements
    #[arg(short, long)]
    sweeps_between_measurements: usize,

    /// specify number of seconds between saves
    #[arg(short, long)]
    measurements_between_saves: usize,
}

fn main() -> Result<()> {
    // initialize the random number generator
    let mut rng = rand::thread_rng();

    // parse the arguments
    let cli = Cli::parse();

    match cli.command {
        Commands::New(settings) => {
            // print settings to user
            println!("Starting new simulation");
            println!("Data will be saved in: {}", settings.name);
            println!("Beta is set to: {}", settings.beta);
            println!("Lattice width is set to {}", settings.lattice_width);
            println!("Ordered start is set to {}", settings.ordered);
            println!(
                "Simulation will perform {} measurements",
                settings.measurements
            );
            println!(
                "Burn in phase is {} sweeps long",
                settings.equilibration_sweeps
            );
            println!(
                "{} sweeps will be performed in between measurements",
                settings.sweeps_between_measurements
            );
            println!(
                "Simulation will be saved every {} measurements",
                settings.measurements_between_saves
            );

            // create the save file, give error if it exists to prevent accidental overwriting of data
            let file = File::create_excl(&settings.name)
                .with_context(|| format!("Failed to create file {}", settings.name))?;

            // create dataset
            let action_dataset = file
                .new_dataset::<f64>()
                .chunk((1, settings.measurements_between_saves))
                .shape((0.., settings.measurements_between_saves))
                .create("action_measurements")?;

            // write attributes
            let beta_attribute = action_dataset.new_attr::<f64>().shape([1]).create("beta")?;
            beta_attribute
                .write(&[settings.beta])
                .with_context(|| format!("failed to write beta"))?;

            let lattice_width_attribute = action_dataset
                .new_attr::<usize>()
                .shape([1])
                .create("lattice-width")?;
            lattice_width_attribute.write(&[settings.lattice_width])?;

            let ordered_attribute = action_dataset
                .new_attr::<bool>()
                .shape([1])
                .create("ordered")?;
            ordered_attribute.write(&[settings.ordered])?;

            let equilibration_sweeps_atttribute = action_dataset
                .new_attr::<usize>()
                .shape([1])
                .create("equilibration_sweeps")?;
            equilibration_sweeps_atttribute.write(&[settings.equilibration_sweeps])?;

            let sweeps_between_measurements_attribute = action_dataset
                .new_attr::<usize>()
                .shape([1])
                .create("sweeps-between-measurements")?;
            sweeps_between_measurements_attribute.write(&[settings.sweeps_between_measurements])?;

            // initialize lattice
            let mut lattice: Lattice;

            if settings.ordered {
                lattice = Lattice::new_uniform(settings.lattice_width);
            } else {
                lattice = Lattice::new_random(settings.lattice_width, &mut rng);
            }

            // burn in phase
            for _ in 0..settings.equilibration_sweeps {
                lattice.heatbath_sweep(settings.beta, &mut rng);
            }

            let mut measurement_vector = Vec::with_capacity(settings.measurements_between_saves);
            let mut save_counter = 0;

            // note that if the amount measurements is not divisible by the amount of measurement between saves some data is lost
            for i in 0..settings.measurements {
                for _ in 0..settings.sweeps_between_measurements {
                    lattice.heatbath_sweep(settings.beta, &mut rng);
                }
                measurement_vector.push(lattice.average_action());

                if (i + 1) % settings.measurements_between_saves == 0 {
                    action_dataset
                        .resize((save_counter + 1, settings.measurements_between_saves))?;
                    action_dataset.write_slice(&measurement_vector, (save_counter, ..))?;
                    measurement_vector.clear();
                    save_counter += 1;
                }
            }

            println!("simulation complete");
            Ok(())
        }
        Commands::Resume(settings) => {
            println!("started with resume");
            todo!();
        }
    }
}
