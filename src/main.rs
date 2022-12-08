pub mod lattice;
pub mod phasevector;

use clap::{Args, Parser, Subcommand};
use hdf5::{File, Result};
use lattice::{sample_theta, Lattice};

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

    /// specify number
    #[arg(short, long)]
    measurements: usize,
    #[arg(short, long)]
    equilibration_sweeps: usize,
    #[arg(short, long)]
    measurements_between_sweeps: usize,
    #[arg(short, long)]
    seconds_between_saves: u64,
}

fn main() {
    let mut rng = rand::thread_rng();

    // let cli = Cli::parse();

    /* match cli.command {
        Commands::New(settings) => {
            println!("started with new");
        }
        Commands::Resume(settings) => {
            println!("started with resume");
        }
    }*/

    /*
    let file = File::create("test.h5").unwrap();
    let group = file.create_group("dir").unwrap();

    let builder = group.new_dataset_builder();

    let ds = builder.with_data(&samples).create("test").unwrap();
    */

    let mut data_ordered = [0f64; 30000];
    let mut data_unorderd = [0f64; 30000];

    let mut lattice_ordered = Lattice::new_uniform(6);
    let mut lattice_unordered = Lattice::new_random(6, &mut rng);

    for i in 0..30000 {
        for _ in 0..20 {
            lattice_ordered.heatbath_update(0.1, &mut rng);
            lattice_unordered.heatbath_update(0.1, &mut rng);
        }

        data_ordered[i] = lattice_ordered.average_action();
        data_unorderd[i] = lattice_unordered.average_action();
    }

    let file_ordered = File::create("ordered.h5").unwrap();
    let file_unordered = File::create("unordered.h5").unwrap();
    let group_ordered = file_ordered.create_group("dir").unwrap();
    let group_unordered = file_unordered.create_group("dir").unwrap();

    let builder_ordered = group_ordered.new_dataset_builder();
    let builder_unordered = group_unordered.new_dataset_builder();

    let ds1 = builder_ordered
        .with_data(&data_ordered)
        .create("ordered")
        .unwrap();
    let ds2 = (builder_unordered.with_data(&data_unorderd))
        .create("unordered")
        .unwrap();

    println!("this worked");
}
