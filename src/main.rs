use clap::{Parser, Subcommand};

mod gray_flips;

mod ising;
use ising::{
    gray_sampler, mcmc_cluster_sampler, run_ising_simulation, AvgTracker, ExhaustiveAggregator,
    MCMCAggregator,
};

mod lattice;
use lattice::{BoundaryConditions, Lattice};

mod plots;
use plots::{plot_histogram, plot_lattice};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Inverse thermodynamic temperature for the simulation
    #[arg(long)]
    beta: f64,

    /// The dimension d of the lattice (lattice will be a n^d hypercube)
    #[arg(long, short = 'd', default_value_t=2)]
    lattice_dimension: usize,

    /// The width n of the lattice (lattice will be a n^d hypercube)
    #[arg(long, short = 'n')]
    lattice_width: Option<usize>,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute histogram based on exhaustive iteration of states
    Exact {
        /// Use the symmetry breaking boundary conditions
        #[arg(long)]
        symmetry_breaking: Option<bool>,
    },
    /// Compute probabilities using MCMC methods
    MCMC {
        /// Number of steps of the Markov chain
        #[arg(long, short)]
        iter: Option<usize>,
    },
}

fn main() {
    let cli = Cli::parse();

    let series = match &cli.command {
        Commands::Exact { symmetry_breaking } => {
            let boundary_conditions = if symmetry_breaking.unwrap_or(false) {
                BoundaryConditions::SymmetryBreaking
            } else {
                BoundaryConditions::Periodic
            };
            let lattice = Lattice::new(cli.lattice_dimension, cli.lattice_width.unwrap_or(4), boundary_conditions);
            run_ising_simulation(
                gray_sampler(lattice),
                vec![],
                ExhaustiveAggregator::new(lattice, cli.beta),
            )
        }
        Commands::MCMC { iter } => {
            let iterations = iter.unwrap_or(25_000);
            let lattice_width = cli.lattice_width.unwrap_or(100);
            let lattice = Lattice::new(cli.lattice_dimension, lattice_width, BoundaryConditions::Periodic);
            let mut tracker = AvgTracker::new(iterations / 100);
            run_ising_simulation(
                mcmc_cluster_sampler(iterations, lattice, cli.beta, &|spins| {
                    plot_lattice(
                        format!(
                            "[d={}, β={}, n={}, it={}] Spins lattice",
                            cli.lattice_dimension, cli.beta, lattice_width, iterations
                        ),
                        lattice,
                        spins,
                    )
                }),
                vec![&mut tracker],
                MCMCAggregator::new(iterations, lattice, cli.beta),
            )
        }
    };
    series.into_iter().for_each(plot_histogram);
}
