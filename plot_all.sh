PLOT_DIR=./plots
BIN=$(realpath ./target/release/ising)

cargo build --release &&

# d = 2

(
    mkdir -p $PLOT_DIR/2d/exact
    cd $PLOT_DIR/2d/exact

    for beta in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    do
        $BIN -n 4 --beta $beta exact
        $BIN -n 4 --beta $beta exact --symmetry-breaking true  # FIXME: Check clap
    done
) &&

(
    mkdir -p $PLOT_DIR/2d/mcmc
    cd $PLOT_DIR/2d/mcmc
    for beta in 0.1 0.42 0.43 0.44 0.45 1.0
    do
        $BIN -n 100 --beta $beta mcmc --iter 50000
    done
) &&

# d = 3

(
    mkdir -p $PLOT_DIR/3d
    cd $PLOT_DIR/3d
    for beta in 0.1 0.2 0.21 0.22 0.23 0.24 1.0
    do
        $BIN -d 3 -n 10 --beta $beta mcmc --iter 50000
    done
)
