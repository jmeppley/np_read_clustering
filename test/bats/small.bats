setup() {
    eval "$(conda shell.bash hook)"
    ENV_DIR=env
    ENV_FILE=conda/snake.yaml
    mkdir -p test/outputs
    if [ "$ENV_FILE" -nt "$ENV_DIR/conda-meta/history" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p ./$ENV_DIR --force --quiet > test/outputs/.conda.$ENV_DIR.log 2>&1
    fi
    conda activate ./$ENV_DIR
}

@test "Testing end to end with a small input" {
    # clear any old output
    rm -rf test/outputs/nprc.small
 
    # run whole thing
    run bash -c "snakemake --configfile config.yaml --use-conda --conda-frontend mamba -p -j 2 > test/outputs/nprc.small.log 2>&1"
    [ "$status" -eq 0 ]
 
    # run again (it should do nothing)
    run bash -c "snakemake --configfile config.yaml --use-conda -p -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}
