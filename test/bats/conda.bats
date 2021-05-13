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

@test "Testing conda setup" {
    # clear any old envs
    rm -rf .snakemake/conda conda/.*.ready

    # make sure folder for log file exists
    mkdir -p test/outputs

    run bash -c "snakemake --use-conda --conda-frontend mamba -p -j 2 condaprep > test/outputs/condaprep.log 2>&1"
    [ "$status" -eq 0 ]
}

