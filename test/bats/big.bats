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

@test "Testing step_1 with larger input" {
    # clear any old output
    #rm -rf test/outputs/nprc

    # run step 1
    run bash -c "snakemake --configfile config.yaml -p --use-conda --conda-frontend mamba --config all_fasta=test/test.fasta work_dir=test/outputs/nprc -j 20 step_1 > test/outputs/nprc.step1.log 2>&1"
    [ "$status" -eq 0 ]

    # run again (it should do nothing)
    run bash -c "snakemake --configfile config.yaml -p --use-conda --config all_fasta=test/test.fasta work_dir=test/outputs/nprc -j 20 -n step_1 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Testing group 1" {
    run bash -c "snakemake  --configfile config.yaml -p --use-conda \
                            --config group=1 all_fasta=test/test.fasta work_dir=test/outputs/nprc \
                            -j 20 step_4 > test/outputs/nprc.step4.g1.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Testing group 0" {
    run bash -c "snakemake  --configfile config.yaml -p --use-conda \
                            --config group=0 all_fasta=test/test.fasta work_dir=test/outputs/nprc \
                            -j 20 step_4 > test/outputs/nprc.step4.g0.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Testing group 4" {
    run bash -c "snakemake  --configfile config.yaml -p --use-conda \
                            --config group=2 all_fasta=test/test.fasta work_dir=test/outputs/nprc \
                            -j 20 step_4 > test/outputs/nprc.step4.g2.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Testing remaining groups" {
    run bash -c "snakemake  --configfile config.yaml -p --use-conda \
                            --config all_fasta=test/test.fasta work_dir=test/outputs/nprc \
                            -j 20 step_4 > test/outputs/nprc.step4.gs.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Testing finish" {
    # run to finish
    run bash -c "snakemake --configfile config.yaml -p \
                           --use-conda --conda-frontend mamba \
                           --config all_fasta=test/test.fasta work_dir=test/outputs/nprc -j 20 finish > test/outputs/nprc.finish.log 2>&1"
    [ "$status" -eq 0 ]
}
