rule ggcaller:
    """Run ggCaller for graph-based gene prediction on input genome assemblies."""
    input:
        samples=(
            config["refs"] if config.get("refs") is not None
            else config["reads"]
        ),
    output:
        directory(config["output_dir"] + "/ggcaller"),
    params:
        cli_args=config["ggcaller_cli_args"],
        input_flag=(
            "--refs" if config.get("refs") is not None else "--reads"
        ),
    conda:
        "envs/ggcaller.yaml"
    threads: 40
    log:
        config["output_dir"] + "/logs/ggcaller.log",
    shell:
        "ggcaller {params.input_flag} {input.samples} {params.cli_args} "
        "--out {output} --threads {threads} --gene-finding-only > {log} 2>&1"