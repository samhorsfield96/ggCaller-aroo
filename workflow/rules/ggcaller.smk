rule ggcaller:
    """Run ggCaller for graph-based gene prediction on input genome assemblies."""
    input:
        samples=config["samples"],
    output:
        directory(config["output_dir"] + "/ggcaller"),
    params:
        cli_args=config["ggcaller"]["cli_args"],
    log:
        config["output_dir"] + "/logs/ggcaller.log",
    shell:
        "ggcaller --reads {input.samples} {params.cli_args} --out {output} > {log} 2>&1"
