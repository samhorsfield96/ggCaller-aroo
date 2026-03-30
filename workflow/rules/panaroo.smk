rule panaroo:
    """Run Panaroo for pangenome gene clustering on ggCaller GFF output."""
    input:
        gff_dir=config["output_dir"] + "/ggcaller",
    output:
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    params:
        cli_args=config["panaroo_cli_args"],
        out_dir=config["output_dir"] + "/panaroo",
    threads: 40
    conda:
        "envs/panaroo.yaml"
    log:
        config["output_dir"] + "/logs/panaroo.log",
    shell:
        "panaroo -i {input.gff_dir}/GFF/* --threads {threads} -o {params.out_dir} {params.cli_args} > {log} 2>&1"
