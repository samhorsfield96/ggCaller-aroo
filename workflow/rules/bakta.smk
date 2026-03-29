rule translate_representatives:
    """Translate representative nucleotide sequences from Panaroo to amino acids."""
    input:
        config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/translated/pan_genome_reference.faa",
    log:
        config["output_dir"] + "/logs/translate_representatives.log",
    script:
        "../scripts/translate_nucleotide.py"


rule bakta_proteins:
    """Annotate representative proteins from Panaroo using bakta_proteins."""
    input:
        proteins=config["output_dir"] + "/translated/pan_genome_reference.faa",
    output:
        config["output_dir"] + "/bakta_proteins/proteins.tsv",
    params:
        db=config["bakta"]["db"],
        cli_args=config["bakta"]["cli_args"],
        out_dir=config["output_dir"] + "/bakta_proteins",
        prefix="proteins",
    log:
        config["output_dir"] + "/logs/bakta_proteins.log",
    shell:
        "bakta_proteins --db {params.db} --output {params.out_dir} "
        "--prefix {params.prefix} {params.cli_args} {input.proteins} > {log} 2>&1"
