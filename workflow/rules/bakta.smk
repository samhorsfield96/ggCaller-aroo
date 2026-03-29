rule translate_representatives:
    """Translate representative nucleotide sequences from Panaroo to amino acids."""
    input:
        config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/translated/pan_genome_reference.faa",
    log:
        config["output_dir"] + "/logs/translate_representatives.log",
    conda:
        "envs/bakta.yaml"
    script:
        "../scripts/translate_nucleotide.py"

rule bakta_proteins:
    """Annotate representative proteins from Panaroo using bakta_proteins."""
    input:
        proteins=config["output_dir"] + "/translated/pan_genome_reference.faa",
    output:
        config["output_dir"] + "/bakta_proteins/proteins.tsv",
    threads: 40
    params:
        db=config["bakta"]["db"],
        out_dir=config["output_dir"] + "/bakta_proteins",
        prefix="proteins",
    log:
        config["output_dir"] + "/logs/bakta_proteins.log",
    conda:
        "envs/bakta.yaml"
    shell:
        "bakta_proteins --force --db {params.db} --output {params.out_dir} --threads {threads} "
        "--prefix {params.prefix} {input.proteins} > {log} 2>&1"
