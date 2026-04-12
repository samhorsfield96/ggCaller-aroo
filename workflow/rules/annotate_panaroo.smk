rule annotate_pan_ref:
    """Propagate bakta_proteins annotations to panaroo pangenome references output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/pan_genome_reference.fa"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="pan_ref",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_pan_ref.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_dna_CDS:
    """Propagate bakta_proteins annotations to panaroo CDS output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/combined_DNA_CDS.fasta"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="dna_cds",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_dna_CDS.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_dna_prot:
    """Propagate bakta_proteins annotations to panaroo protein output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/combined_protein_CDS.fasta"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="prot_cds",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_prot_CDS.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_gene_data:
    """Propagate bakta_proteins annotations to panaroo gene_data output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/gene_data.csv"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="gene_data",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_gene_data.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_gpa:
    """Propagate bakta_proteins annotations to panaroo gene presence/absence output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/gene_presence_absence.csv"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="gpa",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_gpa.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_gpa_roary:
    """Propagate bakta_proteins annotations to panaroo roary gene presence/absence output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/gene_presence_absence_roary.csv"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="gpa_roary",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_gpa_roary.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"

rule annotate_gml:
    """Propagate bakta_proteins annotations to panaroo final gml output files."""
    input:
        bakta_tsv=config["output_dir"] + "/bakta_proteins/proteins.tsv",
        pan_ref=config["output_dir"] + "/panaroo/pan_genome_reference.fa",
    output:
        config["output_dir"] + "/annotated/final_graph.gml"
    params:
        panaroo_dir=config["output_dir"] + "/panaroo",
        output_type="gml",
    log:
        config["output_dir"] + "/logs/annotate_panaroo_gml.log",
    threads: 1
    conda:
        "envs/bakta.yaml"
    shell:
        "python ../scripts/annotate_panaroo.py --bakta_tsv {input.bakta_tsv} --output {output} --logfile {log} --output_type {params.output_type} --panaroo_dir {params.panaroo_dir}"