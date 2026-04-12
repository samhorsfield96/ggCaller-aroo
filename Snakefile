configfile: "config/config.yaml"

include: "workflow/rules/ggcaller.smk"
include: "workflow/rules/panaroo.smk"
include: "workflow/rules/bakta.smk"
include: "workflow/rules/annotate_panaroo.smk"

rule all:
    input:
        config["output_dir"] + "/annotated/pan_genome_reference.fa",
        config["output_dir"] + "/annotated/combined_DNA_CDS.fasta",
        config["output_dir"] + "/annotated/combined_protein_CDS.fasta",
        config["output_dir"] + "/annotated/gene_data.csv",
        config["output_dir"] + "/annotated/gene_presence_absence.csv",
        config["output_dir"] + "/annotated/gene_presence_absence_roary.csv",
        config["output_dir"] + "/annotated/final_graph.gml",
