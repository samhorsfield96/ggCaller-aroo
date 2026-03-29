configfile: "config/config.yaml"

include: "workflow/rules/ggcaller.smk"
include: "workflow/rules/panaroo.smk"
include: "workflow/rules/bakta.smk"


rule all:
    input:
        config["output_dir"] + "/bakta_proteins/proteins.tsv",
