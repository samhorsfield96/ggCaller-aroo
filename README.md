# ggCallaroo

A Snakemake pipeline combining graph-based gene prediction with
[ggCaller](https://github.com/bacpop/ggCaller), pangenome gene clustering
with [Panaroo](https://github.com/gtonkinhill/panaroo), and functional annotation
of representative proteins with [Bakta](https://github.com/oschwengers/bakta).

## Pipeline overview

```
Input assemblies
      │
      ▼
 ggCaller  ── gene prediction (GFF + protein FASTA per sample)
      │
      ▼
  Panaroo   ── pangenome clustering → pan_genome_reference.fa
      │
      ▼
 Translation ── nucleotide → amino-acid FASTA (Biopython)
      │
      ▼
bakta_proteins ── functional annotation of representative proteins
```

1. **ggCaller** – graph-based gene prediction on bacterial genome assemblies.
2. **Panaroo** – pangenome gene clustering; produces `pan_genome_reference.fa`
   (one representative nucleotide sequence per gene cluster).
3. **Translation** – converts representative nucleotide sequences to amino acids.
4. **Bakta** – annotates representative proteins with `bakta_proteins`.

## Installation

### via docker

Install [docker](https://docs.docker.com/get-started/get-docker/), then run:

```
docker pull samhorsfield96/ggcallaroo:main
```

To run within the container, use the below command, replacing `path to output dir` and `path to workdir` with absolute paths and changing other parameters as required:

```
docker run -v <path to output dir>:/output -v <path to workdir>:/data samhorsfield96/ggcallaroo:main snakemake --cores 4 --config refs=/data/refs.txt ggcaller_cli_args="--save" panaroo_cli_args="--clean-mode moderate" bakta_db=bakta_db/db-light output_dir=/output/results
```

### From source

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/ggCallaroo.git
cd ggCallaroo
mamba env create -n ggcallaroo "snakemake>=9.19.0" "bakta>=1.12.0"
mamba activate ggcallaroo
```

## Usage

### 1. Download bakta database

Install Bakta and download the latest bakta database:

`bakta_db download --output <output-path> --type [light|full]`

Add the bakta_db filepath to the config.yaml file.

### 2. Prepare a samples file

Create a plain-text file (e.g. `input.txt`) with one genome assembly FASTA
path per line:

```
/path/to/sample1.fasta
/path/to/sample2.fasta
/path/to/sample3.fasta
```

You can do this using the command:

`ls -d -1 $PWD/*.fasta > input.txt`

### 3. Edit `config/config.yaml`

Set `samples` to your samples file, `output_dir` to the desired results
directory, and supply any extra tool arguments via the `cli_args` fields.
Set `bakta.db` to your local Bakta database path.

### 4. Run the pipeline

```bash
snakemake --cores <N> --use-conda
```

To perform a dry-run first:

```bash
snakemake --cores <N> -n --use-conda
```

## Configuration

All tool-specific parameters are controlled via `config/config.yaml`.

| Key | Description |
|-----|-------------|
| `input.txt` | Path to text file listing genome assembly FASTA paths (one per line) |
| `output_dir` | Base output directory |
| `ggcaller_cli_args` | CLI arguments passed verbatim to `ggcaller` |
| `panaroo_cli_args` | CLI arguments passed verbatim to `panaroo` |
| `bakta_db` | Path to Bakta database directory |

### Example ggCaller arguments

```yaml
ggcaller:
  cli_args: "--save --kmer 31"
```

### Example Panaroo arguments

```yaml
panaroo:
  cli_args: ""--clean-mode moderate -a core --remove-invalid-genes"
```

## Output

| Path | Description |
|------|-------------|
| `{output_dir}/ggcaller/` | ggCaller output (GFF files, protein FASTAs, etc.) |
| `{output_dir}/panaroo/` | Panaroo output (pangenome files) |
| `{output_dir}/panaroo/pan_genome_reference.fa` | Representative gene sequences (nucleotide) |
| `{output_dir}/translated/pan_genome_reference.faa` | Translated representative proteins |
| `{output_dir}/bakta_proteins/` | Bakta annotation output |
| `{output_dir}/bakta_proteins/proteins.tsv` | Annotation results (TSV) |
| `{output_dir}/logs/` | Log files for each pipeline step |
