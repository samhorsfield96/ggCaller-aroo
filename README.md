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
      │
      ▼
Annotate Panaroo ── update Panaroo outputs with functional annotation
```

1. **ggCaller** – graph-based gene prediction on bacterial genome assemblies.
2. **Panaroo** – pangenome gene clustering; produces `pan_genome_reference.fa`
   (one representative nucleotide sequence per gene cluster).
3. **Translation** – converts representative nucleotide sequences to amino acids.
4. **Bakta** – annotates representative proteins with `bakta_proteins`.
5. **Annotate Panaroo** - uses Bakta representatives to update annotations in Panaroo outputs

## From docker

### Installation 

Install [docker](https://docs.docker.com/get-started/get-docker/), then run:

```
docker pull samhorsfield96/ggcallaroo:main
```

### Usage

To run within the container, navigate to the directory containing the FASTA files (e.g. `work_dir`) and run the below code to generate the input files:

```
cd work_dir && ls -d -1 $PWD/*.fa > refs.txt && sed "s|^.*/|/data/|" refs.txt > refs_docker.txt
```

Then run inside the docker container, changing other parameters as required:

```
docker run -v $(pwd):/output -v $(pwd):/data samhorsfield96/ggcallaroo:main snakemake --cores 4 --use-conda --conda-frontend mamba --config refs=/data/refs_docker.txt ggcaller_cli_args="--save" panaroo_cli_args="--clean-mode moderate" bakta_db=bakta_db/db-light output_dir=/output/results
```

The results will be generated in `results`

## From source

### Installation

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/ggCallaroo.git
cd ggCallaroo
mamba env create -n ggcallaroo "snakemake>=9.19.0" "bakta>=1.12.0"
mamba activate ggcallaroo
```

### Usage

#### 1. Download bakta database

Install Bakta and download the latest bakta database:

`bakta_db download --output <output-path> --type [light|full]`

Add the bakta_db filepath to the config.yaml file.

#### 2. Prepare a samples file

Create a plain-text file (e.g. `input.txt`) with one genome assembly FASTA
path per line:

```
/path/to/sample1.fasta
/path/to/sample2.fasta
/path/to/sample3.fasta
```

You can do this using the command:

`ls -d -1 $PWD/*.fasta > input.txt`

#### 3. Edit `config/config.yaml`

Set `samples` to your samples file, `output_dir` to the desired results
directory, and supply any extra tool arguments via the `cli_args` fields.
Set `bakta.db` to your local Bakta database path.

#### 4. Run the pipeline

```bash
snakemake --cores <N> --use-conda
```

If using mamba, run instead:

```bash
snakemake --cores <N> --use-conda --conda-frontend mamba
```

## Configuration

All tool-specific parameters are controlled via `config/config.yaml`.

| Key | Description |
|-----|-------------|
| `refs` or `reads` | Path to reference genome files or reads FASTA files (one per line) |
| `output_dir` | Base output directory |
| `ggcaller_cli_args` | CLI arguments passed verbatim to `ggcaller` |
| `panaroo_cli_args` | CLI arguments passed verbatim to `panaroo` |
| `bakta_db` | Path to Bakta database directory |

## Output

| Path | Description |
|------|-------------|
| `{output_dir}/ggcaller/` | ggCaller output (GFF files, protein FASTAs, etc.) |
| `{output_dir}/panaroo/` | Panaroo output (pangenome files) |
| `{output_dir}/bakta_proteins/` | Bakta annotation output |
| `{output_dir}/annotated/` | Directory of annotated Panaroo outputs |
| `{output_dir}/logs/` | Log files for each pipeline step |

