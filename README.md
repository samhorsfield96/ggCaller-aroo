# ggCalleraroo

A Snakemake pipeline combining graph-based gene prediction with
[ggCaller](https://github.com/samhorsfield96/ggCaller), pangenome gene clustering
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

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0
- [ggCaller](https://github.com/samhorsfield96/ggCaller)
- [Panaroo](https://github.com/gtonkinhill/panaroo)
- [Bakta](https://github.com/oschwengers/bakta) (provides `bakta_proteins`)
- [Biopython](https://biopython.org/)

## Usage

### 1. Prepare a samples file

Create a plain-text file (e.g. `samples.txt`) with one genome assembly FASTA
path per line:

```
/path/to/sample1.fasta
/path/to/sample2.fasta
/path/to/sample3.fasta
```

### 2. Edit `config/config.yaml`

Set `samples` to your samples file, `output_dir` to the desired results
directory, and supply any extra tool arguments via the `cli_args` fields.
Set `bakta.db` to your local Bakta database path.

### 3. Run the pipeline

```bash
snakemake --cores <N>
```

To perform a dry-run first:

```bash
snakemake --cores <N> -n
```

## Configuration

All tool-specific parameters are controlled via `config/config.yaml`.

| Key | Description |
|-----|-------------|
| `samples` | Path to text file listing genome assembly FASTA paths (one per line) |
| `output_dir` | Base output directory |
| `ggcaller.cli_args` | CLI arguments passed verbatim to `ggcaller` |
| `panaroo.cli_args` | CLI arguments passed verbatim to `panaroo` |
| `bakta.db` | Path to Bakta database directory |
| `bakta.cli_args` | CLI arguments passed verbatim to `bakta_proteins` |

### Example ggCaller arguments

```yaml
ggcaller:
  cli_args: "--annotation sensitive --graph --threads 16"
```

### Example Panaroo arguments

```yaml
panaroo:
  cli_args: "--clean-mode moderate -a core --core_threshold 0.98 --threads 16"
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
