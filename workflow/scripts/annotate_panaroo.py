"""Propagate bakta_proteins functional annotations onto all panaroo output files.

Updated files written to the annotated output directory:
  - pan_genome_reference.fa          (FASTA headers gain [gene=] [product=])
  - combined_DNA_CDS.fasta           (same, indexed by seqID via GML)
  - combined_protein_CDS.fasta       (same)
  - gene_data.csv                    (gene_name and description columns filled)
  - gene_presence_absence.csv        (Non-unique Gene name and Annotation filled)
  - gene_presence_absence_roary.csv  (same)
  - final_graph.gml                  (annotation and description fields filled)
"""

import csv
import logging
import re
from pathlib import Path
import argparse

from Bio import SeqIO

# ── Helper: load bakta annotations ────────────────────────────────────

def load_bakta_annotations(path):
    """Parse bakta_proteins TSV → {group_id: {'gene': str, 'product': str}}."""
    result = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0] == "ID":
                continue
            group_id = parts[0]
            gene     = parts[2].strip() if len(parts) > 2 else ""
            product  = parts[3].strip() if len(parts) > 3 else ""
            result[group_id] = {"gene": gene, "product": product}
    return result

# ── Helper: parse GML to get group IDs ────────────────────────────────────

def parse_gml_seqid_to_group(gml_path):
    """Parse panaroo final_graph.gml → {seqID: group_name}.

    Each GML node stores its cluster members via repeated 'seqIDs' keys and its
    group label via the 'name' key (added by panaroo, always after seqIDs).
    """
    seqid_to_group   = {}
    bracket_depth    = 0
    node_start_depth = None
    current_seqids   = []
    current_name     = None

    with open(gml_path) as fh:
        for line in fh:
            stripped = line.strip()

            if stripped.endswith("["):
                bracket_depth += 1
                if stripped == "node [":
                    node_start_depth = bracket_depth
                    current_seqids   = []
                    current_name     = None

            elif stripped == "]":
                if (
                    node_start_depth is not None
                    and bracket_depth == node_start_depth
                ):
                    if current_name:
                        for sid in current_seqids:
                            seqid_to_group[sid] = current_name
                    node_start_depth = None
                bracket_depth -= 1

            elif node_start_depth is not None and bracket_depth == node_start_depth:
                m = re.match(r'seqIDs\s+"([^"]+)"', stripped)
                if m and m.group(1) != "_networkx_list_start":
                    current_seqids.append(m.group(1))
                m = re.match(r'name\s+"([^"]+)"', stripped)
                if m:
                    current_name = m.group(1)

    return seqid_to_group

# ── Helper: build FASTA annotation string ────────────────────────────────────

def fasta_annotation(ann):
    """Return '[gene=X] [product=Y]' annotation string (omits empty fields)."""
    parts = []
    if ann.get("gene"):
        parts.append(f"[gene={ann['gene']}]")
    if ann.get("product"):
        parts.append(f"[product={ann['product']}]")
    return " ".join(parts)

# ── Helper: write FASTA record ───────────────────────────────────────────────

def write_fasta(fh, record_id, annotation_str, sequence):
    """Write a single FASTA record with 60-char wrapped sequence."""
    header = f">{record_id} {annotation_str}" if annotation_str else f">{record_id}"
    fh.write(header + "\n")
    seq = str(sequence)
    for i in range(0, len(seq), 60):
        fh.write(seq[i : i + 60] + "\n")

# ── Main annotation update functions ─────────────────────────────────────────

def update_cds_fasta(in_path, out_path, logging=logging):
    written = 0
    with open(out_path, "w") as fh:
        for record in SeqIO.parse(in_path, "fasta"):
            group = seqid_to_group.get(record.id, "")
            ann   = annotations.get(group, {}) if group else {}
            desc  = fasta_annotation(ann)
            write_fasta(fh, record.id, desc, record.seq)
            written += 1
    logging.info("%s: updated %d records.", out_path.name, written)

# ── Update gene presence/absence ──────────────────────────────────────────

def update_gpa(in_path, out_path, logging=logging):
    updated = 0
    with (
        open(in_path, newline="") as in_fh,
        open(out_path, "w", newline="") as out_fh,
    ):
        reader = csv.DictReader(in_fh)
        writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            ann = annotations.get(row["Gene"], {})
            if ann:
                row["Non-unique Gene name"] = ann.get(
                    "gene", row.get("Non-unique Gene name", "")
                )
                row["Annotation"] = ann.get(
                    "product", row.get("Annotation", "")
                )
                updated += 1
            writer.writerow(row)
    logging.info("%s: updated %d rows.", out_path.name, updated)


# ── Update GML annotations ───────────────────────────────────────────────

def gml_escape(s):
    """Escape backslashes and double-quotes for GML string values."""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def update_gml(in_path, out_path, logging=logging):
    """Buffer each node block, then flush with updated annotation/description."""
    bracket_depth    = 0
    node_start_depth = None
    buffer           = []
    current_name     = None
    updated_nodes    = 0

    _ann_re  = re.compile(r'^\s*annotation\s+"[^"]*"')
    _desc_re = re.compile(r'^\s*description\s+"[^"]*"')
    _name_re = re.compile(r'^\s*name\s+"([^"]+)"')

    with open(in_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            stripped = line.strip()

            if stripped.endswith("["):
                bracket_depth += 1
                if stripped == "node [":
                    node_start_depth = bracket_depth
                    buffer           = [line]
                    current_name     = None
                elif node_start_depth is not None and bracket_depth > node_start_depth:
                    buffer.append(line)
                else:
                    out_fh.write(line)

            elif stripped == "]":
                if (
                    node_start_depth is not None
                    and bracket_depth == node_start_depth
                ):
                    # Flush node buffer with annotation substitutions
                    buffer.append(line)
                    ann     = annotations.get(current_name, {})
                    gene    = gml_escape(ann.get("gene",    ""))
                    product = gml_escape(ann.get("product", ""))
                    for buf_line in buffer:
                        if _ann_re.match(buf_line):
                            indent = re.match(r"^(\s*)", buf_line).group(1)
                            out_fh.write(f'{indent}annotation "{gene}"\n')
                        elif _desc_re.match(buf_line):
                            indent = re.match(r"^(\s*)", buf_line).group(1)
                            out_fh.write(f'{indent}description "{product}"\n')
                        else:
                            out_fh.write(buf_line)
                    if ann:
                        updated_nodes += 1
                    node_start_depth = None
                    current_name     = None
                    buffer           = []
                elif node_start_depth is not None and bracket_depth > node_start_depth:
                    buffer.append(line)
                else:
                    out_fh.write(line)
                bracket_depth -= 1

            else:
                if node_start_depth is not None:
                    m = _name_re.match(line)
                    if m:
                        current_name = m.group(1)
                    buffer.append(line)
                else:
                    out_fh.write(line)

    logging.info("final_graph.gml: updated %d nodes.", updated_nodes)

def get_options():
    description = "Annotates tokenised gene clusters"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python annotate_tokens.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--bakta_tsv',
                    required=True,
                    help='Path to bakta TSV file.')
    IO.add_argument('--panaroo_dir',
            required=True,
            help='Path to panaroo output directory.')
    IO.add_argument('--outdir',
                required=True,
                help='Output directory.')
    IO.add_argument('--output_type',
                required=True,
                choices=['gpa', 'gpa_roary', 'gene_data', 'pan_ref', 'dna_cds', 'prot_cds', 'gml'],
                help='Output type.')
    IO.add_argument('--logfile',
                required=True,
                help='Path to log file.')
    return parser.parse_args()

def main():
    options = get_options()
    bakta_tsv   = Path(options.bakta_tsv)
    panaroo_dir = Path(options.panaroo_dir)
    outdir       = Path(options.outdir)
    logfile      = Path(options.logfile)

    logging.basicConfig(
        filename=str(logfile),
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    out_gpa       = outdir / "gene_presence_absence.csv" if options.output_type == 'gpa' else None
    out_gpa_roary = outdir / "gene_presence_absence_roary.csv" if options.output_type == 'gpa_roary' else None
    out_gene_data = outdir / "gene_data.csv" if options.output_type == 'gene_data' else None
    out_pan_ref   = outdir / "pan_genome_reference.fa" if options.output_type == 'pan_ref' else None
    out_dna_cds   = outdir / "combined_DNA_CDS.fasta" if options.output_type == 'dna_cds' else None
    out_prot_cds  = outdir / "combined_protein_CDS.fasta" if options.output_type == 'prot_cds' else None
    out_gml       = outdir / "final_graph.gml" if options.output_type == 'gml' else None

    outdir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load bakta annotations ─────────────────────────────────────────────────
    annotations = load_bakta_annotations(bakta_tsv)
    logging.info("Loaded bakta annotations for %d groups.", len(annotations))

    # ── 2. Parse GML: seqID → group_name ─────────────────────────────────────────
    seqid_to_group = parse_gml_seqid_to_group(panaroo_dir / "final_graph.gml")
    logging.info("Parsed %d seqID→group mappings from GML.", len(seqid_to_group))

    # ── 3. pan_genome_reference.fa ───────────────────────────────────────────────

    if out_pan_ref:
        written = 0
        with open(out_pan_ref, "w") as fh:
            for record in SeqIO.parse(panaroo_dir / "pan_genome_reference.fa", "fasta"):
                ann  = annotations.get(record.id, {})
                desc = fasta_annotation(ann)
                write_fasta(fh, record.id, desc, record.seq)
                written += 1
        logging.info("pan_genome_reference.fa: updated %d records.", written)


    # ── 4. combined_DNA_CDS.fasta & combined_protein_CDS.fasta ──────────────────

    if out_dna_cds:
        update_cds_fasta(panaroo_dir / "combined_DNA_CDS.fasta", out_dna_cds, logging=logging)
    if out_prot_cds:
        update_cds_fasta(panaroo_dir / "combined_protein_CDS.fasta", out_prot_cds, logging=logging)

    # ── 5. gene_data.csv ─────────────────────────────────────────────────────────

    if out_gene_data:
        updated = 0
        with (
            open(panaroo_dir / "gene_data.csv", newline="") as in_fh,
            open(out_gene_data, "w", newline="") as out_fh,
        ):
            reader = csv.DictReader(in_fh)
            writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                group = seqid_to_group.get(row["clustering_id"], "")
                ann   = annotations.get(group, {}) if group else {}
                if ann:
                    row["gene_name"]   = ann.get("gene",    row.get("gene_name",   ""))
                    row["description"] = ann.get("product", row.get("description", ""))
                    updated += 1
                writer.writerow(row)
        logging.info("gene_data.csv: updated %d rows.", updated)


    # ── 6. gene_presence_absence.csv & gene_presence_absence_roary.csv ───────────

    if out_gpa:
        update_gpa(
            panaroo_dir / "gene_presence_absence.csv", out_gpa, logging=logging
        )
    if out_gpa_roary:
        update_gpa(
            panaroo_dir / "gene_presence_absence_roary.csv",  out_gpa_roary, logging=logging
        )

    # ── 7. final_graph.gml ───────────────────────────────────────────────────────

    if out_gml:
        update_gml(panaroo_dir / "final_graph.gml", out_gml, logging=logging)

    logging.info("Annotation update complete.")
