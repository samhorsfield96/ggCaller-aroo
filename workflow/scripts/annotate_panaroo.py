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

from Bio import SeqIO

logging.basicConfig(
    filename=str(snakemake.log[0]),
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)

bakta_tsv   = Path(str(snakemake.input.bakta_tsv))
panaroo_dir = Path(str(snakemake.params.panaroo_dir))

out_gpa       = Path(str(snakemake.output.gpa))
out_gpa_roary = Path(str(snakemake.output.gpa_roary))
out_gene_data = Path(str(snakemake.output.gene_data))
out_pan_ref   = Path(str(snakemake.output.pan_ref))
out_dna_cds   = Path(str(snakemake.output.dna_cds))
out_prot_cds  = Path(str(snakemake.output.prot_cds))
out_gml       = Path(str(snakemake.output.gml))

out_gpa.parent.mkdir(parents=True, exist_ok=True)


# ── 1. Load bakta annotations ─────────────────────────────────────────────────

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


annotations = load_bakta_annotations(bakta_tsv)
logging.info("Loaded bakta annotations for %d groups.", len(annotations))


# ── 2. Parse GML: seqID → group_name ─────────────────────────────────────────

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


seqid_to_group = parse_gml_seqid_to_group(panaroo_dir / "final_graph.gml")
logging.info("Parsed %d seqID→group mappings from GML.", len(seqid_to_group))


# ── Helper: build FASTA annotation string ────────────────────────────────────

def fasta_annotation(ann):
    """Return '[gene=X] [product=Y]' annotation string (omits empty fields)."""
    parts = []
    if ann.get("gene"):
        parts.append(f"[gene={ann['gene']}]")
    if ann.get("product"):
        parts.append(f"[product={ann['product']}]")
    return " ".join(parts)


def write_fasta(fh, record_id, annotation_str, sequence):
    """Write a single FASTA record with 60-char wrapped sequence."""
    header = f">{record_id} {annotation_str}" if annotation_str else f">{record_id}"
    fh.write(header + "\n")
    seq = str(sequence)
    for i in range(0, len(seq), 60):
        fh.write(seq[i : i + 60] + "\n")


# ── 3. pan_genome_reference.fa ───────────────────────────────────────────────

written = 0
with open(out_pan_ref, "w") as fh:
    for record in SeqIO.parse(panaroo_dir / "pan_genome_reference.fa", "fasta"):
        ann  = annotations.get(record.id, {})
        desc = fasta_annotation(ann)
        write_fasta(fh, record.id, desc, record.seq)
        written += 1
logging.info("pan_genome_reference.fa: updated %d records.", written)


# ── 4. combined_DNA_CDS.fasta & combined_protein_CDS.fasta ──────────────────

def update_cds_fasta(in_path, out_path):
    written = 0
    with open(out_path, "w") as fh:
        for record in SeqIO.parse(in_path, "fasta"):
            group = seqid_to_group.get(record.id, "")
            ann   = annotations.get(group, {}) if group else {}
            desc  = fasta_annotation(ann)
            write_fasta(fh, record.id, desc, record.seq)
            written += 1
    logging.info("%s: updated %d records.", out_path.name, written)


update_cds_fasta(panaroo_dir / "combined_DNA_CDS.fasta",     out_dna_cds)
update_cds_fasta(panaroo_dir / "combined_protein_CDS.fasta", out_prot_cds)


# ── 5. gene_data.csv ─────────────────────────────────────────────────────────

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

def update_gpa(in_path, out_path):
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


update_gpa(
    panaroo_dir / "gene_presence_absence.csv",       out_gpa
)
update_gpa(
    panaroo_dir / "gene_presence_absence_roary.csv",  out_gpa_roary
)


# ── 7. final_graph.gml ───────────────────────────────────────────────────────

def gml_escape(s):
    """Escape backslashes and double-quotes for GML string values."""
    return s.replace("\\", "\\\\").replace('"', '\\"')


def update_gml(in_path, out_path):
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


update_gml(panaroo_dir / "final_graph.gml", out_gml)

logging.info("Annotation update complete.")
