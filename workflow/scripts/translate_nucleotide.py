"""Translate nucleotide FASTA sequences to amino acid sequences.

Used to convert Panaroo representative gene sequences (pan_genome_reference.fa)
to protein sequences for annotation with bakta_proteins.
"""

import logging

from Bio import SeqIO

logging.basicConfig(
    filename=str(snakemake.log[0]),
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)

translated_count = 0
skipped_count = 0

with open(str(snakemake.output[0]), "w") as out_fh:
    for record in SeqIO.parse(str(snakemake.input[0]), "fasta"):
        seq = record.seq

        # Trim to a multiple of 3 if the sequence length is not exact
        remainder = len(seq) % 3
        if remainder != 0:
            logging.warning(
                f"Sequence {record.id} length ({len(seq)}) is not a multiple of 3; "
                f"trimming {remainder} base(s) from the end."
            )
            seq = seq[:-remainder]

        protein = seq.translate(to_stop=True)

        if len(protein) == 0:
            logging.warning(
                f"Sequence {record.id} translated to an empty protein; skipping."
            )
            skipped_count += 1
            continue

        out_fh.write(f">{record.id}\n{protein}\n")
        translated_count += 1

logging.info(
    f"Finished: translated {translated_count} sequences, skipped {skipped_count}."
)
