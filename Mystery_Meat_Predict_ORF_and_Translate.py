# Import libraries
import pandas as pd

# Reading in a Fasta File
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sequences = list(SeqIO.parse("combined_reps.fasta", "fasta"))

# List to store output amino-acid SeqRecord objects
aa_sequences = []

# Loop through each nucleotide sequence in FASTA
for record in sequences:
    # Extracting the nucleotide sequence from FASTA record

    dna_seq = record.seq

    # Create forward and reverse-complement strands
    strands = [
        (dna_seq, "forward"),
        (dna_seq.reverse_complement(), "reverse")
    ]

    best_orf = ""          # Will store the longest ORF found across all frames
    # Will store which frame (0,1,2) contains the longest ORF
    best_frame = None
    best_strand = None

    # Examine both strands
    for strand_seq, strand_name in strands:

        # Translate sequence in all three forward reading frames
        for frame in range(3):

            # Translate starting at this frame; keep stop codons (to_stop=False)
            protein = strand_seq[frame:].translate(to_stop=False)

            # Split translation into fragments separated by stop codons ("*")
            # Each fragment is a potential ORF (no internal stops)
            fragments = str(protein).split("*")

            # Find the longest continuous amino-acid stretch (longest ORF) in this frame
            longest_in_frame = max(fragments, key=len)

            # If this ORF is longer than the best one we've seen so far, store it
            if len(longest_in_frame) > len(best_orf):
                best_orf = longest_in_frame
                best_frame = frame
                best_strand = strand_name

    # Convert the longest ORF string back into a Seq object
    best_protein = Seq(best_orf)

    # Create a new SeqRecord containing the longest ORF for this sequence
    aa_record = SeqRecord(
        best_protein,
        id=record.id,                      # Keep the same FASTA ID
        description=(
            f"{record.description} "
            f"strand={best_strand} "
            f"frame={best_frame}"
        ),
        annotations={
            "type": "longest ORF",
            "frame": best_frame,           # Which reading frame the ORF came from
            "strand": best_strand
        }
    )

    aa_sequences.append(aa_record)         # Add to output list

# Write all translated longest-ORF sequences to a FASTA file
SeqIO.write(aa_sequences, "translated.fas", "fasta")
