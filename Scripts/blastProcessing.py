####################################################################################################
# This script:
# 1. Adds reference gene sequences to gene-specific alignment FASTA files
# 2. Parses BLAST tabular output files
# 3. Extracts matching gene regions from plastid genomes
# 4. Handles strand orientation correctly
# 5. Appends extracted sequences to alignment files
####################################################################################################

# ----------------------------------------------------------
# Import packages
# ----------------------------------------------------------

import os
import re
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict 

# ---------------------------------------------------------------------------------------------------
# Command line arguments
# ---------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract gene sequences from plastid genomes using blast results "
            "and append them to gene-specific alignment FASTA files."))

    parser.add_argument("--blast-dir", default="Blast/Results", help="Directory containing blast result subfolders")

    parser.add_argument("--reference-dir", default="Blast/ReferenceGeneSequences", help="Directory containing reference gene FASTA files")

    parser.add_argument("--genome-dir", default="Blast/PlastidSequences", help="Directory containing test plastid genome FASTA files")

    parser.add_argument("--output-dir", required=True, help="Directory for output alignment FASTA files")

    parser.add_argument(
        "--flanking-region",
        type=int,
        default=0,
        help="Number of base pairs to extract on each side of blast hit (default: 0)")

    parser.add_argument(
        "--ir-cutoff",
        type=int,
        default=5000,
        help="Maximum expected gene length; larger hits are flagged (default: 5000)")

    return parser.parse_args()


args = parse_args()

# ---------------------------------------------------------------------------------------------------
# Get gene list from reference FASTA filenames
# ---------------------------------------------------------------------------------------------------

# Each reference gene FASTA file is named <gene>.fasta
# This automatically generates the list of reference genes to process
geneList = []
for f in os.listdir(args.reference_dir):
    if f.endswith(".fasta"):
        gene_name = os.path.splitext(f)[0]  # Strip ".fasta"
        geneList.append(gene_name.lower())  # Lowercase for consistent matching



# ---------------------------------------------------------------------------------------------------
# Add one copy of each reference sequence to the alignment file
# ---------------------------------------------------------------------------------------------------

os.makedirs(args.output_dir, exist_ok=True)

for gene in geneList:
    alignment_path = os.path.join(args.output_dir, f"{gene}-alignment-unaligned.fasta")
    reference_fasta = os.path.join(args.reference_dir, f"{gene}.fasta")
 
    # Only add reference sequences if the alignment file does not already exist
    if os.path.exists(reference_fasta) and not os.path.exists(alignment_path):
        with open(alignment_path, "w") as alignment_file:
            for record in SeqIO.parse(reference_fasta, "fasta"):
                SeqIO.write(record, alignment_file, "fasta")


# ---------------------------------------------------------------------------------------------------
# Map genome accessions to FASTA file paths
# ---------------------------------------------------------------------------------------------------

genome_files = {}
for f in os.listdir(args.genome_dir):
    if f.endswith(".fasta"):
        name_no_ext = os.path.splitext(f)[0]  # remove ".fasta"
        genome_files[name_no_ext] = os.path.join(args.genome_dir, f)

# ---------------------------------------------------------------------------------------------------
# Process blast results
# ---------------------------------------------------------------------------------------------------

bigNoHitsList = []
bigProblemFilesList = []

blast_folders = [
    os.path.join(args.blast_dir, d)
    for d in os.listdir(args.blast_dir)
    if os.path.isdir(os.path.join(args.blast_dir, d))
]

for directory in blast_folders:

    somethingWrongWithTheseFiles = []
    noHitsFiles = []

    for filename in os.listdir(directory):
        if not filename.endswith(".txt"):
            continue

        blastResults = []

        # Read BLAST output lines
        with open(os.path.join(directory, filename)) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                blastResults.append(line)

        # Skip files with no hits
        if not blastResults:
            noHitsFiles.append(os.path.join(os.path.basename(directory), filename))
            continue

        # Split each line into columns
        tabularResults = [line.split("\t") for line in blastResults]

        geneName = os.path.splitext(filename)[0].lower()

		# --------------------------------------------------------------------------------------
        # The blast file may contain hits from multiple genomes - process them independently
        # --------------------------------------------------------------------------------------

        hits_by_genome = defaultdict(list)

        for row in tabularResults:
            genomeID = row[1]  # qseqid = plastid genome
            hits_by_genome[genomeID].append(row)

        # --------------------------------------------------------------------------------------
        # Loop over each genome in the blast file
        # --------------------------------------------------------------------------------------
        for genomeID, genome_hits in hits_by_genome.items():

            # Clean genome ID from blast formatting (gb|XXXX.1|) and remove version numbers
            if "|" in genomeID:
                genomeID = genomeID.split("|")[-2]
            genomeID = genomeID.split(".")[0]

            # Skip if genome FASTA not found
            if genomeID not in genome_files:
                print(f"WARNING: Genome file for {genomeID} not found in {args.genome_dir}. Skipping.")
                continue

            # Extract reference ID (sseqid) from first hit
            refSeqID = genome_hits[0][1]

            # Initialize lists for parsing
            geneBounds = []
            genomeBounds = []
            strands = []

            # Parse each BLAST hit for THIS genome only
            for row in genome_hits:
                q_start = int(row[3])
                q_end   = int(row[4])
                s_start = int(row[5])
                s_end   = int(row[6])

                geneBounds.append([q_start, q_end])
                genomeBounds.append([s_start, s_end])

                strand = "+" if s_start < s_end else "-"
                strands.append(strand)

            # -------------------------------------------------
            # Handle duplicated hits
            # -------------------------------------------------
            sorter = {}
            for q, s, strand in zip(geneBounds, genomeBounds, strands):
                q_tuple = tuple(q)
                if q_tuple not in sorter:
                    sorter[q_tuple] = [(s, strand)]
                else:
                    sorter[q_tuple].append((s, strand))

            duplicateGeneBounds = [[]]
            duplicateGenomeBounds = [[]]

            for k, v in sorter.items():
                for i, (coords, strand) in enumerate(v):
                    if i >= len(duplicateGeneBounds):
                        duplicateGeneBounds.append([])
                        duplicateGenomeBounds.append([])
                    duplicateGeneBounds[i].append(list(k))
                    duplicateGenomeBounds[i].append((coords, strand))

            genomeBoundList = []

            for numberSet in range(len(duplicateGeneBounds)):
                gene_array = np.array(duplicateGeneBounds[numberSet])
                genome_array = np.array([x[0] for x in duplicateGenomeBounds[numberSet]])
                strand_list = [x[1] for x in duplicateGenomeBounds[numberSet]]

                min_idx = np.unravel_index(np.argmin(gene_array), gene_array.shape)
                max_idx = np.unravel_index(np.argmax(gene_array), gene_array.shape)

                s1 = genome_array[min_idx]
                s2 = genome_array[max_idx]

                start = min(s1, s2)
                end = max(s1, s2)
                strand = strand_list[min_idx[0]]

                genomeBoundList.append((start, end, strand))

                if end - start > args.ir_cutoff:
                    somethingWrongWithTheseFiles.append(os.path.join(os.path.basename(directory), filename))

            # -------------------------------------------------
            # IR duplicate hit check
            # -------------------------------------------------
            if len(genomeBoundList) == 2:
                (start1, end1, strand1), (start2, end2, strand2) = genomeBoundList
                len1 = end1 - start1
                len2 = end2 - start2
                if abs(len1 - len2) < 30:
                    # Keep only one copy if the hits are almost identical (this will happen if they are in the IR)
                    genomeBoundList = [genomeBoundList[0]]

            # -------------------------------------------------
            # Read genome fasta
            # -------------------------------------------------
            genome_record = SeqIO.read(genome_files[genomeID], "fasta")
            genome_seq = genome_record.seq
            
            # -------------------------------------------------
            # Extract species name from FASTA header
            # -------------------------------------------------            
            
            # Expected format: >KU588419.1 Orthilia secunda chloroplast genome ...
            
            # split the parts of the name up by spaces
            description_parts = genome_record.description.split()
            # assign species name in case parsing does not work
            speciesName = "Unknown_species"
			
			# get the second and third items in the list to make the species name
            if len(description_parts) >= 3:
                genus = description_parts[1]
                species = description_parts[2]
                # remove non-letter characters from species epithet
                species = re.sub(r"[^A-Za-z]", "", species)
                speciesName = f"{genus}_{species}"
            

            # -------------------------------------------------
            # Extract sequences and append to alignment files
            # -------------------------------------------------
            
            # make a set to avoid sequences identical in name and content being written out more than once
            written_seqs = set()
            
            for start, end, strand in genomeBoundList:
                start0 = max(start - args.flanking_region - 1, 0)
                end0 = min(end + args.flanking_region, len(genome_seq))

                subseq = genome_seq[start0:end0]

                if strand == "-":
                    subseq = subseq.reverse_complement()

                header = (
                    f"{genome_record.id}|{speciesName}|{geneName}|"
                    f"{start}-{end}|refSeq:{refSeqID}"
                )
                
                # convert DNA sequence to a string
                seq_str = str(subseq)
                
                # do not write out sequences that have already been written out
                key = (header, seq_str)
                if key in written_seqs:
                    continue
                written_seqs.add(key)

                new_record = SeqRecord(subseq, id=header, description="")

                aln_path = os.path.join(
                    args.output_dir, f"{geneName}-alignment-unaligned.fasta"
                )

                with open(aln_path, "a") as out:
                    SeqIO.write(new_record, out, "fasta")

        # End of genome loop

    # Collect problem/no-hit files
    bigProblemFilesList.append(somethingWrongWithTheseFiles)
    bigNoHitsList.append(noHitsFiles)

# -----------------------------------------------------------------------------------------------------------------------
# Write summary files for manual inspection
# -----------------------------------------------------------------------------------------------------------------------
problem_file_path = os.path.join(args.output_dir, "FilesToCheckAgain.txt")
nohits_file_path = os.path.join(args.output_dir, "NoHitsFiles.txt")

with open(problem_file_path, "w") as out:
    for f in sorted(set(sum(bigProblemFilesList, []))):
        out.write(f + "\n")

with open(nohits_file_path, "w") as out:
    for f in sorted(set(sum(bigNoHitsList, []))):
        out.write(f + "\n")
