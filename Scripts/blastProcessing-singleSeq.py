####################################################################################################
# This script:
# 1. Adds reference gene sequences to gene-specific alignment FASTA files
# 2. Parses BLAST tabular output files
# 3. Merges multiple hits into gene regions
# 4. Extracts matching gene regions from plastid genomes
# 5. Handles strand orientation correctly
# 6. Appends extracted sequences to alignment files
####################################################################################################

# ---------------------------------------------------------------------------------------------------
# Import packages
# ---------------------------------------------------------------------------------------------------


import os
import re
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------------------------------
# Command line arguments
# ---------------------------------------------------------------------------------------------------

# whole IR handling section/input option is totally redundant because this only returns one hit for each gene/sequence combo
# get rid of this at a later date


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract gene sequences from plastid genomes using blast results "
            "and append them to gene-specific alignment FASTA files."))

    # Directory containing BLAST result subfolders
    parser.add_argument("--blast-dir", default="Blast/Results")

    # Directory containing reference gene FASTA files
    parser.add_argument("--reference-dir", default="Blast/ReferenceGeneSequences")

    # Directory containing plastid genome FASTA files
    parser.add_argument("--genome-dir", default="Blast/PlastidSequences")

    # Output directory for gene-specific alignment FASTA files
    parser.add_argument("--output-dir", required=True)

    # Optional flanking region to extract around BLAST hits
    parser.add_argument("--flanking-region", type=int, default=0,
                        help="Number of bases to extract upstream and downstream of merged hit.")

    # Maximum allowed genomic distance between BLAST hits when merging
    parser.add_argument("--merge-gap", type=int, default=1000,
                        help="Max genomic gap allowed when merging HSPs (default=1000).")


    # Allowed length difference when determining IR duplicates
    parser.add_argument("--ir-length-tolerance", type=int, default=50,
                        help="Length difference tolerance for IR duplicate removal.")

    return parser.parse_args()


args = parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)


# ---------------------------------------------------------------------------------------------------
# Get list of genes from reference FASTA files
# ---------------------------------------------------------------------------------------------------

# Each reference gene FASTA file is named <gene>.fasta
# This automatically generates the list of reference genes to process
geneList = [os.path.splitext(f)[0].lower()
    for f in os.listdir(args.reference_dir)
    if f.endswith(".fasta")]


# ---------------------------------------------------------------------------------------------------
# Initialize alignment files with reference sequences
# ---------------------------------------------------------------------------------------------------

# For each gene, create an alignment file if it does not already exist.
# We write the reference sequence first so that all test sequences
# will later be appended below it.

for gene in geneList:
    
    # make a fasta file containing just reference sequences
    # make a fasta file for reference sequences and "test" sequences
    aln_path = os.path.join(args.output_dir, f"{gene}-alignment-unaligned.fasta")
    ref_path = os.path.join(args.reference_dir, f"{gene}.fasta")

    # Only create alignment file if it does not already exist
    if os.path.exists(ref_path) and not os.path.exists(aln_path):
        # alignment files are being written not appended
        with open(aln_path, "w") as out:
            # parse the fasta file and write out sequence
            for record in SeqIO.parse(ref_path, "fasta"):
                SeqIO.write(record, out, "fasta")


# ---------------------------------------------------------------------------------------------------
# Map genome accession IDs to genome FASTA file paths
# ---------------------------------------------------------------------------------------------------

# make a dictionary where keys = genbank ID (without decimal) and values = full file path
# means blast hits can be linked to the correct file for sequence extraction
genome_files = {os.path.splitext(f)[0]: os.path.join(args.genome_dir, f)
    for f in os.listdir(args.genome_dir)
    if f.endswith(".fasta")}


# ---------------------------------------------------------------------------------------------------
# Process BLAST result folders
# ---------------------------------------------------------------------------------------------------

# blast outputs are in blast directory specified on command line
blast_folders = [os.path.join(args.blast_dir, d)
    for d in os.listdir(args.blast_dir)
    if os.path.isdir(os.path.join(args.blast_dir, d))]

# loop through the blast results by subfolder (organised by accession)
for directory in blast_folders:

    for filename in os.listdir(directory):

        # only process tabular blast results files
        if not filename.endswith(".txt"):
            continue

        filepath = os.path.join(directory, filename)

        # Read BLAST tabular lines
        with open(filepath) as fh:
            lines = [l.strip().split("\t") for l in fh if l.strip()]

        # Skip if no hits
        if not lines:
            continue
        # Gene name is inferred from BLAST filename
        geneName = os.path.splitext(filename)[0].lower()

        # Group BLAST hits by genome (each BLAST file may contain multiple genomes)
        hits_by_genome = defaultdict(list)
        for row in lines:
            genomeID = row[1]   # qseqid = plastid genome
            hits_by_genome[genomeID].append(row)

        # Process each genome independently
        for genomeID_raw, genome_hits in hits_by_genome.items():

            # Clean genome accession formatting
            if "|" in genomeID_raw:
                genomeID = genomeID_raw.split("|")[-2]
            else:
                genomeID = genomeID_raw

            genomeID = genomeID.split(".")[0]

            if genomeID not in genome_files:
                print(f"WARNING: Genome file for {genomeID} not found.")
                continue

            # --------------------------------------------------------------------------------------
            # Parse BLAST hits into structured dictionaries
            # --------------------------------------------------------------------------------------

            hits = []

            for row in genome_hits:

                # Query coordinates (reference gene)
                q_start = int(row[3])
                q_end   = int(row[4])

                # Subject coordinates (plastid genome)
                s_start = int(row[5])
                s_end   = int(row[6])

                # Determine strand orientation
                strand = "+" if s_start < s_end else "-"

                hits.append({
                    "q_start": min(q_start, q_end),
                    "q_end": max(q_start, q_end),
                    "s_start": min(s_start, s_end),
                    "s_end": max(s_start, s_end),
                    "strand": strand})

            # --------------------------------------------------------------------------------------
            # Merge hits by strand
            # --------------------------------------------------------------------------------------

            merged_regions = []

            for strand in ["+", "-"]:

                # filter hits by strand
                strand_hits = [h for h in hits if h["strand"] == strand]
                if not strand_hits:
                    continue

                # Sort hits along genome by genome start
                strand_hits.sort(key=lambda x: x["s_start"])
                current = strand_hits[0]

                # loop through hits to merge hits that are close on the same strand
                for next_hit in strand_hits[1:]:
                    genomic_gap = next_hit["s_start"] - current["s_end"]
                    query_gap   = next_hit["q_start"] - current["q_end"]

                    # Merge if the hits are close and reasonably contigunous
                    if genomic_gap <= args.merge_gap and query_gap >= -100:
                        current["s_end"] = max(current["s_end"], next_hit["s_end"])
                        current["q_end"] = max(current["q_end"], next_hit["q_end"])
                    else:
                        # save the current region and start again
                        merged_regions.append(current)
                        current = next_hit

                # append the last current region
                merged_regions.append(current)

            # --------------------------------------------------------------------------------------
            # Remove IR duplicates
            # --------------------------------------------------------------------------------------

            # if merged regions appear to be very similar then only keep one
            if len(merged_regions) > 1:
                lengths = [r["s_end"] - r["s_start"] for r in merged_regions]
                if max(lengths) - min(lengths) < args.ir_length_tolerance:
                    merged_regions = [merged_regions[0]]

            # --------------------------------------------------------------------------------------
            # Remove remaining duplicates
            # --------------------------------------------------------------------------------------

            # if there is still more than one hit per genome, keep the longest one
            if len(merged_regions) > 1:
                merged_regions.sort(
                    key=lambda r: (r["s_end"] - r["s_start"]),
                    reverse=True
                )
                merged_regions = [merged_regions[0]]

            # --------------------------------------------------------------------------------------
            # Extract sequences from genome
            # --------------------------------------------------------------------------------------

            genome_record = SeqIO.read(genome_files[genomeID], "fasta")
            genome_seq = genome_record.seq

            # extract species name from FASTA header
            # split fasta header by space
            description_parts = genome_record.description.split()
            # give default name of unknown species in case header can't be parsed
            speciesName = "Unknown_species"

            # usually fasta file headers from genbank go *accession*space*genus*space*specificEpithet*
            if len(description_parts) >= 3:
                genus = description_parts[1]
                species = re.sub(r"[^A-Za-z]", "", description_parts[2])
                speciesName = f"{genus}_{species}"

            # loop through regions and extract the sequences
            for region in merged_regions:
                start = max(region["s_start"] - args.flanking_region - 1, 0)
                end   = min(region["s_end"] + args.flanking_region, len(genome_seq))
                subseq = genome_seq[start:end]

                # reverse complement sequences on the minus strand
                if region["strand"] == "-":
                    subseq = subseq.reverse_complement()

                # write out the fasta header
                header = (f"{genome_record.id}|{speciesName}|{geneName}|"
                    f"{region['s_start']}-{region['s_end']}")

                # define what the new sequence will be (header and best merged hit)
                new_record = SeqRecord(subseq, id=header, description="")

                # append sequence to the appropriate gene's multifasta file
                aln_path = os.path.join(args.output_dir,
                    f"{geneName}-alignment-unaligned.fasta")
                with open(aln_path, "a") as out:
                    SeqIO.write(new_record, out, "fasta")


print("Finished.")