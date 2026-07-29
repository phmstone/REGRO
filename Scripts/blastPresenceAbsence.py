####################################################################################################
# This script:
# 1. Reads the presence/absence TSV produced by the gene profiling script
# 2. Identifies taxa missing functional copies of one or more genes
# 3. Downloads reference GenBank files
# 4. Extracts reference gene sequences from these files
# 5. Downloads plastid fasta files for taxa with missing genes
# 6. Makes BLAST databases for these plastid sequences
# 7. Runs BLAST to recover missing genes
####################################################################################################

# ----------------------------------------------------------
# Import packages
# ----------------------------------------------------------

import argparse
import os
import time
import subprocess
import certifi
import re
from Bio import Entrez, SeqIO

# allows python to interact with internet? Maybe fix this before I send it out
# at least flag that could be an issue for mac users
os.environ["SSL_CERT_FILE"] = certifi.where()

# ----------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------
parser = argparse.ArgumentParser(description="Recover missing plastid genes using BLAST")
parser.add_argument("--input", required=True, help="Presence/absence TSV file")
parser.add_argument("--email", required=True, help="Email address for NCBI")
parser.add_argument("--reference-ids", required=False, default=None, help="Text file of reference GenBank accessions")
parser.add_argument("--blast-type", choices=["blastn", "tblastx"], default="blastn", help="BLAST program to use")
parser.add_argument("--outdir", default="Blast", help="Base output directory (subdirectories will be created and named automatically)")
parser.add_argument("--fastaMode", required=False, help="FASTA file of reference gene sequences. Headers must start with gene name matching TSV (before underscore).")
args = parser.parse_args()

# ----------------------------------------------------------
# Define fixed output directory structure
# ----------------------------------------------------------

# Subdirectories are named automatically but user can name the base output directory
base_outdir = args.outdir
args.reference_outdir = os.path.join(base_outdir, "ReferenceGeneSequences")
args.reference_gbk_dir = os.path.join(base_outdir, "ReferenceGenomes")
args.plastid_fasta_dir = os.path.join(base_outdir, "PlastidSequences")
args.blast_db_dir = os.path.join(base_outdir, "Databases")
args.blast_out = os.path.join(base_outdir, "Results")

# ----------------------------------------------------------
# NCBI set up
# ----------------------------------------------------------

# need an account with NCBI to download things
Entrez.email = args.email

# ----------------------------------------------------------
# Read presence/absence TSV
# ----------------------------------------------------------

# read in the presence/absence TSV file
with open(args.input) as fh:
    lines = [line.rstrip("\n").split("\t") for line in fh]

# separate the header of the TSV file from the presence/absence profiles
# header row contains gene names starting from column 3
header = lines.pop(0)
gene_names = header[2:]  # canonical names from TSV

# make empty lists
latin_names, genbank_ids, presence_absence = [], [], []

# populate lists with relevant fields from TSV
for row in lines:
    latin_names.append(row[0])
    genbank_ids.append(row[1])
    presence_absence.append(row[2:])

print(f"{len(genbank_ids)} taxa total")
print(f"Using {len(gene_names)} genes from TSV header as canonical names")

# ----------------------------------------------------------
# Identify taxa missing genes
# ----------------------------------------------------------
missing_taxa_ids = []
complete_taxa_ids = []

# if all genes are present in  a plastid genome, record that so it can be used as a reference
for i, profile in enumerate(presence_absence):
    if all(value == "0" for value in profile):
        complete_taxa_ids.append(genbank_ids[i])
    else:
        missing_taxa_ids.append(genbank_ids[i])
        
# print summary
print(f"{len(complete_taxa_ids)} taxa with all genes present")
print(f"{len(missing_taxa_ids)} taxa missing ≥1 gene")

# Identify which taxa are missing which genes
missing_by_gene = [[] for _ in range(len(gene_names))]
for taxon_idx, profile in enumerate(presence_absence):
    for gene_idx, value in enumerate(profile):
        if value in {"1", "2"}:  # 1=missing, 2=pseudogene
            missing_by_gene[gene_idx].append(genbank_ids[taxon_idx])


# ----------------------------------------------------------
# FASTA MODE: Use user-provided reference sequences
# ----------------------------------------------------------
if args.fastaMode:
    print("FASTA mode enabled: using user-provided reference sequences")

    os.makedirs(args.reference_outdir, exist_ok=True)

    # track which genes have at least one sequence
    genes_found = set()

    # parse FASTA
    for record in SeqIO.parse(args.fastaMode, "fasta"):
        header = record.id

        # extract gene name (before underscore)
        if "_" in header:
            gene = header.split("_")[0]
        else:
            gene = header

        # validate gene name
        if gene not in gene_names:
            continue

        # validate header format (no forbidden characterss)
        # use the "|" character for parsing blast hits later on down the line
        forbiddenCharacters = {"|"}
        invalidHeader = [c for c in header if c in forbiddenCharacters]
        if invalidHeader:
            raise ValueError(
                f"Invalid FASTA header: {header} \n"
                f"Forbidden characters found: {sorted(set(invalidHeader))}")

        outfile = os.path.join(args.reference_outdir, f"{gene}.fasta")
        with open(outfile, "a") as out:
            out.write(f">{header}\n{str(record.seq)}\n")

        genes_found.add(gene)

    # check all genes are represented
    missing = set(gene_names) - genes_found
    if missing:
        raise ValueError(
            "ERROR: Missing reference sequences for genes:\n" +
            "\n".join(sorted(missing)) + "\n" +
            "Make sure that you have reference IDs that have all of the genes in the gene list annotated as present \n" +
            "If you are having trouble finding such accessions then try using --fastaMode instead")

    print(f"Loaded reference sequences for {len(genes_found)} genes")



# ----------------------------------------------------------
# GENBANK MODE (default)
# ----------------------------------------------------------

if not args.fastaMode:

	# ----------------------------------------------------------
	# Read reference accessions and make new references
	# ----------------------------------------------------------
	
	# empty list for reference IDs to be filled in from user input reference ID .txt file or from .TSV
	reference_ids = []
	
	# Add in reference IDs from the user input reference ID file if present
	if args.reference_ids:
		with open(args.reference_ids) as fh:
			reference_ids = [line.strip() for line in fh if line.strip()]
	
	# Add taxa from the TSV that had all genes present to reference IDs
	reference_ids += complete_taxa_ids
	reference_ids = list(set(reference_ids))  # remove duplicates
	
	# Print error message if the reference ID list is empty
	if not reference_ids:
		raise ValueError("Error: No reference sequences available. \n"
						"Provide a --reference-ids .txt file containing at least taxon with all genes present or remove genes from the gene list so at least one taxon has all genes present. \n"
						"Alternatively try running in --fastaMode using a hybrid reference file. \n"
						"Find instructions for --fastaMode on the GitHub.")
	
	# ----------------------------------------------------------
	# Download reference GenBank files
	# ----------------------------------------------------------
	
	# make the directory for reference genbank files if it doesn't already exist
	os.makedirs(args.reference_gbk_dir, exist_ok=True)
	reference_gbk_files = []
	
	# loop over the accession numbers fo the reference sequences
	for accession in reference_ids:
		gbk_file = os.path.join(args.reference_gbk_dir, f"{accession}.gbk")
		reference_gbk_files.append(gbk_file)
		if os.path.exists(gbk_file):
			continue
		print(f"Downloading reference genome {accession}")
		# download the reference sequence .gbk file
		with Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text") as handle, open(gbk_file, "w") as out:
			out.write(handle.read())
		# add in a 1s delay between each request for genbank
		time.sleep(1)
	
	# ----------------------------------------------------------
	# Define gene lists and alternative spellings
	# ----------------------------------------------------------
	
	# make lowercase set of canonical genes (the genes from the TSV file)
	gene_names_set = set([g.lower() for g in gene_names])
	
	# separate lists for rRNA and tRNA if present in header
	rrna_list = [g for g in gene_names if re.search(r"rrn|rRNA", g, re.IGNORECASE)]
	trna_list = [g for g in gene_names if g.lower().startswith("trn")]
	
	# create normalized tRNA variants for matching
	trna_hyphen, trna_underscore, trna_bracket, trna_no_punct = [], [], [], []
	for gene in trna_list:
		gene_l = gene.lower()
		trna_hyphen.append(gene_l)
		trna_underscore.append(gene_l.replace('-', '_'))
		trna_bracket.append(gene_l.replace('-', '(') + ")")
		trna_no_punct.append(gene_l.replace('-', '').replace('(', '').replace(')', ''))
	
	
	# ----------------------------------------------------------
	# Extract genes from reference GenBank files
	# ----------------------------------------------------------
	
	# make the fasta file directory if it doesn't already exist
	os.makedirs(args.reference_outdir, exist_ok=True)
	# track the number of times genes are written per reference genome to avoid duplicates (from the IR)
	genes_written = set()
	
	for gbk in reference_gbk_files:
		# parse .gbk file into a seq record object
		record = SeqIO.read(gbk, "genbank")
		for feature in record.features:
			# skip features that aren't genes, CDS, tRNA, or rRNA
			if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
				continue
	
			# get gene name from actual gene feature or if not from the product
			gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("product", [""]))[0].lower()
			# more complicated for rrna
			rrna_match = re.search(r'(\d+(?:\.\d+)?)s\s+ribosomal', gene_name)
			if rrna_match:
				number = rrna_match.group(1)
				gene_name_norm = f"rrn{number}"
			else:
				gene_name_norm = re.sub(r'[\s_\-()]', '', gene_name)
				gene_name_norm = gene_name_norm.replace("rna", "")
				gene_name_norm = re.sub(r'rrn(\d+(?:\.\d+)?)s$', r'rrn\1', gene_name_norm)
	
			gene_id = None
	
			# figure out which gene the gene links to in the lists of gene names generated above
			if gene_name_norm in gene_names_set:
				gene_id = gene_name
			# more complicated for trna
			elif gene_name in trna_hyphen + trna_underscore + trna_bracket + trna_no_punct:
				for idx, canonical in enumerate(trna_list):
					variants = [trna_hyphen[idx], trna_underscore[idx], trna_bracket[idx], trna_no_punct[idx]]
					if gene_name in variants:
						gene_id = canonical
						break
						
			# skip if that gene is already dealt with (IR) or not recognised                    
			if gene_id is None:
				continue
			seq_combo = (record.id, gene_id)
			if seq_combo in genes_written:
				continue
	
			# write sequence to a fasta file
			outfile = os.path.join(args.reference_outdir, f"{gene_id}.fasta")
			with open(outfile, "a") as out:
				out.write(f">{record.id} {gene_id}\n")
				out.write(str(feature.extract(record.seq)) + "\n")
	
			genes_written.add(seq_combo)


	# ----------------------------------------------------------
	# Check that all genes have at least one reference sequence
	# ----------------------------------------------------------
	
	missing_reference_genes = []

	# Compare filenames case-insensitively because GenBank gene names
	# may be normalized to lowercase during extraction
	# it's a problem that happens with \linux (not MacOS)
	existing_reference_genes = {
		os.path.splitext(f)[0].lower()
		for f in os.listdir(args.reference_outdir)
		if f.endswith(".fasta")
	}
	
	for gene in gene_names:
		if gene.lower() not in existing_reference_genes:
			missing_reference_genes.append(gene)
			
	if missing_reference_genes:
		for gen in missing_reference_genes:
			raise ValueError("ERROR: the following genes have no reference sequences:\n" +
			"\n".join(missing_reference_genes) +
			"\nAdd reference genomes containing these genes or remove them from the gene list." +
			"\nAlternatively, use fasta mode and provide your own reference gene sequences in a single fasta file.")
	

# ----------------------------------------------------------
# Download plastid FASTA for missing taxa
# ----------------------------------------------------------

# back to common script for both modes 

# make the directory for plastid fasta sequences for accessions with missing genes
os.makedirs(args.plastid_fasta_dir, exist_ok=True)

# make the file name
for i, accession in enumerate(missing_taxa_ids, start=1):
    outfile = os.path.join(args.plastid_fasta_dir, f"{accession}.fasta")
    if os.path.exists(outfile):
        continue
    # print summary to output
    print(f"[{i}/{len(missing_taxa_ids)}] Downloading plastid {accession}")
    # download fasta files for sequences with missing genes
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle, open(outfile, "w") as out:
        out.write(handle.read())
    # wait 1s before each download request from NCBI
    time.sleep(1)

# ----------------------------------------------------------
# Make BLAST databases for missing taxa
# ----------------------------------------------------------
os.makedirs(args.blast_db_dir, exist_ok=True)
for fasta_file in os.listdir(args.plastid_fasta_dir):
    if not fasta_file.endswith(".fasta"):
        continue
    accession = fasta_file.replace(".fasta", "")
    fasta_path = os.path.join(args.plastid_fasta_dir, fasta_file)
    db_prefix = os.path.join(args.blast_db_dir, accession)
    if os.path.exists(db_prefix + ".nhr"):
        continue
    print(f"Making BLAST database for {accession}")
    subprocess.run([
        "makeblastdb",
        "-in", fasta_path,
        "-out", db_prefix,
        "-dbtype", "nucl",
        "-parse_seqids"], 
    check=True)

# ----------------------------------------------------------
# Run BLAST to recover missing genes
# ----------------------------------------------------------

# make directory for blast databases to be stored in
os.makedirs(args.blast_out, exist_ok=True)

# loop through each gene and all taxa missing that gene
for gene_idx, taxa in enumerate(missing_by_gene):
    gene = gene_names[gene_idx]
    query_fasta = os.path.join(args.reference_outdir, f"{gene}.fasta")
    # skip genes with no reference sequences
    if not os.path.exists(query_fasta):
        continue

	# run BLAST using reference gene as query and test taxon's plastid genome as database
    for accession in taxa:
        genome_outdir = os.path.join(args.blast_out, accession)
        os.makedirs(genome_outdir, exist_ok=True)
        out_file = os.path.join(genome_outdir, f"{gene}.txt")
        blast_cmd = [args.blast_type,
            "-db", os.path.join(args.blast_db_dir, accession),
            "-query", query_fasta,
            "-outfmt", "6 qseqid sseqid evalue qstart qend sstart send",
            "-out", out_file]
        subprocess.run(blast_cmd, check=True)

print("BLAST complete.")