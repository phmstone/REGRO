# BLAST-based recovery of missing plastid genes.
# Need to have run a bunch of other scripts before

# Pipeline:
# 1. Read presence/absence TSV
# 2. Identify taxa missing one or more genes
# 3. Download reference GenBank files (input list needed)
# 4. Extract reference gene sequences
# 5. Download plastid FASTAs for missing taxa
# 6. Make BLAST databases
# 7. BLAST missing genes against plastomes


import argparse
import os
import time
import subprocess
import certifi
from Bio import Entrez, SeqIO

# allows python to interact with internet? Maybe fix this before I send it out
os.environ["SSL_CERT_FILE"] = certifi.where()

# ----------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Recover missing plastid genes using BLAST")

parser.add_argument(
    "--input",
    required=True,
    help="Presence/absence TSV file")

parser.add_argument(
    "--email",
    required=True,
    help="Email address for NCBI Entrez")

parser.add_argument(
    "--reference-ids",
    required=True,
    help="Text file of reference GenBank accessions (one per line)")

parser.add_argument(
    "--reference-outdir",
    default="Blast/ReferenceGeneSequences",
    help="Directory to store reference gene FASTAs")

parser.add_argument(
    "--reference-gbk-dir",
    default="Blast/ReferenceGenomes",
    help="Directory to store reference GenBank files")

parser.add_argument(
    "--plastid-fasta-dir",
    default="Blast/PlastidSequences",
    help="Directory to store plastid FASTAs")

parser.add_argument(
    "--blast-db-dir",
    default="Blast/Databases",
    help="Directory where BLAST databases will be created")

parser.add_argument(
    "--blast-out",
    default="Blast/Results",
    help="Directory for BLAST output")

args = parser.parse_args()

# ----------------------------------------------------------
# NCBI set up
# ----------------------------------------------------------

# need an account with NCBI to download things
Entrez.email = args.email

# ----------------------------------------------------------
# Read in presence/ absence CSV file
# ----------------------------------------------------------

# read in the presence/absence TSV file
with open(args.input) as fh:
    lines = [line.rstrip("\n").split("\t") for line in fh]

# separate the header of the TSV file from the presence/absence profiles
header = lines.pop(0)
gene_names = header[2:]

# make empty lists
latin_names = []
genbank_ids = []
presence_absence = []

# populate lists with relevant fields from TSV
for row in lines:
    latin_names.append(row[0])
    genbank_ids.append(row[1])
    presence_absence.append(row[2:])

print(f"{len(genbank_ids)} taxa total")

# ----------------------------------------------------------
# Identify taxa that are missing genes
# ----------------------------------------------------------

missing_taxa_ids = []
complete_taxa_ids = []

# if all genes are present in a plastid genome, then all entries will be identical 
# technically this means plastid genomes with all missing genes or all pseudogenes would go into "complete taxa" here
for i, profile in enumerate(presence_absence):
    if len(set(profile)) == 1:
        complete_taxa_ids.append(genbank_ids[i])
    else:
        missing_taxa_ids.append(genbank_ids[i])

# print summary
print(f"{len(complete_taxa_ids)} taxa with all genes present")
print(f"{len(missing_taxa_ids)} taxa missing ≥1 gene")

# ----------------------------------------------------------
# Identifying which taxa are missing which genes
# ----------------------------------------------------------

missing_by_gene = [[] for _ in range(len(gene_names))]

for taxon_idx, profile in enumerate(presence_absence):
    for gene_idx, value in enumerate(profile):
    	# "1" represents a missing gene and "2" represents a pseudogene
        if value in {"1", "2"}:
            missing_by_gene[gene_idx].append(genbank_ids[taxon_idx])

# ----------------------------------------------------------
# Read in reference accessions and make new references
# ----------------------------------------------------------

# make a list of refrence genbank IDs from the input txt file
with open(args.reference_ids) as fh:
    reference_ids = [line.strip() for line in fh if line.strip()]
    
# Add taxa from the TSV that had all genes present to reference IDs
reference_ids += complete_taxa_ids
reference_ids = list(set(reference_ids))  # remove duplicates just in case


# ----------------------------------------------------------
# Download reference genbank files
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
    with Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gbwithparts",
        retmode="text"
    ) as handle, open(gbk_file, "w") as out:
        out.write(handle.read())
    
    # add in a 1s delay between each request for genbank
    time.sleep(1)

# ----------------------------------------------------------
# Define gene lists and alternative spellings
# ----------------------------------------------------------

geneString = "ndhA	ndhB	ndhC	ndhD	ndhE	ndhF	ndhG	ndhH	ndhI	ndhJ	ndhK	ccsA	cemA	petA	petB	petD	petG	petL	petN	psaA	psaB	psaC	psaI	psaJ	psbA	psbB	psbC	psbD	psbE	psbF	psbH	psbI	psbJ	psbK	psbL	psbM	psbN	psbT	psbZ	rbcL	ycf3	ycf4	rpoA	rpoB	rpoC1	rpoC2	atpA	atpB	atpE	atpF	atpH	atpI	infA	rpl2	rpl14	rpl16	rpl20	rpl22	rpl23	rpl32	rpl33	rpl36	rps2	rps3	rps4	rps7	rps8	rps11	rps12	rps14	rps15	rps16	rps18	rps19	accD	clpP	matK	ycf1	ycf2	rrn4.5	rrn5	rrn16	rrn23	trnA-UGC	trnC-GCA	trnD-GUC	trnE-UUC	trnF-GAA	trnfM-CAU	trnG-GCC	trnG-UCC	trnH-GUG	trnI-CAU	trnI-GAU	trnK-UUU	trnL-CAA	trnL-UAA	trnL-UAG	trnM-CAU	trnN-GUU	trnP-UGG	trnQ-UUG	trnR-ACG	trnR-UCU	trnS-GCU	trnS-GGA	trnS-UGA	trnT-GGU	trnT-UGU	trnV-GAC	trnV-UAC	trnW-CCA	trnY-GUA"
gene_list = geneString.split('\t')
trna_list = gene_list[83:]
rrna_list = gene_list[79:83]
protein_genes = gene_list[:79]

trna_sesamum_list = ['tRNA-Ala (UGC)', 'tRNA-Cys (GCA)', 'tRNA-Asp (GUC)', 'tRNA-Glu (UUC)', 'tRNA-Phe (GAA)', 'tRNA-fM (CAU)', 'tRNA-Gly (GCC)', 'tRNA-Gly (UCC)', 'tRNA-His (GUG)', 'tRNA-Ile (CAU)', 'tRNA-Ile (GAU)', 'tRNA-Lys (UUU)', 'tRNA-Leu (CAA)', 'tRNA-Leu (UAA)', 'tRNA-Leu (UAG)', 'tRNA-Met (CAU)', 'tRNA-Asn (GUU)', 'tRNA-Pro (UGG)', 'tRNA-Gln (UUG)', 'tRNA-Arg (ACG)', 'tRNA-Arg (UCU)', 'tRNA-Ser (GCU)', 'tRNA-Ser (GGA)', 'tRNA-Ser (UGA)', 'tRNA-Thr (GGU)', 'tRNA-Thr (UGU)', 'tRNA-Val (GAC)', 'tRNA-Val (UAC)', 'tRNA-Trp (CCA)', 'tRNA-Tyr (GUA)']

rrna_word_list = ['23S rRNA', '4.5S rRNA', '5S rRNA', '16S rRNA']
rrn_alternate_list = ['rrn23s', 'rrn4.5s', 'rrn5s', 'rrn16']

# make lowercase sets for case-insensitive matching
protein_genes_set = set([g.lower() for g in protein_genes])
rrna_genes_set = set([g.lower() for g in rrna_list + rrna_word_list + rrn_alternate_list])
trna_genes_set = set([g.lower() for g in trna_list + trna_sesamum_list])

# create alternate tRNA spellings
trna_hyphen, trna_underscore, trna_bracket, trna_no_punct = [], [], [], []
for gene in trna_list:
    gene_l = gene.lower()
    trna_hyphen.append(gene_l)
    trna_underscore.append(gene_l.replace('-', '_'))
    trna_bracket.append(gene_l.replace('-', '(') + ")")
    trna_no_punct.append(gene_l.replace('-', '').replace('(', '').replace(')', ''))

# ----------------------------------------------------------
# Extract genes from GenBank files
# ----------------------------------------------------------

# make the fasta file directory if it doesn't already exist
os.makedirs(args.reference_outdir, exist_ok=True)

genes_written = set()

for gbk in reference_gbk_files:
	# parse .gbk file into a seq record object
    record = SeqIO.read(gbk, "genbank")

    for feature in record.features:
    	# skip features that aren't genes, CDS, tRNA, or rRNA
        if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
            continue

        gene_name = None
        # get gene name from actual gene feature or if not from the product
        if "gene" in feature.qualifiers:
            gene_name = feature.qualifiers["gene"][0].lower()
        elif "product" in feature.qualifiers:
            gene_name = feature.qualifiers["product"][0].lower()
        else:
            continue

        gene_id = None
		
		# figure out which gene the gene links to in the lists of gene names generated above
        if gene_name in protein_genes_set:
            gene_id = gene_name
        elif any(gene_name in v for v in rrna_genes_set):
            for idx, canonical in enumerate(rrna_list):
                if rrna_list[idx].lower() in gene_name or rrna_word_list[idx].lower() in gene_name or rrn_alternate_list[idx].lower() in gene_name:
                    gene_id = canonical
                    break
        elif any(gene_name in v for v in trna_genes_set):
            for idx, canonical in enumerate(trna_list):
                variants = [trna_hyphen[idx], trna_underscore[idx], trna_bracket[idx], trna_no_punct[idx], trna_sesamum_list[idx].lower()]
                if gene_name in variants:
                    gene_id = canonical
                    break
		
		# skip if that gene is already dealt with (IR) or not recognised
        if gene_id is None or gene_id in genes_written:
            continue

        # write sequence
        seq = feature.extract(record.seq)
        outfile = os.path.join(args.reference_outdir, f"{gene_id}.fasta")
        with open(outfile, "w") as out:
            out.write(f">{record.id} {gene_id}\n")
            out.write(str(seq) + "\n")

        genes_written.add(gene_id)

# make a lowercase set of genes_written for comparison
genes_written_lc = set(g.lower() for g in genes_written)

# compare to lowercase gene_list to make warning for genes not found
missing_refs = set(g.lower() for g in gene_list) - genes_written_lc
if missing_refs:
    print("WARNING: no reference sequence found for:")
    for g in sorted(missing_refs):
        print(f"  {g}")


# ----------------------------------------------------------
# Download plastid fasta files for taxa with missing genes
# ----------------------------------------------------------


# make the directory for plastid fasta sequences for accessions with missing genes
os.makedirs(args.plastid_fasta_dir, exist_ok=True)


for i, accession in enumerate(missing_taxa_ids, start=1):
    outfile = os.path.join(args.plastid_fasta_dir, f"{accession}.fasta")

    if os.path.exists(outfile):
        continue
	# print summary to output
    print(f"[{i}/{len(missing_taxa_ids)}] Downloading plastid {accession}")
	# download fasta files for sequences with missing genes
    with Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    ) as handle, open(outfile, "w") as out:
        out.write(handle.read())
    # wait 1s before each download request from NCBI
    time.sleep(1)



# ----------------------------------------------------------
# Make blast databases for taxa with missing genes
# ----------------------------------------------------------


# make directory for blast databases to be stored in
os.makedirs(args.blast_db_dir, exist_ok=True)

# for fasta file in the fasta file directory make a blast database
for fasta_file in os.listdir(args.plastid_fasta_dir):
    if not fasta_file.endswith(".fasta"):
        continue

	# establish file paths and file names
    accession = fasta_file.replace(".fasta", "")
    fasta_path = os.path.join(args.plastid_fasta_dir, fasta_file)
    db_prefix = os.path.join(args.blast_db_dir, accession)

    # Check whether database already exists
    if os.path.exists(db_prefix + ".nhr"):
        continue

    print(f"Making BLAST database for {accession}")
	
	# make the blast database
    makeblastdb_cmd = [
        "makeblastdb",
        "-in", fasta_path,
        "-out", db_prefix,
        "-dbtype", "nucl",
        "-parse_seqids"
    ]

    subprocess.run(makeblastdb_cmd, check=True)




# ----------------------------------------------------------
# run blast with reference gene sequences on missing data
# ----------------------------------------------------------


# make the directory for blast results to be stored in
os.makedirs(args.blast_out, exist_ok=True)

for gene_idx, taxa in enumerate(missing_by_gene):
    gene = gene_names[gene_idx]
    query_fasta = os.path.join(args.reference_outdir, f"{gene}.fasta")

    if not os.path.exists(query_fasta):
        continue

    for accession in taxa:
        out_file = os.path.join(args.blast_out, f"{gene}-{accession}.txt")

        # use blastn for nucleotide blast
        blast_cmd = [
            "blastn",
            "-db", os.path.join(args.blast_db_dir, accession),
            "-query", query_fasta,
            "-outfmt", "7 qaccession saccession evalue qstart qend sstart send",
            "-out", out_file
        ]

        subprocess.run(blast_cmd, check=True)

#print("Pipeline complete.")
