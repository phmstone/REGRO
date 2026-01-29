########################################################################################################################################################
# This script:
# 1. Reads in a multi-sequence .gbk file
# 2. Parses annotations in the .gbk file by taxon
# 3. Creates one output multifasta file per gene (specified in gene_list) containing sequences for all taxa where that gene was annotated as present
# 4. The multifasta files are not aligned 
########################################################################################################################################################


import argparse #command line function
import os #enables making new directories and moving around the directory structure
from Bio import SeqIO #for parsing genbank files


# ----------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------

# set up argument parsing for the command line
parser = argparse.ArgumentParser(
    description="Extract chloroplast genes from a GenBank file into a multifasta file per gene")

# input GenBank file (required)
parser.add_argument(
 	"--input",
    required=True,
    help="Input GenBank file")

# output directory for fasta files (optional)
# if no output directory name is given it will be called PresentGeneMultiFastas
parser.add_argument(
    "--outdir",
    default="PresentGeneAlignments",
    help="Output directory (default: PresentGeneMultiFastas)")

# parse command-line arguments
args = parser.parse_args()

# create the output directory if it does not already exist
os.makedirs(args.outdir, exist_ok=True)



# ----------------------------------------------------------
# Define reference gene lists
# ----------------------------------------------------------

# genes come in lots of different spelling variations so need to catch them all 
# input the list of genes as a string and break it up into a list by the tabs in between 
geneString = "ndhA	ndhB	ndhC	ndhD	ndhE	ndhF	ndhG	ndhH	ndhI	ndhJ	ndhK	ccsA	cemA	petA	petB	petD	petG	petL	petN	psaA	psaB	psaC	psaI	psaJ	psbA	psbB	psbC	psbD	psbE	psbF	psbH	psbI	psbJ	psbK	psbL	psbM	psbN	psbT	psbZ	rbcL	ycf3	ycf4	rpoA	rpoB	rpoC1	rpoC2	atpA	atpB	atpE	atpF	atpH	atpI	infA	rpl2	rpl14	rpl16	rpl20	rpl22	rpl23	rpl32	rpl33	rpl36	rps2	rps3	rps4	rps7	rps8	rps11	rps12	rps14	rps15	rps16	rps18	rps19	accD	clpP	matK	ycf1	ycf2	rrn4.5	rrn5	rrn16	rrn23	trnA-UGC	trnC-GCA	trnD-GUC	trnE-UUC	trnF-GAA	trnfM-CAU	trnG-GCC	trnG-UCC	trnH-GUG	trnI-CAU	trnI-GAU	trnK-UUU	trnL-CAA	trnL-UAA	trnL-UAG	trnM-CAU	trnN-GUU	trnP-UGG	trnQ-UUG	trnR-ACG	trnR-UCU	trnS-GCU	trnS-GGA	trnS-UGA	trnT-GGU	trnT-UGU	trnV-GAC	trnV-UAC	trnW-CCA	trnY-GUA"
gene_list = geneString.split('\t')
# extract the tRNA list
trna_list = gene_list[83:]
# extract the rrna list
rrna_list = gene_list[79:83]
# extract the protein coding gene list
geneList = gene_list[:79]

# have only seen this style in one genbank file (Sesamum) but it's strange and must be spelled out in full
trna_sesamum_list = ['tRNA-Ala (UGC)', 'tRNA-Cys (GCA)', 'tRNA-Asp (GUC)', 'tRNA-Glu (UUC)', 'tRNA-Phe (GAA)', 'tRNA-fM (CAU)', 'tRNA-Gly (GCC)', 'tRNA-Gly (UCC)', 'tRNA-His (GUG)', 'tRNA-Ile (CAU)', 'tRNA-Ile (GAU)', 'tRNA-Lys (UUU)', 'tRNA-Leu (CAA)', 'tRNA-Leu (UAA)', 'tRNA-Leu (UAG)', 'tRNA-Met (CAU)', 'tRNA-Asn (GUU)', 'tRNA-Pro (UGG)', 'tRNA-Gln (UUG)', 'tRNA-Arg (ACG)', 'tRNA-Arg (UCU)', 'tRNA-Ser (GCU)', 'tRNA-Ser (GGA)', 'tRNA-Ser (UGA)', 'tRNA-Thr (GGU)', 'tRNA-Thr (UGU)', 'tRNA-Val (GAC)', 'tRNA-Val (UAC)', 'tRNA-Trp (CCA)', 'tRNA-Tyr (GUA)']

# write out alternative rrna gene name formats
rrna_word_list = ['23S rRNA', '4.5S rRNA', '5S rRNA', '16S rRNA']
rrna_alternate_list = ['rrn16', 'rrn4.5s', 'rrn5s', 'rrn23s']

# make everything lowercase for case insensitive matching
gene_list = [g.lower() for g in gene_list]
rrna_list = [g.lower() for g in rrna_list]
trna_list = [g.lower() for g in trna_list]
trna_sesamum_list = [g.lower() for g in trna_sesamum_list]


# ----------------------------------------------------------
# Generate alternative tRNA spellings
# ----------------------------------------------------------

# make empty lists for alternate spellings of tRNA genes
trna_hyphen = []
trna_underscore = []
trna_bracket = []
trna_no_punct = []

# make the alternate trna spellings
for gene in trna_list:
    trna_hyphen.append(gene)
    trna_underscore.append(gene.replace('-', '_'))
    trna_bracket.append(gene.replace('-', '(') + ")")
    trna_no_punct.append(gene.replace('-', '').replace('(', '').replace(')', ''))

# make everything lowercase for case insensitive matching
trna_hyphen = [g.lower() for g in trna_hyphen]
trna_underscore = [g.lower() for g in trna_underscore]
trna_bracket = [g.lower() for g in trna_bracket]
trna_no_punct = [g.lower() for g in trna_no_punct]

# ----------------------------------------------------------
# Read the genbank file
# ----------------------------------------------------------

# parse the GenBank file into a dictionary of SeqRecord objects
# each key in the dictionary is one genbank accession
records = SeqIO.to_dict(SeqIO.parse(args.input, "genbank"))


# ----------------------------------------------------------
# Loop to extract genes into multifasta files
# ----------------------------------------------------------

# loop over each accession in the GenBank file
for recordID in records:
    record = records[recordID]

    # loop over all annotated features in the record
    for feature in record.features:

        # only deal with features that are annotated as genes
        if feature.type != "gene":
            continue

        # skip features without a gene qualifier
        if "gene" not in feature.qualifiers:
            continue

        # a feature can have multiple gene names; loop over them too
        for gene_name in feature.qualifiers["gene"]:
        
        	# convert gene name to lowercase for case-insensitive matching
            gene_name_lc = gene_name.lower()
            
            # protein coding gene names
            if gene_name_lc in gene_list:
            	gene = gene_name_lc
                
            # rrna gene names
            elif (
                gene_name_lc in rrna_list or
                gene_name_lc in rrna_word_list or
                gene_name_lc in rrna_alternate_list
            ):
                gene = gene_name_lc

			# trna gene names
            elif (
                gene_name_lc in trna_hyphen or
                gene_name_lc in trna_underscore or
                gene_name_lc in trna_bracket or
                gene_name_lc in trna_no_punct or
                gene_name_lc in trna_sesamum_list
            ):
                gene = gene_name_lc

            # if the gene name does not match any known spelling, skip it
            else:
                continue

            # extract the DNA sequence for the gene.
            sequence = ""
            for part in feature.location.parts:
                sequence += str(part.extract(record.seq))

            # get the start and end points for reporting
            start = int(feature.location.start)
            end = int(feature.location.end)

            # make the output fasta filename (one file per gene)
            out_file = os.path.join(
                args.outdir,
                f"{gene}_alignment_unaligned.fasta"
            )

            # append the sequence to the gene's fasta file
            with open(out_file, "a") as fh:
                fh.write(f">{recordID} | {gene} {start}-{end}\n")
                fh.write(sequence + "\n")

            # print progress to the terminal
            # maybe cut because this is a bit much
            print(f"{recordID} | {gene} {start}-{end}")
