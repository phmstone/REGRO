################################################################################################################
# Combined Script:
# 1. Reads in a multi-sequence .gbk file
# 2. Parses annotations in the .gbk file by taxon
# 3. Creates a TSV stating whether a gene is present, missing, or pseudogenized
# 4. Optionally makes a nexus-style text block for character state mapping
# 5. Outputs multifastas by gene for taxa that have genes annotated as being present
################################################################################################################

# --------------------------------------------------------------------------------------------------------------
# Import packages
# --------------------------------------------------------------------------------------------------------------

import argparse
import re
import os
import sys 
from Bio import SeqIO

# --------------------------------------------------------------------------------------------------------------
# Command line arguments
# --------------------------------------------------------------------------------------------------------------

# This allows the script to be run from the command line
parser = argparse.ArgumentParser(description="Extract chloroplast genes: presence/absence TSV, optional Nexus, and multifasta per gene")

# Required arguments
# input GenBank file
parser.add_argument("--input", required=True, help="Input GenBank (.gbk) file")
# TSV output
parser.add_argument("--tsv", required=True, help="Output TSV file for gene presence/absence")


# Optional arguments
# Optional gene list file
parser.add_argument("--gene_file", help="Optional file containing one gene per line (overrides default list)")
# nexus block file 
parser.add_argument("--nexus", help="Optional Nexus-style text block output")
# multifasta output directories
parser.add_argument("--outdir", default="PresentGeneMultiFastas",
                    help="Directory to store full gene sequences (default: PresentGeneMultiFastas)")
parser.add_argument("--coding_outdir", default="PresentCodingSeqMultiFastas",
                    help="Directory to store coding sequences (default: PresentCodingSeqMultiFastas)")
parser.add_argument("--pseudo_outdir", default="PseudogeneMultiFastas",
                    help="Directory to store pseudogene sequences (default: PseudogeneMultiFastas)")
# gene alias file
parser.add_argument("--alias_file", help="Optional gene alias file: canonical_name <tab> synonym")

# Parse the arguments
args = parser.parse_args()

# Create output directories if they don't exist
os.makedirs(args.outdir, exist_ok=True)
if args.coding_outdir:
    os.makedirs(args.coding_outdir, exist_ok=True)    
if args.pseudo_outdir:
    os.makedirs(args.pseudo_outdir, exist_ok=True)

# --------------------------------------------------------------------------------------------------------------
# Define gene list (default or user-provided)
# --------------------------------------------------------------------------------------------------------------

# Add in the default genes tested (canonical angiosperm plastid genes)
default_gene_string = (
    "ndhA\tndhB\tndhC\tndhD\tndhE\tndhF\tndhG\tndhH\tndhI\tndhJ\tndhK\t"
    "ccsA\tcemA\tpetA\tpetB\tpetD\tpetG\tpetL\tpetN\t"
    "psaA\tpsaB\tpsaC\tpsaI\tpsaJ\t"
    "psbA\tpsbB\tpsbC\tpsbD\tpsbE\tpsbF\tpsbH\tpsbI\tpsbJ\tpsbK\tpsbL\tpsbM\tpsbN\tpsbT\tpsbZ\t"
    "rbcL\tycf3\tycf4\t"
    "rpoA\trpoB\trpoC1\trpoC2\t"
    "atpA\tatpB\tatpE\tatpF\tatpH\tatpI\t"
    "infA\t"
    "rpl2\trpl14\trpl16\trpl20\trpl22\trpl23\trpl32\trpl33\trpl36\t"
    "rps2\trps3\trps4\trps7\trps8\trps11\trps12\trps14\trps15\trps16\trps18\trps19\t"
    "accD\tclpP\tmatK\tycf1\tycf2\t"
    "rrn4.5\trrn5\trrn16\trrn23\t"
    "trnA-UGC\ttrnC-GCA\ttrnD-GUC\ttrnE-UUC\ttrnF-GAA\ttrnfM-CAU\t"
    "trnG-GCC\ttrnG-UCC\ttrnH-GUG\ttrnI-CAU\ttrnI-GAU\ttrnK-UUU\t"
    "trnL-CAA\ttrnL-UAA\ttrnL-UAG\ttrnM-CAU\ttrnN-GUU\ttrnP-UGG\t"
    "trnQ-UUG\ttrnR-ACG\ttrnR-UCU\ttrnS-GCU\ttrnS-GGA\ttrnS-UGA\t"
    "trnT-GGU\ttrnT-UGU\ttrnV-GAC\ttrnV-UAC\ttrnW-CCA\ttrnY-GUA")

# Load default or user-supplied gene list if it's been supplied
if args.gene_file:
    with open(args.gene_file) as gf:
        gene_list = [line.strip() for line in gf if line.strip()]
    print(f"Using custom gene list from {args.gene_file}")
else:
    gene_list = default_gene_string.split("\t")
    print("Using default angiosperm plastid gene list")
    
# empty output fasta files if they already exist by making new ones 
# stops endless appending of sequences
for gene in gene_list:
    out_file = os.path.join(args.outdir, f"{gene}_alignment_unaligned.fasta")
    open(out_file, "w").close()  # empties file

    if args.coding_outdir:
        coding_file = os.path.join(args.coding_outdir, f"{gene}_coding_unaligned.fasta")
        open(coding_file, "w").close()  # empties coding file

# normalise gene names for matching
normalised_targets = {}
for g in gene_list:
    norm = g.lower()
    norm = re.sub(r'[\s_\-()]', '', norm)  # remove punctuation, underscores, hyphens, parentheses
    norm = norm.replace("rna", "")          # remove "rna" from names
    normalised_targets[norm] = g            # mapping normalised back to the canonical names

# deal with the alias file if it was included
if args.alias_file:
    with open(args.alias_file) as af:
        for line in af:
            # remove trailing whitespace 
            line = line.strip()
            # ignore the comments
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # if the line is not separated by one tab then just skip and ignore it
            if len(parts) != 2:
                print(f"Skipping malformed alias line: {line}")
                continue
            # in the file the canonical name is separated by the synonym by one line
            canonical, synonym = parts
            # normalise the canonical names
            norm_canon = canonical.lower()
            norm_canon = re.sub(r'[\s_\-()]', '', norm_canon)
            norm_canon = norm_canon.replace("rna", "")
            # throws an error if the normalised canonical name is not in the normalised gene list
            if norm_canon not in normalised_targets:
                print(f"Warning: canonical '{canonical}' not in gene list")
                continue
            canonical = normalised_targets[norm_canon] 
            # normalise synonym exactly like GenBank names (is how they will be normalised later on when reading .gbk file)
            norm_syn = synonym.lower()
            norm_syn = re.sub(r'[\s_\-()]', '', norm_syn)
            norm_syn = norm_syn.replace("rna", "")
            # what happens if a gene name is already included in the gene list or included twice in the alias list?
            if norm_syn in normalised_targets:
                existing = normalised_targets[norm_syn]
                # error messge if normalised synonyms from alias file map to two canonical gene names
                if existing != canonical:
                    sys.exit(f"\nERROR: Alias conflict detected.\n"
                        f"synonym '{synonym}' maps to both "
                        f"'{existing}' and '{canonical}'. \n"
                        f"Check and fix the alias file.")
                else:
                    # synonym given normalises to the gene name, then just carry on 
                    print(f"Warning: synonym '{synonym}' overwrites existing mapping, check the alias file for duplicate synonyms")
                    print(f"This usually means that an alternative spelling/format was simply normalised to a gene name already included in the gene list")
                    continue
            # Add synonym into main lookup dictionary
            normalised_targets[norm_syn] = canonical
    print(f"Alias file loaded.")

# maybe cut this
# the number in {len(normalised_targets)} is higher if an alias file is used because it includes synonyms
print(f"Tracking {len(gene_list)} genes.")

# --------------------------------------------------------------------------------------------------------------
# Read GenBank file
# --------------------------------------------------------------------------------------------------------------

# make a dictionary out of the genbank files
# the genbank IDs are the keys and the other features of the genbank files as values

try:
    records = SeqIO.to_dict(SeqIO.parse(args.input, "genbank"))
except Exception as e:
    print(f"Error reading GenBank file: {e}")
    records = {}

print(f"Number of sequences read: {len(records)}")

# --------------------------------------------------------------------------------------------------------------
# Keep sequences without gene annotations, but print a warning
# --------------------------------------------------------------------------------------------------------------

# if a .gbk file has no genes present as features then there is a problem

# make a list of sequences to remove
no_gene_list = []

# figure out if a gene is present
for key, seq_record in records.items():
    has_gene = any("gene" in f.qualifiers for f in seq_record.features)
    if not has_gene:
        no_gene_list.append(key)

# print the genbank IDs that have no gene annotations if there are any
if len(no_gene_list) > 0:
    print(f"There are {len(no_gene_list)} sequences with no gene annotations: {no_gene_list}")
    print("Check the GenBank file to see if an alias file is needed to accmmodate alternative gene names for these taxa.")
# print the number of sequences moving forward
print(f"Remaining sequences: {len(records)}")


# --------------------------------------------------------------------------------------------------------------
# Process each GenBank record
# --------------------------------------------------------------------------------------------------------------

# make empty lists to store values in
taxa_names = []
genbank_ids = []
all_present_genes = []
all_pseudogenes = []
multi_accession_log = []

# loop over taxa in the dictionary 
for recordID, record in records.items():
    skip_record = False
    # extract taxon name from genabnk record
    taxon_name = "_".join(record.description.split()[:2])


    # establish empty lists to be filled in
    present_genes = set()
    pseudogenes = set()

    # Dictionaries to store sequences for multifastas
    full_gene_candidates = {}
    coding_gene_candidates = {}
    pseudogene_candidates = {}

    # Process features
    for feature in record.features:
        # only consider sequences that could be relevant
        if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
            continue
        gene_name = feature.qualifiers.get("gene", [""])[0]
        product_name = feature.qualifiers.get("product", [""])[0]
        combined_name = gene_name if gene_name else product_name
        # if no name match then move on
        if not combined_name:
            continue

        # normalise the spelling names in the same way as above to find a match
        norm_name = combined_name.lower()
        # special case for some rrna spellings where whole word is used 
        rrna_match = re.search(r'(\d+(?:\.\d+)?)s\s+ribosomal', norm_name)
        if rrna_match:
            norm_name = f"rrn{rrna_match.group(1)}"
        # standard normalisation for everything else
        else:
            norm_name = re.sub(r'[\s_\-()]', '', norm_name)
            # removes s from the end of e.g. rrn16s
            norm_name = re.sub(r'rrn(\d+(?:\.\d+)?)s$', r'rrn\1', norm_name)
            norm_name = norm_name.replace("rna", "")
        # if no name match then move on
        if norm_name not in normalised_targets:
            continue
        canonical_name = normalised_targets[norm_name]

        # if there's pseudo in the gene information, classify it as a pseudogene
        is_pseudo = "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers
        if is_pseudo:
            pseudogenes.add(canonical_name)
        else:
            present_genes.add(canonical_name)

        # -----------------------------------------------------------------------------------------------------------
        # Storing sequences for profiles and writing out to multifastas
        # -----------------------------------------------------------------------------------------------------------

        # Pre-define start and end safely for all location types
        start, end = sorted([int(feature.location.start), int(feature.location.end)])

        # Make a multi-accession flag (reset for each new feature)
        is_multi_accession = False

        try:
            extracted_seq = str(feature.extract(record.seq))
            
        except ValueError as error:
            if "another sequence" in str(error) or "references mandatory" in str(error):
                # Attempt to manually resolve multi-accession / trans-spliced exons
                extracted_pieces = []
                try:
                    parts = getattr(feature.location, "parts", [feature.location])
                    sub_seq_failed = False
                    
                    for part in parts:
                        part_seq = ""
                        target_seq_obj = None
                        
                        # Determine if this part references another accession
                        ref_acc = None
                        if hasattr(part, "ref") and part.ref:
                            ref_acc = part.ref
                        else:
                            # Fallback: check string representation of the part for accession pattern
                            part_str = str(part)
                            match = re.search(r'([A-Z]+\d+\.\d+):', part_str)
                            if match:
                                ref_acc = match.group(1)
                                
                        if ref_acc and ref_acc != record.id:
                            if ref_acc in records:
                                target_seq_obj = records[ref_acc].seq
                            else:
                                sub_seq_failed = True
                                break
                        else:
                            target_seq_obj = record.seq
                            
                        # Extract coordinates safely using integer start/end of the part
                        p_start = int(part.start)
                        p_end = int(part.end)
                        
                        sub_seq = str(target_seq_obj[p_start:p_end])
                        
                        # Handle strand orientation
                        if part.strand == -1:
                            from Bio.Seq import Seq
                            sub_seq = str(Seq(sub_seq).reverse_complement())
                            
                        extracted_pieces.append(sub_seq)
                        
                    # If the other half of a multi accession trans-spliced sequence isn't downloaded
                    if sub_seq_failed:
                        raise ValueError("Remote accession not found in input file")
                        
                    extracted_seq = "".join(extracted_pieces)
                    start, end = sorted([int(feature.location.start), int(feature.location.end)])
                    
                    # Raise flag for successfully joining exons from different accessions
                    is_multi_accession = True
                    
                    # Get information on the multi accession joins for a log file
                    # Empty list to collect all individual accession identifiers referenced by  split exons/parts of multi-accession feature.
                    used_accessions = []
                    
                    for p in parts:
                        # Check if the location part stores a cross-reference accession
                        part_ref = getattr(p, "ref", None)
                        if part_ref:
                            # If it does, append it
                            used_accessions.append(part_ref)
                        else:
                            # If can't find a cross-ref, then use a regular expression to search for a standard GenBank accession pattern.
                            match_acc = re.search(r'([A-Z]+\d+(?:\.\d+)?)', str(p))
                            # make sure match is genuine and that it is not the same as the current GenBank ID
                            if match_acc and match_acc.group(1) != recordID:
                                used_accessions.append(match_acc.group(1))
                            else:
                                # Default back to the primary record ID if no external accession identifier is detected
                                used_accessions.append(recordID)
                    
                    # remove duplicate accessions and sort the remaining ones
                    used_accessions = sorted(list(set(used_accessions)))
                    
                    # information for the sequence names and multi accession log file
                    multi_accession_log.append({
                        "species": taxon_name,
                        "genbankID": recordID,
                        "gene": canonical_name,
                        "accessions_joined": ", ".join(used_accessions)
                    })
                    
                except Exception as inner_err:
                    print(f"\nWARNING: Could not resolve cross-accession feature for {recordID} ({canonical_name}): {inner_err}")
                    skip_record = True
                    break
            else:
                raise

         # Pseudogenes
        if is_pseudo:
            if canonical_name not in pseudogene_candidates or len(extracted_seq) > len(pseudogene_candidates[canonical_name][0]):
                pseudogene_candidates[canonical_name] = (extracted_seq, start, end, is_multi_accession)
            continue

        # Full gene sequence (only type "gene")
        if feature.type == "gene":
            # store the gene if first time seen, if the name has been seen before then keep the longer copy
            if canonical_name not in full_gene_candidates or len(extracted_seq) > len(full_gene_candidates[canonical_name][0]):
                full_gene_candidates[canonical_name] = (extracted_seq, start, end, is_multi_accession)

        # Coding sequence (CDS/tRNA/rRNA)
        if feature.type in {"CDS", "tRNA", "rRNA"}:
            # store both sequence and feature for exon extraction
            if canonical_name not in coding_gene_candidates or len(extracted_seq) > len(coding_gene_candidates[canonical_name][0]):
                coding_gene_candidates[canonical_name] = (extracted_seq, feature, is_multi_accession)

    # skip the sequences that refer to other sequences
    if skip_record:
        continue
    
    # append names once they have been deemed safe (no exons in other accessions)
    taxa_names.append(taxon_name)
    genbank_ids.append(record.name)
    
    # if a gene is included in pseudogenes and present genes, remove it from pseudogenes
    # potentially could be an issue if partial copies are present at the end of one IR
    pseudogenes = pseudogenes - present_genes
    all_present_genes.append(present_genes)
    all_pseudogenes.append(pseudogenes)

    # ------------------------------------------------------------------------------------------------------------------------
    # Write multifasta files for full gene sequences, CDS, and pseudogenes
    # ------------------------------------------------------------------------------------------------------------------------
    
    # make a list for genes that have been written out to file
    written_full = set()
    # loop through the dictionary created earlier
    for cname, (seq, start, end, multi_flag) in full_gene_candidates.items():
        # unique ID for each gene
        key = (cname, start, end)
        # skip exact duplicates (sometimes tRNA and rRNA are written as two features)
        if key in written_full:
            continue  
        written_full.add(key)
        # write out the gene to file
        out_file = os.path.join(args.outdir, f"{cname}_alignment_unaligned.fasta")
        with open(out_file, "a") as fh:
            # include the genbank ID, species name, and gene coordinates
            if multi_flag:
                fh.write(f">{recordID} : {cname} {start}-{end} [MULTI-ACCESSION]\n")
            else:
                fh.write(f">{recordID} : {cname} {start}-{end}\n")
            fh.write(seq + "\n")


    # make a list for CDS that have been written out to file
    written_coding = set()
    # CDS is optional but full multifastas (above) is innate
    if args.coding_outdir:
        # loop through the dictionary created earlier
        for cname, (seq, feature_obj, multi_flag) in coding_gene_candidates.items():
            start, end = sorted([int(feature_obj.location.start), int(feature_obj.location.end)])
            # unique ID for each gene
            key = (cname, start, end)
            # skip exact duplicates (sometimes tRNA and rRNA are written as two features)
            if key in written_coding:
                continue
            written_coding.add(key)
            # get exon boundaries
            exon_ranges = []
            if hasattr(feature_obj.location, "parts") and feature_obj.location.parts:
                for part in feature_obj.location.parts:
                    exon_ranges.append(f"{int(part.start)+1}-{int(part.end)}")
            else:
                exon_ranges.append(f"{int(feature_obj.location.start)+1}-{int(feature_obj.location.end)}")
            exon_str = ";".join(exon_ranges)
            # write out the file
            out_file = os.path.join(args.coding_outdir, f"{cname}_coding_unaligned.fasta")
            with open(out_file, "a") as fh:
                if multi_flag:
                    fh.write(f">{recordID} : {cname} exons={exon_str} [MULTI-ACCESSION]\n")
                else:
                    fh.write(f">{recordID} : {cname} exons={exon_str}\n")
                fh.write(seq + "\n")

    # writing out pseudogenes is also optional
    if args.pseudo_outdir:   
        written_pseudo = set()
        for cname, (seq, start, end) in pseudogene_candidates.items():
            # skip if a functional copy of the gene exists
            if cname in present_genes:
                continue
            # skip exact duplicates
            key = (cname, start, end)
            if key in written_pseudo:
                continue
            written_pseudo.add(key)
            # write out the file
            out_file = os.path.join(args.pseudo_outdir, f"{cname}_pseudogene_unaligned.fasta")
            with open(out_file, "a") as fh:
                fh.write(f">pseudogene_{recordID} : {cname} {start}-{end}\n")
                fh.write(seq + "\n")
# --------------------------------------------------------------------------------------------------------------
# Generate gene presence/absence profile for TSV/Nexus
# --------------------------------------------------------------------------------------------------------------

gene_profiles = []

for i in range(len(taxa_names)):
    profile = []
    for gene in gene_list:
        if gene in all_pseudogenes[i]:
            profile.append("2")
        elif gene in all_present_genes[i]:
            profile.append("0")
        else:
            profile.append("1")
    gene_profiles.append(profile)

print("Gene presence/absence profiles generated.")

# --------------------------------------------------------------------------------------------------------------
# Optional nexus-style output
# --------------------------------------------------------------------------------------------------------------

# the nexus file style output is optional and will only be generated if specified on the command line
# useful if using the same species (with the same names) to make a phylogeny
# can use this block for character state mapping e.g. in Mesquite

if args.nexus:
    longest_taxon = max(taxa_names, key=len)
    with open(args.nexus, "w") as f:
        for i, profile in enumerate(gene_profiles):
            f.write(taxa_names[i].ljust(len(longest_taxon) + 5) + "".join(profile) + "\n")
    print(f"Nexus-style text file written: {args.nexus}")

# --------------------------------------------------------------------------------------------------------------
# Write TSV output
# --------------------------------------------------------------------------------------------------------------

if args.tsv:
    with open(args.tsv, "w") as f:
        f.write("speciesName\tgenbankID\t" + "\t".join(gene_list) + "\n")
        for i in range(len(taxa_names)):
            f.write(f"{taxa_names[i]}\t{genbank_ids[i]}\t" + "\t".join(gene_profiles[i]) + "\n")
    print(f"TSV file written: {args.tsv}")

# --------------------------------------------------------------------------------------------------------------
# Write multi-accession join log
# --------------------------------------------------------------------------------------------------------------

if multi_accession_log:
    join_log_file = "multiAccessionJoins.tsv"
    import csv
    
    # Filter out entries where the genbankID is the only thing listed in accessions_joined and only keep unique rows.
    seen = set()
    unique_log = []
    for row in multi_accession_log:
        # Split the comma-separated accessions back into an individual list of strings
        joined_items = [acc.strip() for acc in row["accessions_joined"].split(",")]
        # See what was joined besides the main ID
        other_accessions = [acc for acc in joined_items if acc != row["genbankID"]]
        
        # Only log if there are actual external accessions linked
        if other_accessions:
            row["accessions_joined"] = ", ".join(other_accessions)
            row_tuple = tuple(sorted(row.items()))
            # Only log each unique combination once
            if row_tuple not in seen:
                seen.add(row_tuple)
                unique_log.append(row)
    # write out the file
    if unique_log:
        with open(join_log_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["species", "genbankID", "gene", "accessions_joined"], delimiter="\t")
            writer.writeheader()
            writer.writerows(unique_log)
        print(f"Multi-accession join report written: {join_log_file}")
    else:
        print("No external multi-accession joins found to log.")

print("Processing complete.")