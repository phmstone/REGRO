####################################################################################################
# This script:
# 1. Takes a .txt file of genbank IDs as input (one per line)
# 2. Downloads a .gbk file containing all .gbk records for the IDs specified in the .txt file
####################################################################################################


import os # for interacting with operating system of computer
import time # allows 1s delay between downloading genbank files
import certifi # security certificate problems need this 
import argparse # allows command line inputs to be included
from Bio import Entrez # for interacting with genbank

# allows python to interact with internet? Maybe fix this before I send it out
os.environ["SSL_CERT_FILE"] = certifi.where()

# -----------------------------
# Command line arguments
# -----------------------------

parser = argparse.ArgumentParser(
    description="Download GenBank records from a list of accession IDs")

# user has to input their own email address in order to download from genbank
parser.add_argument(
    "--email",
    required=True,
    help="Email address (required by NCBI Entrez)")

# user has to give input file
parser.add_argument(
    "--input",
    required=True,
    help="Text file containing GenBank accession IDs (one per line)")

# user has to give output file name
parser.add_argument(
    "--output",
    required=True,
    help="Output GenBank file (.gbk)")

# user can give their own time delay between requesting .gbk downloads from genbank if they want
# 1s is the default 
parser.add_argument(
    "--delay",
    type=float,
    default=1.0,
    help="Delay in seconds between NCBI requests (default: 1)")

# enables above inputs to be parsed 
args = parser.parse_args()

# -----------------------------
# Load and normalize GenBank IDs
# -----------------------------

# always identify yourself to NCBI
Entrez.email = args.email

# some accessions are entered multiple times with IDs that end in decimals after the first upload
# here I take the first upload only by ignoring uploads with decimal places after them
with open(args.input) as f:
    genbank_ids = {line.strip().split(".")[0] for line in f if line.strip()}

print(f"Total unique IDs accessions found: {len(genbank_ids)}")


# -----------------------------
# Download records
# -----------------------------
print("1 s delay between requests so genbank does not crash")

# putting the output file in the output directory made in the above chunk (if does not already exist)
output_file = os.path.join(args.output)

# write all the .gbk files from each accession into one big .gbk file which will be the output file
with open(output_file, "w") as out_handle:
	# loop over all the accessions found in the previous chunk
    for i, accession in enumerate(genbank_ids, start=1):
        try:
        	# download the file from genbank (only works if ID is valid because of above "try")
            with Entrez.efetch(
                db="nucleotide", # using the genbank nucleotide database
                id=accession, # accession number
                rettype="gbwithparts", # the type of .gbk file being downloaded
                retmode="text" # text format
            ) as net_handle:
                out_handle.write(net_handle.read())

            print(f"{i}: downloaded {accession}") # prints the number accession it is and accession ID
            # wait the specified number (or 1) second before downloading another file
            time.sleep(args.delay) 
		
		# if the genbank ID is invalid then this fail message prints
        except Exception as e:
            print(f"Failed for {accession}: {e}")

print("Saved")
