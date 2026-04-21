#!/usr/bin/env bash

# a bash script to run the whole pipeline in one go on a local computer
# assumes that the pipeline and dependencies are already installed

# Configuration

# path to scripts
scriptDir="../../Scripts"
# email for your NCBI account
email="your.NCBI.account@email.address"
# where all outputs from this script will be placed
projectDir="ericalesExample"

# Inputs
inputGenBankIDs="ericaceaeGenBankIDs.txt"
referencesID="referenceIDs.txt"
aliasFile="gene_alias.txt"

# Outputs
multiGenBankFile="$projectDir/Ericales.gbk"
presenceAbsenceTSV="$projectDir/EricalesPresenceAbsence.tsv"
presentGeneSequences="$projectDir/presentGeneSequences"
presentGeneCodingSequences="$projectDir/presentGeneCodingSequences"
originalHeatmap="$projectDir/EricalesHeatMap.png"
blastOutDir="$projectDir/Blast"
blastResults="$projectDir/Blast/Results"
plastidFastas="$projectDir/Blast/PlastidSequences"
referenceGenes="$projectDir/Blast/ReferenceGeneSequences"
unalignedMultiFastas="$projectDir/unalignedMultiFastas"
alignedMultiFastas="$projectDir/alignedMultiFastas"
updatedPresenceAbsenceTSV="$projectDir/EricalesPresenceAbsence-updated.tsv"
changeLog="$projectDir/changeLog.tsv"
updatedHeatmap="$projectDir/EricalesHeatMap-updated.png"

# Create output directory
mkdir -p "$projectDir"

# Downloading files from GenBank
python3 "$scriptDir/fetchGenBank.py" \
  --email "$email" \
  --input "$inputGenBankIDs" \
  --output "$multiGenBankFile"

# Reading annotations to determine which genes are present or absent
python3 "$scriptDir/presenceAbsence.py" \
  --input "$multiGenBankFile" \
  --tsv "$presenceAbsenceTSV" \
  --alias_file "$aliasFile" \
  --outdir "$presentGeneSequences" \
  --coding_outdir "$presentGeneCodingSequences"

# Generating a heat map based on GenBank annotations
python3 "$scriptDir/heatMapPlot.py" \
  "$presenceAbsenceTSV" \
  -o "$originalHeatmap"

# Using BLAST to find genes/gene fragments that may not have been annotated
python3 "$scriptDir/blastPresenceAbsence.py" \
  --input "$presenceAbsenceTSV" \
  --email "$email" \
  --reference-ids "$referencesID" \
  --blast-type blastn \
  --outdir "$blastOutDir"

# Cutting out these genes/gene fragments from the genome sequences using the coordinates from the BLAST results
python3 "$scriptDir/blastProcessing.py" \
  --output-dir "$unalignedMultiFastas" \
  --ir-cutoff 4000 \
  --blast-dir "$blastResults"\
  --reference-dir "$referenceGenes"\
  --genome-dir "$plastidFastas"

# Aligning the sequences found as hits with BLAST
python3 "$scriptDir/aligner.py" \
  --input "$unalignedMultiFastas" \
  --output "$alignedMultiFastas"

# Updating the presence/absence TSV file and producing a change log file
python3 "$scriptDir/updateTSV.py" \
  --ogTSV "$presenceAbsenceTSV" \
  --alignDir "$alignedMultiFastas" \
  --outTSV "$updatedPresenceAbsenceTSV" \
  --changeLog "$changeLog"

# Making a new heatmap based on the updated TSV
python3 "$scriptDir/heatMapPlot.py" \
  "$updatedPresenceAbsenceTSV" \
  -o "$updatedHeatmap"