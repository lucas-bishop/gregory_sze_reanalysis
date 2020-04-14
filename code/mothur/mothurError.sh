#! /bin/bash
# mothurError.sh
# William L. Close
# modified by Lucas Bishop
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export WORKDIR=$1
export REFDIR=$2



######################
# Run Error Analysis #
######################

# Calculating error rates compared to Zymo reference sequences
echo PROGRESS: Calculating error rate.

mothur "#get.groups(fasta="${WORKDIR}"/errorinput.fasta, count="${WORKDIR}"/errorinput.count_table, groups=Mock1-Mock2-Mock3-Mock4);
	seq.error(fasta=current, count=current, reference="${REFDIR}"/zymo_mock.align, aligned=F)"



# Moving error analysis files to error directory
echo PROGRESS: Storing error logs.

mkdir -p "${WORKDIR}"/error_analysis

mv "${WORKDIR}"/errorinput.* "${WORKDIR}"/error_analysis/


