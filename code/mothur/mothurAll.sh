#! /bin/bash
# mothurAll.sh
# William L. Close
# modified by: Lucas Bishop
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SAMPLEDIR=$1
export WORKDIR=$2
export REFDIR=$3



####################
# Running Pipeline #
####################

# # Replacing hyphens in group names to prevent errors when running mothur later
# echo PROGRESS: Fixing sample names by replacing hypens.

# # Replacing hyphens with the letter "N" for "Negative"
# bash code/mothur/mothurNames.sh data/mothur/raw/ N


# Running the mothur pipeline
echo PROGRESS: Running mothur.

# Generating data files
bash code/mothur/mothurShared.sh "${SAMPLEDIR}" "${WORKDIR}" "${REFDIR}"
bash code/mothur/mothurCountShared.sh "${WORKDIR}"
bash code/mothur/mothurSubsampleShared.sh "${WORKDIR}"

# Analyzing data files
bash code/mothur/mothurError.sh "${WORKDIR}" "${REFDIR}"
bash code/mothur/mothurRarefaction.sh "${WORKDIR}"/sample.final.shared
bash code/mothur/mothurAlphaBeta.sh "${WORKDIR}"/sample.final.shared
bash code/mothur/mothurPCoA.sh "${WORKDIR}"/sample.final.braycurtis.0.03.lt.ave.dist
bash code/mothur/mothurNMDS.sh "${WORKDIR}"/sample.final.braycurtis.0.03.lt.ave.dist

# Cleaning up
mkdir -p "${WORKDIR}"/intermediate/
mv "${WORKDIR}"/$(ls -p "${WORKDIR}" | grep -v "/" | grep -v "final") "${WORKDIR}"/intermediate/





