#! /bin/bash
# mothurShared.sh
# Author: William L. Close
# Modified by: Lucas Bishop
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SAMPLEDIR=$1
export WORKDIR=$2
export REFDIR=$3



###################
# Run QC Analysis #
###################

# Making contigs from fastq.gz files, aligning reads to references, removing any non-bacterial sequences, calculating distance matrix, making shared file, and classifying OTUs
echo PROGRESS: Assembling, quality controlling, clustering, and classifying sequences.

mothur "#make.file(inputdir="${SAMPLEDIR}", outputdir="${WORKDIR}");
	make.contigs(file=current);
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference="${REFDIR}"/silva.v4.fasta);
	screen.seqs(fasta=current, count=current, start=1968, end=11550);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference="${REFDIR}"/trainset16_022016.pds.fasta, taxonomy="${REFDIR}"/trainset16_022016.pds.tax, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	dist.seqs(fasta=current, cutoff=0.03);
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03);
	classify.otu(list=current, count=current, taxonomy=current, label=0.03)"



# Renaming output files for use later
mv "${WORKDIR}"/*.precluster.pick.pick.fasta "${WORKDIR}"/errorinput.fasta
mv "${WORKDIR}"/*.vsearch.pick.pick.count_table "${WORKDIR}"/errorinput.count_table
mv "${WORKDIR}"/*.opti_mcc.shared "${WORKDIR}"/final.shared
mv "${WORKDIR}"/*.cons.taxonomy "${WORKDIR}"/final.taxonomy



####################################
# Make Group-Specific Shared Files #
####################################

# Sample shared file
echo PROGRESS: Creating sample shared file.

# Removing all mock and control groups from shared file leaving only samples
mothur "#remove.groups(shared="${WORKDIR}"/final.shared, groups=PCRwater_neg_1-PCRwater_neg_2-PCRwater_neg_3-empty_1-empty_2-empty_3-empty_4-empty_5-empty_6-empty_7-empty_8-Z_mock_ext-Zymo_mockpcr_1-Zymo_mockpcr_2-Zymo_mockpcr_3)"

# Renaming output file
mv "${WORKDIR}"/final.0.03.pick.shared "${WORKDIR}"/sample.final.shared



# Mock shared file
echo PROGRESS: Creating mock shared file.

# Removing non-mock groups from shared file
mothur "#get.groups(shared="${WORKDIR}"/final.shared, groups=z_mock_ext-Zymo_mockpcr_1-Zymo_mockpcr_2-Zymo_mockpcr_3)"

# Renaming output file
mv "${WORKDIR}"/final.0.03.pick.shared "${WORKDIR}"/mock.final.shared



# Control shared file
echo PROGRESS: Creating control shared file.

# Removing any non-control groups from shared file
mothur "#get.groups(shared="${WORKDIR}"/final.shared, groups=PCRwater_neg_1-PCRwater_neg_2-PCRwater_neg_3)"

# Renaming output file
mv "${WORKDIR}"/final.0.03.pick.shared "${WORKDIR}"/control.final.shared


