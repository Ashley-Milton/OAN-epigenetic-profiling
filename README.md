# OAN-epigenetic-profiling

## General overview

This repository contains scripts used to process and visualise ChIP-seq, Hi-C and RNA-seq data.

insert image here

## ChIP-seq
##### Chromatin ImmunoPrecipitation sequencing

ChIP-seq data processing scripts available here include those for quality control (FastQC and MultiQC), read trimming (Trimmomatic), read mapping to a reference genome (Subread), and ChIP peak calling (macs2). 

ChIP-seq data visualisation/analysis scripts available here include those for intersecting ChIP data with transcription start site (TSS) information, generating TSS metagene plots, calculating ChIP peak densities and performing non-parametric statistical tests.

## Hi-C
##### High-throughput Chromosome conformation capture 

Hi-C data processing scripts available here include those for quality control (FastQC and MultiQC), read trimming (Trimmomatic), reference genome indexing (bwa), read mapping to a reference genome (bwa) and Hi-C matrix construction and processing (hicexplorer).

Hi-C data visualisation/analysis scripts available here include those for generating distance-dependent Hi-C contact plots, counting intra- and inter-chromosomal interactions, generating inter-chromosomal interaction plots and generating a top ten chromosome pair plot.

## RNA-seq
##### RNA sequencing

RNA-seq data processing scripts available here include those for quality control (FastQC and MultiQC), reference genome indexing (Subread), read mapping to a reference genome (Subread) and feautre counting (Subread).

RNA-seq visualisation/analysis script is for processing, normalising and visualising gene count data.
