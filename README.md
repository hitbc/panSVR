PanSVR
======

PanSVR: pan-genome augmented short read realignment for sensitive detection of structural variations

## Table of Contents
1. [Dependency](#dependency)
2. [Quick start](#Quick-start)
3. [Introduction](#Introduction)
4. [Memory usage](#Memory-usage)
5. [Build project](#build-project)
6. [Simple run](#build-index)
7. [Run using self-defined parameters](#run-classifation)
8. [Run analysis](#run-analysis)
9. [Demo data](#Demo-data)

## Dependency

Building of PanSVR depends on g++, make and zlib.
Running of PanSVR depends on samtools and pigz.

Run following commands to install them (Ubuntu).
```
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install zlib1g-dev
sudo apt-get install samtools
sudo apt-get install pigz

```
## Quick start
```
git clone https://github.com/hitbc/PanSVR.git
cd ./PanSVR
bash ./build
#SV calling
./PanSVR classify -t 4 ./demo_index ./demo/ERR1050068.fastq -o ./ERR1050068.sam
```

## Introduction

PanSVR (Pan-genome augmented Structure Variation calling tool with read Re-alignment), is a novel pan-genome-based SV calling approach. PanSVR uses several tailored methods to implement precise re-alignment for SV-spanning reads against well-organized pan-genome reference with plenty of known SVs. PanSVR enables to greatly improve the quality of short read alignments and produce clear and homogenous SV signatures which facilitate SV calling. Benchmark results on real sequencing data suggest that PanSVR is able to largely improve the sensitivity of SV calling than that of state-of-the-art SV callers, especially for the SVs from repeat-rich regions and/or novel insertions which are difficult to existing tools.

## Memory usage

Normally, PanSVR used less than 4 Gigabytes memory in all steps. Besides, 100 Gigabytes space in hard disk is needed to store tmp files. 

## Build project

To build PanSVR
```
git clone https://github.com/hitbc/PanSVR.git
cd ./PanSVR
bash ./build
```
