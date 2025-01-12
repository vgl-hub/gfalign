# gfalign

Graph alignment and analysis.

## Installation

To install run `git clone https://github.com/vgl-hub/gfalign.git --recursive` and `make -j` in the `gfalign` folder.
If you wish to generate alignments using gfalign (`align` option), you need to have GraphAligner installed. You can either install [GraphAligner](https://github.com/maickrau/GraphAligner) directly (recommended) or run `make GraphAligner` (requires Conda).

## Usage

`gfalign [options] [tool] [arguments]`

Gfalign uses a graph-aware aligner, e.g. GraphAligner, to align reads to an assembly graph and then evaluates the graph based on the alignment.

First generate an alignment, e.g.:

`gfalign -p hifi align -f testFiles/random2.reads.fq -g testFiles/random2.gfa -a aln.gaf`

Then compute statistics:

`gfalign eval -g aln.gaf`

Or decorate the graph with information from the alignment:

`gfalign eval -g aln.gaf -f testFiles/random2.gfa -o newGFA.gfa`

To check out all tool and options available use `gfalign -h`.

## How to cite

If you use **gfalign** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
