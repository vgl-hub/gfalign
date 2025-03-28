# gfalign

Graph alignment and analysis.

## Prerequisites
An installation of conda is required.

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

### Solve tangles in graph

Gfalign can find optimal paths between two nodes in assembly graphs. This can be useful particularly when trying to resolve complex tangles in genome assembly. To achieve this, you can run gfalign in `search` mode, e.g.:
```
gfalign filter -g all_reads.gaf -n nodelist.ls -o alignment.gaf --min-nodes 3 # it is good practice to remove reads that only align between two nodes
gfalign search -f assembly.gfa -n nodelist.tsv --source utig4-1 --destination utig4-2 -g alignment.gaf -m 100000000
```
`nodelist` is a two-column tab-separated file with all nodes in the tangle that the search cna visit and should try to include (i.e. for which the search maximizes for).  If a Hamiltonian path exists, it will be reported.

The output is a multi-column tab-separated file structured as follows:
```
path	# paths	#bad alignments	#good aligments	#diff	#score	#unique
```

Where `path` is the actual path, `# paths` is the number of paths from source to destination found so far, `#bad alignments` is the number of alignments that are incosistent with this path, `#good alignments` is the number of alignments that are fully consistent with the path, `diff` is `bad-good`, `score` is gfalign's score for the alignment, and `unique` is the number of unique nodes in the path.

Only paths that during the search improve on the number of nodes included in the path and in the good/bad alignments will be outputted.

Gfalign does only performs a node pseudoalignment based on the gaf alignment. Paths can be futher validated by realigning the reads to the path, e.g.:

```
cat assembly.gfa tangle.path.gfa > assembly.tangle.gfa # add a P line to the gfa with the path through the tangle
gfastats assembly.tangle.gfa -o assembly.tangle.fasta # extract the sequence in the tangle

rdeval reads.fastq.gz --include <(cut alignment.gaf) --homo 0 -o reads_subset.hc.fastq # extract reads for this region, homopolymer compress if the graph is hc
minimap2 -x map-ont -a -t 32 assembly.tangle.fasta reads_subset.hc.fastq | samtools sort -OBAM -o tangle.bam # align reads to the resolved tangle
samtools view -q1 -F 0x100 tangle.bam -o tangle.no-secondary.no-q0.bam # filter secondary alignments and q0 alignments
```

You can also align the node sequences to the path:

```
gfastats --include nodelist.ls assembly.gfa -o nodes.fasta --discover-paths # get fasta sequences for all nodes in the tangle
minimap2 -x map-hifi -a -t 32 assembly.tangle.fasta nodes.fasta | samtools sort -OBAM -o tangle_utig4-23_utig4-386.nodelist.utgs.bam # align nodes to path
```

Or align the path through the tangle back to the graph (or to another graph) with GraphAligner to inspect the quality of the alignment.

## How to cite

If you use **gfalign** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
