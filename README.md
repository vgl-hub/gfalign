# gfalign

Graph alignment and analysis.

## Installation

To install run `git clone https://github.com/vgl-hub/gfalign.git --recursive` and `make -j` in the `gfalign` folder.

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
