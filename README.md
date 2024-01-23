# Pore_c_py
This package provides Python scripts for working with Pore-C data. It
is not intended to be used directly by end uses, but rather as part of
our Nextflow workflow [wf-pore-c](https://github.com/epi2me-labs/wf-pore-c).
Hence the terse nature of this documentation.

## Installation

A package is available to install through either `pip`:
```
pip install pore-c-py
```
or through conda:
```
conda install -c nanoporetech pore-c-py
```

## Usage

```
$ pore-c-py --help
usage: pore_c_py [OPTIONS] COMMAND [ARGS].

Available subcommands are:
    digest       Digest concatemer sequences into monomers using a restriction enzyme.      
    annotate     Annotate alignments with a "walk", which simply
                 enumerates the alignment coordinates of the monomers comprising the
                 concatemer.
```

### Example

The following is indicative use, similar to that performed by
[wf-pore-c](https://github.com/epi2me-labs/wf-pore-c).

```
INPUT="myreads.bam"
ENZYME="NlaIII"
REF="myref.fasta"
OUTPUT="all"

pore-c-py digest "${INPUT}" "${ENZYME}" \
    | samtools fastq -T '*' \
    | minimap2 -ay -t 8 -x map-ont "${REF}" - \
    | pore-c-py annotate - "${OUTPUT}" --monomers --stdout --summary --chromunity \
    | tee "${OUTPUT}.ns.bam" \
    | samtools sort --write-index -o "${OUTPUT}.cs.bam" -
samtools index "${OUTPUT}.ns.bam"
```

The `digest` program can read its input from standard input, so it can be used with
[bamindex](https://github.com/epi2me-labs/fastcat?tab=readme-ov-file#bamindex) in
order to process a subset of a file. This is particularly useful for
distributing the workload on a cluster.
