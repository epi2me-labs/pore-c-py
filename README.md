# Pore_c_py
This package provides python scripts for working with Pore-C data

Use `--help` with to find detailed usage instructions.

## Installation

Via pip:
```
pip install pore-c-py
```
Or via conda:
```
conda install -c epi2me-labs pore-c-py
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

### Examples

*Digest concatemer sequences in to monomers:*
```
pore-c-py digest <input concatemers bam> <restriction enzyme> <output monomers bam file name>
```

*Parse a bam of aligned monomers:*
```
pore-c-py annotate <input bam> <output bam file name>
```

