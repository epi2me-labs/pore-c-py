# Porecpy
This package provides python scripts for working with Pore-c data

Use `--help` with to find detailed usage instructions.

## Installation
Via pip:
```
pip install porecpy
```
Or via conda:
```
conda install -c epi2me-labs porecpy
```

## Usage

```
$ porecpy --help
usage: porecpy [OPTIONS] COMMAND [ARGS].

Available subcommands are:
    digest       Digest concatemer sequences into monomers using a restriction enzyme.      
    parse-bam    Parse a BAM file of aligned monomers.
```

### Examples

*Digest concatemer sequences in to monomers:*
```
porecpy digest <input concatemers bam> <restriction enzyme> <output monomers bam file name>
```

*Parse a bam of aligned monomers:*
```
porecpy parse-bam <input bam> <output bam file name>
```

