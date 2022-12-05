## 0.7.2 (2022-12-05)

### Fix

- **digest**: don't strip mod base tags by default

## 0.7.1 (2022-12-02)

### Feat

- **create-chunked-ubam**: break an unaligned bam into chunks

## 0.7.0 (2022-12-01)

### Feat

- **bam**: add PG tag to bam files produced
- **export**: add separate export subcommand
- **paired-end-export**: enable filtering paired end read export to cis only and distance bounds
- **chromunity**: add option to merge fragments colinear in the read and genome into a single 'monomer'
- **cli.py**: add option to remove SAM tags from input, rename and move subcommands
- **digest**: move digest to top-level command, add option to remove sam tags from output
- **aligns.py**: add functions to group alignments colinear in the genome

### Fix

- **model.py**: fix mixed up sam tags

### Refactor

- remove kw_only from dataclass (requires python 3.10 or higher)
- **overlaps.py**: remove fragment overlap code
- **dataclasses**: replace attrs uses with dataclasses
- **mappy**: remove genome index and mappy uses

## 0.6.0 (2022-11-15)

### Feat

- **io.py**: write basic summary stats in json format
- **io.py**: write a chromunity-compatible parquet file

### Fix

- **model.py**: use read start,end coords as monomer id
- **Snakefile**: fix integration tests

## 0.5.1 (2022-10-31)

### Fix

- **aligns.py**: fix bug caused by change in lru_cache signature

### Refactor

- **sam_tags.py**: rename to sam_utils to reflect broader functoin
- **utils.py**: remove utils.py put functionaliy into io and sam_tags
- **sam_tags.py**: move sam tag processing to separate module
- **model.py**: make pore-s-specific tags into global settings

## 0.5.0 (2022-10-21)

### Fix

- **variants.py**: remove partial phasing code
- **testing.py**: fix phased vcf generation

## 0.4.0 (2022-10-15)

### Feat

- **model.py**: set paired end fields and add walk metadata
- **model.py**: downgrade mm tags to work with pysam, ease roundtrip of tag data
- **cli.py**: analyse the first n reads

### Fix

- **hatch.toml**: fix type

## 0.3.0 (2022-10-11)

### Fix

- **model.py**: handle case where read is 1-base long and has '*' as quak

### Refactor

- refactor model.py

## 0.2.0 (2022-10-04)

### Feat

- **aligns.py**: Further annotation of alignments
- **aligns.py**: extract pairwise alignments
- **aligns.py**: sort alignments by concatemer_index
- **cli.py**: add cli tool to generate test data
- **align.py**: group algnments by concatemer id
- **porec**: refactor to use generators to avoid memory issues
- **utils.py**: code to handle sam flags
- **io.py**: write aligments to sorted bam
- **log.py**: remove loguru dependency add remora-inspired logging setup
- **map.py**: find fragment overlaps for alignments
- **map.py**: add data structure to hold alignments for a concatemer
- **index.py**: move indexing code to separate file
- **test_digest**: initial implementation of virtual digest
- **Initial-commit**: iniital commit

### Fix

- **monomers.py**: fix edge case for empty cuts
- **pyproject.toml**: fix broken toml

### Refactor

- **pyproject.toml**: remove testw command, can't get it to work
- **model.py**: move some method logic to functions
- **map.py**: add multithreaded alignment
- **cli.py**: update mapping cli
- **map.py**: move test code into place
- **digest.py**: use filecollection to manage multiple output
- **digest.py**: refactor to support mapping
- **digest.py**: dataclass and pandera schema for digest
- **digest.py**: move functions out of test code
- **pyscaffold**: remove old pyscaffold files
