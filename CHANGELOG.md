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
