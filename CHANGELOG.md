# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.1.1]
### Changed
- Updated requirements of conda package.

## [v2.1.0]
### Added
- Option `max_monomers` to `digest` command to remove reads containing more
  than max monomers. Excluded reads can be output to second BAM file and/or
  a text file of read names.
### Removed
- `--max_reads` argument from `digest` command, as it was non-functional. 

## [v2.0.6]
### Fixed
- Chromunity parquet file output contains columns required for use with Chromunity.

## [v2.0.5]
### Changed
- Updated LICENSE to Oxford Nanopore Technologies PLC. Public License Version 1.0.

## [v2.0.4]
### Fixed
- Chromunity writer group colinear output error.
- Modified bases digest step and associated tests.

## [v2.0.3]
### Added
- Digest sub command will accept stdin input.

## [v2.0.2]
### Added
- Threads parameter to chunk_bam sub command.

## [v2.0.1]
### Fixed
- Annotate sub command check of existing files.

## [v2.0.0]
### Added
- Stdout options for annotate and digest sub commands.
### Changed
- Simplified logging for easier parsing.
- Rewrote CLI, removing utils and keeping only those programs which are used.
- Parse bam subcommand renamed annotate.
### Removed
- A lot of code and dependencies.

## [v1.0.0]
### Changed
- Tidy up for release as a conda package in the EPI2MELabs channel 

## [v0.7.2]
### Fixed
- digest: don't strip mod base tags by default

## [v0.7.1]
### Added
- create-chunked-ubam: break an unaligned bam into chunks

## [v0.7.0]
### Added
- bam: add PG tag to bam files produced
- export: add separate export subcommand
- paired-end-export: enable filtering paired end read export to cis only and distance bounds
- chromunity: add option to merge fragments colinear in the read and genome into a single 'monomer'
- cli.py: add option to remove SAM tags from input, rename and move subcommands
- digest: move digest to top-level command, add option to remove sam tags from output
- aligns.py: add functions to group alignments colinear in the genome
### Fixed
- model.py: fix mixed up sam tags
### Changed
- remove kw_only from dataclass (requires python 3.10 or higher)
- overlaps.py: remove fragment overlap code
- dataclasses: replace attrs uses with dataclasses
- mappy: remove genome index and mappy uses

## [v0.6.0]
### Added
- io.py: write basic summary stats in json format
- io.py: write a chromunity-compatible parquet file
### Fixed
- model.py: use read start,end coords as monomer id
- Snakefile: fix integration tests

## [v0.5.1]
### Fixed
- aligns.py: fix bug caused by change in lru_cache signature
### Changed
- sam_tags.py: rename to sam_utils to reflect broader functoin
- utils.py: remove utils.py put functionaliy into io and sam_tags
- sam_tags.py: move sam tag processing to separate module
- model.py: make pore-s-specific tags into global settings

## [v0.5.0]
### Fixed
- variants.py: remove partial phasing code
- testing.py: fix phased vcf generation

## [v0.4.0]
### Changed
- model.py: set paired end fields and add walk metadata
- model.py: downgrade mm tags to work with pysam, ease roundtrip of tag data
- cli.py: analyse the first n reads
### Fixed
- hatch.toml: fix type

## [0.3.0]
### Fixed
- model.py: handle case where read is 1-base long and has '*' as quak
### Changed
- refactor model.py

## [0.2.0]
### Added
- aligns.py: Further annotation of alignments
- aligns.py: extract pairwise alignments
- aligns.py: sort alignments by concatemer_index
- cli.py: add cli tool to generate test data
- align.py: group algnments by concatemer id
- porec: refactor to use generators to avoid memory issues
- utils.py: code to handle sam flags
- io.py: write aligments to sorted bam
- log.py: remove loguru dependency add remora-inspired logging setup
- map.py: find fragment overlaps for alignments
- map.py: add data structure to hold alignments for a concatemer
- index.py: move indexing code to separate file
- test_digest: initial implementation of virtual digest
- Initial-commit: iniital commit
### Fixed
- monomers.py: fix edge case for empty cuts
- pyproject.toml: fix broken toml
### Changed
- pyproject.toml: remove testw command, can't get it to work
- model.py: move some method logic to functions
- map.py: add multithreaded alignment
- cli.py: update mapping cli
- map.py: move test code into place
- digest.py: use filecollection to manage multiple output
- digest.py: refactor to support mapping
- digest.py: dataclass and pandera schema for digest
- digest.py: move functions out of test code
- pyscaffold: remove old pyscaffold files
