from loguru import logger


def test_map(scenario):
    logger.debug(scenario)
    fastq = scenario.concatemer_fastq
    logger.debug(fastq)
    # res = simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme='EcoRI')
    # logger.debug(res)
