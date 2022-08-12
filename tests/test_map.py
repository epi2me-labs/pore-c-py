from pathlib import Path

import pytest

from pore_c2.cli import app
from pore_c2.index import create_index
from pore_c2.map import map_concatemers
from pore_c2.testing import Scenario


def test_mappy(scenario: Scenario, tmp_path):

    index_files = create_index(
        fasta=scenario.reference_fasta,
        enzyme=scenario.enzyme,
        prefix=Path(tmp_path / "index"),
    )
    index_metadata = index_files.load_metadata()

    map_concatemers(
        enzyme=index_metadata.enzyme,
        fastq=scenario.concatemer_fastq,
        mmi=index_files.mmi,
        minimap_settings=index_metadata.mappy_settings,
        fragment_pq=index_files.fragments,
    )
    raise ValueError(index_files)


@pytest.mark.skip
def test_map_cli(runner, scenario, tmp_path):
    outfile = tmp_path / "results.pq"
    result = runner.invoke(
        app,
        [
            "map",
            str(scenario.concatemer_fastq),
            scenario.enzyme,
            scenario.reference_fasta,
            str(outfile),
        ],
    )
    assert result.exit_code == 0


# @pytest.mark.skip
# def test_map(scenario):
#    logger.debug(scenario)
#    ref_fasta = scenario.reference_fasta
#    idx_file = Path(str(ref_fasta) + ".bwt")
#    if not idx_file.exists():
#        sp.run(["bwa", "index", str(ref_fasta)], capture_output=True, check=True)
#        # logger.debug(list(ref_fasta.parent.glob("*.*")))
#    fastq = scenario.concatemer_fastq
#    threads = 1
#    bam_file = fastq.with_suffix(".sam")
#    sp.run(
#        f"bwa bwasw -t {threads} -b 5 -q 2 -r 1 -T 15 -z 10 "
#        f"{ref_fasta} {fastq} | samtools sort -n -o {bam_file} -",
#        shell=True,
#        capture_output=True,
#    )
#
#    ol = FragmentOverlapper.from_fragments_df(scenario.fragments_df)
#    af = AlignmentFile(str(bam_file))
#    for query_name, aligns in groupby(af, lambda x: x.query_name):
#        aligns = list(aligns)
#        concatemer_id, fragments = query_name.split(":")
#        fragments = [int(_) for _ in fragments.split("_")]
#        d = ConcatemerAlignData(
#            concatemer_id, sorted(aligns, key=lambda x: x.query_alignment_start)
#        )
#        for a in d.raw_alignments:
#            logger.debug(
#                f"{a.query_name}:{a.query_alignment_start}-{a.query_alignment_end} "
#                f"{a.reference_start}-{a.reference_end}"
#            )
#            try:
#                olap = ol.find_overlap(
#                    a.reference_name, a.reference_start, a.reference_end
#                )
#            except ZeroDivisionError:
#                print(a)
#    # res = simulate_fasta_with_cut_sites(fasta, lengths, expected, enzyme='EcoRI')
#    # logger.debug(res)
