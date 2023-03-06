"""Conf test."""
# flake8: noqa

# Fixtures are basically test data, you see some of these things
# below as arguments to test functions. We made need equivalents
# of some of these. Though I suggest keeping fixtures with the
# tests in which they are used.


# @pytest.fixture(scope="session")
# def mock_read_no_qual() -> ReadSeq:
#     """Mock read no qual."""
#     return ReadSeq(name="read_no_qual", sequence="ACTG", quality=None)
#
#
# @pytest.fixture(scope="session")
# def mock_read_w_qual() -> ReadSeq:
#     """Mock read with qual."""
#     return ReadSeq(name="read_w_qual", sequence="ACTGACTG", quality="!" * 8)
#
#
# @pytest.fixture(scope="session")
# def mock_read_w_tags() -> ReadSeq:
#     """Mock read with tags."""
#     return ReadSeq(
#         name="read_w_tags",
#         sequence="AACGTTCGAAC",
#         quality="!!00{}22[]]",
#         tags={
#             "RG": "RG:Z:RG01",
#             "Mm": "Mm:Z:C+m,0,1;",
#             "Ml": "Ml:B:C,122,128",
#         },
#     )
#
#
# @pytest.fixture(scope="session")
# def mock_read_w_rg() -> ReadSeq:
#     """Mock read w_rg."""
#     return ReadSeq(
#         name="read_w_rg",
#         sequence="AACGTTCGAAC",
#         quality="!!00{}22[]]",
#         tags={
#             "RG": "RG:Z:RG01",
#         },
#     )
#
#
# @pytest.fixture(scope="session")
# def mock_reads(
#     mock_read_no_qual: ReadSeq,
#     mock_read_w_qual: ReadSeq,
#     mock_read_w_tags,
#     mock_read_w_rg,
# ) -> Dict[str, ReadSeq]:
#     """Mock reads."""
#     reads = [
#         mock_read_no_qual, mock_read_w_qual, mock_read_w_tags, mock_read_w_rg]
#     return {r.name: r for r in reads}
#
#
# @pytest.fixture
# def monomer_read_seqs():
#     """Monomer read sequences."""
#     num_concatemers, aligns_per_concatemer = 2, 10
#     res = []
#     for concat_idx in range(num_concatemers):
#         concatemer_id = f"CONCAT{concat_idx}"
#         for monomer_idx in range(aligns_per_concatemer):
#             start = monomer_idx * 10
#             end = start + 10
#             monomer_id = f"{concatemer_id}:{start}:{end}"
#             m = MonomerReadSeq(
#                 concatemer_id=concatemer_id,
#                 monomer_id=monomer_id,
#                 coords=ConcatemerCoords(
#                     start=start,
#                     end=end,
#                     read_length=aligns_per_concatemer * 10,
#                     subread_idx=monomer_idx,
#                     subread_total=aligns_per_concatemer,
#                 ),
#                 read_seq=ReadSeq(
#                     name=monomer_id,
#                     sequence="A" * 10,
#                     quality="!" * 10,
#                     align_info=AlignInfo(
#                         ref_name=f"chr{monomer_idx + 1}",
#                         ref_pos=100 * monomer_idx,
#                         cigar="10M",
#                         map_quality=20,
#                         length=10,
#                     ),
#                     tags={MOLECULE_TAG: f"{MOLECULE_TAG}:Z:{concatemer_id}"},
#                 ),
#             )
#             res.append(m)
#     return res
#