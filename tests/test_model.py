from pore_c2.model import ConcatemerMetaData


def test_align_data_md():
    d = ConcatemerMetaData.from_tags(["MI:Z:READ10", "Xc:B:i,10,300,1,10"])
    assert d == ConcatemerMetaData(
        concatemer="READ10", start=10, end=300, subread_idx=1, subread_total=10
    )
