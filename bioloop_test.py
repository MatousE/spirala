import bioloopcalc as bio


def test_calculate_cg():
    assert bio.calculatecg('CATGCATG') == '50.0'

