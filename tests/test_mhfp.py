def test_calc_mhfp_smoke():
    from mqn_search.mhfp_512 import calc_mhfp
    fp = calc_mhfp("CCO")
    assert fp is not None

def test_mhfp_consistent():
    from mqn_search.mhfp_512 import check_env
    consistent = check_env()
    assert consistent is True
