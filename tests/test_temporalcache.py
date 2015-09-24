import TemporalCache


def test_basic_temporal_cache():
    tc=TemporalCache.TemporalCache('tests/data/circle')
    assert len(tc.data) == 2


def test_temporal_cache_with_range():
    tc=TemporalCache.TemporalCache('tests/data/circle',0,1)
    assert len(tc.data) == 2
    assert tc.lower==0
    assert tc.upper==1
