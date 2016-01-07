""" Unit tests for the temporal cache."""
import particle_model.TemporalCache as TC

def test_basic_temporal_cache():
    """Test if a simple cache will create itself."""
    temp_cache = TC.TemporalCache('particle_model/tests/data/circle')
    assert len(temp_cache.data) == 3

def test_temporal_cache_with_range():
    """Test if we can create a range."""
    temp_cache = TC.TemporalCache('particle_model/tests/data/circle', 0, 1)
    assert len(temp_cache.data) == 3
    assert temp_cache.lower == 0
    assert temp_cache.upper == 1

    count = 0
    for data in temp_cache.data:
        if data[2]:
            count += 1

    assert temp_cache.lower == 0
    assert temp_cache.upper == 1
    assert count == 2


def test_temporal_cache_new_range():
    """Test if we can reset the range."""
    temp_cache = TC.TemporalCache('particle_model/tests/data/circle', 0, 1)

    temp_cache.range(6, 8)

    assert temp_cache.lower == 1
    assert temp_cache.upper == 2

    count = 0
    for data in temp_cache.data:
        if data[2]:
            count += 1

    assert count == 2



def test_temporal_cache_call():
    """Test if we can make a call."""
    temp_cache = TC.TemporalCache('particle_model/tests/data/circle', 0, 1)

    data, alpha, names = temp_cache(0.5)

    assert data[0][0] == 0.0
    assert data[1][0]-5.0 < 1.e-8
    assert alpha-0.1 < 1e-8
