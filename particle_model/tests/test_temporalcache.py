import particle_model.TemporalCache as TC
import vtk


def test_basic_temporal_cache():
    tc=TC.TemporalCache('particle_model/tests/data/circle')
    assert len(tc.data) == 3


def test_temporal_cache_with_range():
    tc=TC.TemporalCache('particle_model/tests/data/circle',0,1)
    assert len(tc.data) == 3
    assert tc.lower==0
    assert tc.upper==1

    count=0
    for d in tc.data:
        if d[2]: count+=1

    assert tc.lower==0
    assert tc.upper==1

    print tc.data

    assert count==2


def test_temporal_cache_new_range():
    tc=TC.TemporalCache('particle_model/tests/data/circle',0,1)

    tc.range(6,8)

    assert tc.lower==1
    assert tc.upper==2

    count=0
    for d in tc.data:
        if d[2]: count+=1

    print tc.data

    assert count==2



def test_temporal_cache_call():

    tc=TC.TemporalCache('particle_model/tests/data/circle',0,1)

    d,alpha=tc(0.5)

    print d,alpha


    assert d[0][0]==0.0
    assert d[1][0]-5.0 < 1.e-8
    assert alpha-0.1<1e-8
    
    
