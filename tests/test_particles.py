import Particles

def test_tests():
    assert 1


def test_basic_particle_initialization():
    from numpy import zeros

    p=zeros(3)
    v=zeros(3)

    pt=Particles.particle(p,v)

    assert all(pt.p==p) and all(pt.v==v)



def test_basic_particle_bucket_initialization(tmpdir):
    from numpy import zeros

    N=10

    p=zeros((N,3))
    v=zeros((N,3))

    p=Particles.particle_bucket(p,v)


