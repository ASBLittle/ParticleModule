""" Test parallel execution."""
from particle_model import Parallel
import pytest

@pytest.mark.skipif(Parallel.get_size()>1,reason='Serial test')
def test_serial():

    assert not Parallel.is_parallel()


@pytest.mark.skipif(Parallel.get_size() == 1,reason='Parallel test')
def test_paralle():

    assert Parallel.is_parallel()
