""" Unit test for options."""

import libspud

import numpy
from particle_model import Options


def test_options_reader():

    reader = Options.OptionsReader("particle_model/tests/data/test.flml")

    assert reader.get_name() == "python_state_to_vtk_test"
    assert all(reader.get_gravity() == [0.,0.,0.])
    assert all(reader.get_rotation()[0] == numpy.zeros(3))
    assert all(reader.get_rotation()[1] == numpy.zeros(3))
    assert reader.get_outlet_ids() == None
    assert reader.get_inlets() == []
    assert reader.get_mesh_filename() == "Unstructured.msh"
    assert reader.get_current_time() == 1.0
    assert reader.get_timestep() == 1.0
    assert reader.get_finish_time() == 1.0
    assert reader.get_adapts_at_first_timestep() == 0
    assert reader.get_dump_period()[0] == 0.0
    assert reader.get_dump_period()[1] == False
    

