"""Test drag models."""
import particle_model.DragModels as DM
from numpy import array

def test_stokes_drag():
    """Test Stokes drag law"""
    assert DM.stokes_drag(fluid_velocity=1.0, particle_velocity=0.0,
                          diameter=1.0, rho=1.0, fluid_viscosity=1.0) == 18.0


def test_turbulent_drag():
    """Test turbulent drag law"""
    fluid_vel = array((1., 0., 0.))
    vel = array((0., 0., 0.))
    emax = array((1.0e-8, 1.0e-8, 1.0e-8))
    drag = 0.44 * 3.0 / 32.0 * fluid_vel
    assert all(DM.turbulent_drag(fluid_velocity=fluid_vel, particle_velocity=vel,
                                 diameter=1.0, rho=1.0, rho_f=1.0, fluid_viscosity=1.0) - drag < emax)

    fluid_vel = array((0., 5., 0.))
    vel = array((0., 1., 0.))
    drag = 0.44 * 3.0 / 20.0 * vel

    assert all(DM.turbulent_drag(fluid_velocity=fluid_vel, particle_velocity=vel,
                                 diameter=10.0, rho=1.0, rho_f=1.0, fluid_viscosity=1.0) - drag < emax)
