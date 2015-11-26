""" Module deadling with initialising the code from an xml format options file using spud"""


import libspud
import numpy

class OptionsReader(object):
    """ Handle processing the XML into python"""
    def __init__(self,options_file=None):

        if options_file:
            libspud.load_options(options_file)


            self.simulation_name=libspud.get_option('/simulation_name')
            self.dimension=libspud.get_option('/geometry/dimension')


    def get_name(self):
        """ Return the simulation name (for output/cross communication with fluidity)."""
        return libspud.get_option('/simulation_name')

    def get_gravity(self):
        """ Return the acceleration due to gravity."""
        options_base='/physical_parameters/gravity'
        magnitude=0
        direction=numpy.zeros(3)

        if libspud.have_option(options_base):
            magnitude=libspud.get_option(options_base+'/magnitude')
            direction[:self.dimension]=libspud.get_option(options_base
                        +'/vector_field[0]/prescribed/value[0]/constant')

        return magnitude*direction

    def get_rotation(self):
        """ Return the rotation vector."""
        options_base='/physical_parameters/coriolis/specified_axis'
        magnitude=0
        omega=numpy.zeros(3)
        origin=numpy.zeros(3)

        if libspud.have_option(options_base):
            omega[2]=libspud.get_option(options_base+'/rotational_velocity')
            origin[:self.dimension]=libspud.get_option(options_base+'point_on_axis')

        return omega, origin

    def get_model_option(self, option_name):
        """ interogate the model specific options """

        options_base = '/embedded_models/particle_model/'
        if libspud.have_option:
            return libspud.get_option(options_base+option_name)
        else:
            return None 

    def get_mesh_filename(self):
        """Return the mesh file name"""
        return libspud.get_option('/geometry/mesh::CoordinateMesh/from_file/file_name')+'.msh'
    
    def get_current_time(self):
        """Return the current time from the options file."""
        return libspud.get_option('/timestepping/timestep')

    def get_timestep(self):
        """Return the timestep from the options file."""
        return libspud.get_option('/timestepping/timestep')

    def get_finish_time(self):
        """Return the finish time from the options file."""
        return libspud.get_option('/timestepping/finish_time')
