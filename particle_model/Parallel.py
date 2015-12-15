""" Module containing helper routines for use with MPI parallelisation."""

from mpi4py import MPI

def isparallel():
    """ Check if this is a parallel run."""
    comm=MPI.COMM_WORLD

    return comm.Get_size()>1

def get_rank():
    """ Get MPI rank."""
    comm=MPI.COMM_WORLD

    return comm.Get_rank()

def get_size():
    """ Get MPI size. """
    comm=MPI.COMM_WORLD

    return comm.Get_size()

def get_world_comm():
    """ Get the world communicator."""
    comm=MPI.COMM_WORLD

    return comm
