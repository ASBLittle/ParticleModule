""" Module containing helper routines for use with MPI parallelisation."""

try:
    from mpi4py import MPI
except:
    MPI = None

import numpy
import copy
import itertools 

def is_parallel():
    """ Check if this is a parallel run."""

    if MPI is None:
        return False
    
    comm = MPI.COMM_WORLD

    return comm.Get_size()>1

def barrier():
    """ Set up an MPI barrier."""

    if MPI is None:
        return None

    comm = MPI.COMM_WORLD
    
    comm.Barrier()

    return None

def get_rank():
    """ Get MPI rank."""

    if MPI is None:
        return 0
    
    comm = MPI.COMM_WORLD

    return comm.Get_rank()

def get_size():
    """ Get MPI size. """

    if MPI is None:
        return 1
    
    comm = MPI.COMM_WORLD

    return comm.Get_size()

def get_world_comm():
    """ Get the world communicator."""
    comm = MPI.COMM_WORLD

    return comm

def is_root(root=0):
    """ Is this process rank 0?"""

    return get_rank()==root
    

def point_in_bound(pnt, bound):
    """Check whether a point is inside the bounds""" 
    if pnt[0]<bound[0]:
        return False
    if pnt[0]>bound[1]:
        return False
    if pnt[1]<bound[2]:
        return False
    if pnt[1]>bound[3]:
        return False
    if pnt[2]<bound[4]:
        return False
    if pnt[2]>bound[5]:
        return False
    return True

def gather_bounds(bounds):
    """ Exchange bounds across multiple processors """
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    all_bounds = np.empty([size, 6], dtype=bounds.dtype)

    comm.Allgather(bounds, all_bounds)

    return all_bounds

def distribute_particles(particle_list, system, time=0.0):
    """ Handle exchanging particles across multiple processors """

    try:
        block = system.temporal_cache.block
    except:
        data = system.temporal_cache(time)
        block = data[0][1][-2]
    bounds = system.temporal_cache.get_bounds(time)
    

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if not is_parallel():
        return set(particle_list)

    all_bounds = numpy.empty([size, 6], dtype=float)

    comm.Allgather(bounds, all_bounds)

    odata=[]

    for i in range(size):
        if i==rank:
            odata.append([])
            continue
        plist=[]
        for par in particle_list:
            if not point_in_bound(par.pos, all_bounds[i]):
                continue
            par=par.copy()
            plist.append(par)
        odata.append(plist)

    data = comm.alltoall(sendobj=odata)

    output=set(particle_list)
    for i in range(size):
        if i==rank: continue

        for par in data[i]:
            par.system = system
            output.add(par)

    b = system.particle_in_system(output, time, rank)
    output =  list(itertools.compress(output,
                                   b))

    return output

newid = itertools.count().next

class particle_id(object):

    def __init__(self,hash=None):
        global newid
        if hash is None:
            self.creator_id = newid()
            self.creator_rank=get_rank()
        else:
            self.creator_id, self.creator_rank = divmod(hash, get_size())
            if self.creater_id> particle_id.newid():
                ### we should really update the id creator here
                newid = itertools.count(self.creater_id+1).next
            

    def __call__(self):
        return self.creator_id*get_size()+self.creator_rank

    def __hash__(self):
        return self()

    def __eq__(self,obj):
        return self()==obj()


def cell_owned(block, ele):
    if not is_parallel():
        return True
    if not block.IsA("vtkMultiblockDataSet"):
        return True
    if not block.GetBlock(0).GetCellData().HasArray("ElementOwner"):
        return True
    return  block.GetBlock(0).GetCellData().GetScalar("ElementOwner").GetValue(ele) == get_rank()+1


def point_owned(block,X):
    out = numpy.empty(X.shape[0],bool)
    for k, x in enumerate(X):
        ele = find_cell(block, x)
        if ele > -1:
            out[k] = cell_owned(block, ele)
        else:
            out[k] = False
    
