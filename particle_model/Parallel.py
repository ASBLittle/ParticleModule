""" Module containing helper routines for use with MPI parallelisation."""

try:
    from mpi4py import MPI
except:
    MPI = None

import itertools 

def is_parallel():
    """ Check if this is a parallel run."""

    if MPI is None:
        return False
    
    comm = MPI.COMM_WORLD

    return comm.Get_size()>1

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
    if pnt[1]<bound[3]:
        return False
    if pnt[2]<bound[4]:
        return False
    if pnt[2]<bound[5]:
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

def distribute_particles(particle_list,bounds):
    """ Handle exchanging particles across multiple processors """
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if not is_parallel():
        return set(particle_list)

    all_bounds = np.empty([size, 6], dtype=float)

    comm.Allgather(bounds, all_bounds)

    data=[]

    for i in range(size):
        if i==rank: continue
        data.append([par for par in particle_list 
                     if point_in_bound(par.pos, all_bounds[i])])

    data = comm.alltoall(data)

    output=set()
    for i in range(size):
        set.add(data[i])

    return set

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
            out[i] = cell_owned(block, ele)
        else:
            out[i] = False
    
