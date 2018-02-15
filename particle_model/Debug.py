""" Debugging routines for the particle_model code."""

import logging
import atexit

try:
    import builtins
except ImportError:
    import __builtin__ as builtins

try:
    profile = builtins.__dict__['profile']
except KeyError:
    def profile(obj):
        """Do nothing profiler."""
        return obj
    builtins.__dict__['profile'] = profile

def make_line_profiler(*args, **kwargs):
    """Register and return a line profiler."""
    import line_profiler
    lprof = line_profiler.LineProfiler()
    atexit.register(profile.print_stats, *args, **kwargs)
    return lprof


logging.basicConfig(format='%(module)s:%(funcName)s:%(lineno)d - %(message)s')
logger = logging.getLogger('particle_model')
