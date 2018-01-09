""" Debugging routines for the particle_model code."""

import logging

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

logging.basicConfig(format='%(module)s:%(funcName)s:%(lineno)d - %(message)s')
logger = logging.getLogger('particle_model')
