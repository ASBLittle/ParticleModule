try:
    import builtins
except:
    import __builtin__ as builtins

try:
    profile = builtins.__dict__['profile']
except:
    def profile(obj):
        return obj

import logging
logging.basicConfig(format='%(module)s:%(funcName)s:%(lineno)d - %(message)s')
logger = logging.getLogger('particle_model')


