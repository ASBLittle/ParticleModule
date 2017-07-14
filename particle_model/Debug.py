try:
    import builtins
except:
    import __builtin__ as builtins

try:
    profile = builtins.__dict__['profile']
except:
    def profile(obj):
        return obj
