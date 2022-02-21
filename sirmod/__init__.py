'''
Initialize HPVmod by importing all the modules
'''

# import requirements
#
# # Import the version and print the license unless verbosity is disabled, via e.g. os.environ['COVASIM_VERBOSE'] = 0
# from .version import __version__, __versiondate__, __license__
# if settings.options.verbose:
#     print(__license__)

# Import the actual model
from .parameters    import * # Depends on
from .utils         import * # Depends on defaults
from .model         import * # Depends on almost everything


