#
#  Copyright (C) 2015-2018 Huw Jenkins & University of York
#
import os

from fragon import version
__version__ = version.__version__

if 'CCP4' not in os.environ:
  raise RuntimeError('could not find a valid CCP4 installation!')
try:
  import clipper
except ImportError:
  raise RuntimeError('Clipper Python module not found!')
