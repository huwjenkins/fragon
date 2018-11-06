from __future__ import absolute_import, division, print_function
import os

from fragon import version
__version__ = version.__version__

if 'CCP4' not in os.environ:
  raise RuntimeError('could not find a valid CCP4 installation!')
