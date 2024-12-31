import sys

sys.modules['problem'] = None
sys.modules['optimizer'] = None

from .problem import *
from .optimizer import *