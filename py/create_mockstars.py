'''select DM particles with the same E,L distribution as stars'''
from dynio import *
import sys

halo=sys.argv[1]
seed=100
F=Tracer_t()
lib.mock_stars(halo, seed, ctypes.byref(F))
print 'saving...'
lib.save_mockstars(halo, ctypes.byref(F))
