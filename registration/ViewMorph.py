# -*- coding: utf-8 -*-
import numpy as np
from enthought.mayavi.mlab import *

NC = np.ma.load('TempElasNodes_Iter55_TimeThuOct21')
Tri = np.ma.load('BSkullFTtri')

figure(fgcolor=(0,0,0),bgcolor=(1,1,1))
triangular_mesh(NC[:,0],NC[:,1],NC[:,2],Tri,color=(0.8,0.8,0.6),opacity=1)