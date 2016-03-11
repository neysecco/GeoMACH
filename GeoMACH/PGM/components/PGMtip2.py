"""
GeoMACH cone class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components.PGMinterpolant import PGMinterpolant
from GeoMACH.PGM.core.PGMface import PGMface


class PGMtip2(PGMinterpolant):
    """ Cone component """

    def __init__(self, config, comp, side, weight=0.1):
        """
        Parameters
        ----------
        config : ``PGMconfiguration``
           Pointer to the configuration containing this component
        comp : ``str``
           Name of the wing component
        side : ``str``
           The side of the wing being closed off:
           'right' for the v=0 side,
           'left' for the v=1 side
        weight : ``float``
           The weight given to the tangent vectors
           of the tip of the wing when interpolating
        """
        super(PGMtip2, self).__init__()

        self._comp = config.comps[comp]
        self._side = side
        self._weight = weight

        nx = self._comp.faces['upp']._num_surf['u']
        ny = self._comp.faces['le']._num_surf['u']

        self._num_surf['x'] = nx
        self._num_surf['y'] = ny

        self.faces[''] = PGMface(nx, ny)

    def set_diff(self):
        # Get reference to the tip face
        face = self.faces['']

        # First set everything to C1 continuity
        face.set_diff_surf(True)

        # Then check if we need to turn off continuity at the trailing edge
        if not self._comp._round_te:
            #face.set_diff_surf(False, ind_i=-1, ind_u=2)
            #face.set_diff_corner(False, ind_i=-1, ind_j=0)
            #face.set_diff_corner(False, ind_i=-1, ind_j=-1)
            pass

    def compute(self, name):
        face = self.faces['']
        #num_u = self._comp.faces['upp']._num_cp_total['u']
        #num_v = self._comp.faces['le']._num_cp_total['u']
        num_u = face._num_cp_total['u']
        num_v = face._num_cp_total['v']

        if name == 'cp_prim': #If we are at the cp_prim step...
            return super(PGMtip2, self).compute(name) #Call the function that sets up the normal properties
        elif name == 'cp_bez': #If we are at the cp_bez step...
            rgt = self._comp.faces['upp']
            top = self._comp.faces['le']
            lft = self._comp.faces['low']
            bot = self._comp.faces['te']
            if self._side == 'right':
                W = rgt.vec_inds['cp_prim'][::-1,:2,:]
                N = top.vec_inds['cp_prim'][:,:2,:]
                E = lft.vec_inds['cp_prim'][:,:2,:]
                S = bot.vec_inds['cp_prim'][::-1,:2,:]
            elif self._side == 'left':
                E = rgt.vec_inds['cp_prim'][::-1,-1:-3:-1,:]
                N = top.vec_inds['cp_prim'][::-1,-1:-3:-1,:]
                W = lft.vec_inds['cp_prim'][:,-1:-3:-1,:]
                S = bot.vec_inds['cp_prim'][:,-1:-3:-1,:]
            nD = 0
            nD += 3 * 2 * num_v + 3 * 2 * (num_u - 2)
            nD += 3 * 8 * (num_v - 3 + num_u - 3)
            nD += 3 * 8

            print E.shape
            print N.shape
            print W.shape
            print S.shape

            Da, Di, Dj = PGMlib.computeconewireframe(nD, num_u, num_v, 1.0, self._weight, W, E, N, S, face.vec_inds['cp_bez'])
            Das, Dis, Djs = super(PGMtip2, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            return Das + [Da], Dis + [Di], Djs + [Dj]
        elif name == 'cp_coons': #If we are at the cp_coons step...
            nD = 3 * 8 * 4 * ((num_u-1)/2 -1) * ((num_v-1)/2 -1)
            Da, Di, Dj = PGMlib.computeconecoons(nD, num_u, num_v, face.vec_inds['cp_coons'])
            Das, Dis, Djs = super(PGMtip2, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            return Das + [Da], Dis + [Di], Djs + [Dj]
            
