"""
GeoMACH junction class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components.PGMinterpolant import PGMinterpolant
from GeoMACH.PGM.core.PGMface import PGMface


class PGMjunction2(PGMinterpolant):
    """ Junction component """

    def __init__(self, config, fcomp, face, 
                 vdir, loc, mcomp, side,
                 fweight=0, mweight=0):
        """
        Parameters
        ----------
        config : ``PGMconfiguration``
           Pointer to the configuration containing this component
        fcomp : ``str``
           Name of the 'female' component
           (the one to which the wing is attached)
        face : ``str``
           Name of the face on the female component
           to which the wing is attached
        vdir : ``str``
           Direction of the v axis (of the face of the female component)
           is pointing: 'E', 'N', 'W', or 'S'
        loc : ``list`` of 2 ``int``s
           The 'i' and 'j' indices of the northwest-most face
           of the junction when viewed
           with the wing's upper surface facing up
        mcomp : ``str``
           Name of the 'male' component
           (the one being attached to the female component)
        side : ``str``
           The side of wing being attached:
           'right' for the v=0 side,
           'left' for the v=1 side
        fweight : ``float``
           The weight applied to the tangent vectors
           of the female component when interpolating
        mweight : ``float``
           The weight given to the tangent vectors
           of the male component when interpolating
        """
        super(PGMjunction2, self).__init__()

        self._fcomp = config.comps[fcomp].faces[face]
        self._loc = {'u': loc[0], 'v': loc[1]}
        self._mcomp = config.comps[mcomp]
        self._side = side

        #Check if the user gave a single weight for all the edges or if he specified one weight for each edge
	self._fweight = numpy.zeros(6)
	self._mweight = numpy.zeros(6)
	self._fweight[:] = fweight
	self._mweight[:] = mweight

        # Get the number of surfaces in two consecutive faces of the male component
        keys = self._mcomp.faces.keys()
        self._num_u = self._mcomp.faces[keys[1]]._num_surf['u']
        self._num_v = self._mcomp.faces[keys[0]]._num_surf['u']
        # Now set the number of faces of the junction
        #
        #
        #    +---> v
        #    | +--------+--------------+---------------+-----------+  ---
        #    | |        |              |               |           |   |
        #    V |        |              |               |           |   |
        #    u |        |              |               |           |   |
        #      +--------+--------------+---------------+-----------+   |
        #      |        |         first face           |           |   |
        #      |        |                              |           |   |
        #      |        |                              |           |   |
        #      +--------+ 2nd        male              +-----------+   |
        #      |        | face     component           |           |   |  # of male faces + 2
        #      |        |           outline            |           |   |
        #      +--------+                              +-----------+   |
        #      |        |                              |           |   |
        #      |        |                              |           |   |
        #      |        |                              |           |   |
        #      +--------+--------------+---------------+-----------+   |
        #      |        |              |               |           |   |
        #      |        |              |               |           |   |
        #      |        |              |               |           |   |
        #      |        |              |               |           |   |
        #      +--------+--------------+---------------+-----------+  ---
        #
        #      |---------------------------------------------------|
        #                       # of male faces + 2

        self.faces[''] = PGMface(2 + self._num_u, 2 + self._num_v)
        # Blank interior surfaces (which would be inside the male component)
        self.faces['']._surf_indices[1:-1, 1:-1] = -1

        if vdir == 'E':
            self._rotate = lambda P: P
            self._flip = lambda nu, nv: [nu,nv]
        elif vdir == 'N':
            self._rotate = lambda P: numpy.swapaxes(P,0,1)[::-1,:]
            self._flip = lambda nu, nv: [nv[::-1],nu]
        elif vdir == 'W':
            self._rotate = lambda P: P[::-1,::-1]
            self._flip = lambda nu, nv: [nu[::-1],nv[::-1]]
        elif vdir == 'S':
            self._rotate = lambda P: numpy.swapaxes(P,0,1)[:,::-1]
            self._flip = lambda nu, nv: [nv,nu[::-1]]

    def set_diff(self):

        # Here we assign continuity to the edges

        face = self.faces['']
        face.set_diff_surf(True)
        
        '''
        for ind_j in xrange(1, 1 + self._num_surf_wing):
            face.set_diff_surf(False, ind_i=0, ind_j=ind_j, ind_u=2)
            face.set_diff_surf(False, ind_i=-1, ind_j=ind_j, ind_u=0)

        face.set_diff_surf(False, ind_i=0, ind_j=0, ind_u=2, ind_v=2)
        face.set_diff_surf(False, ind_i=-1, ind_j=0, ind_u=0, ind_v=2)
        face.set_diff_surf(False, ind_i=0, ind_j=-1, ind_u=2, ind_v=0)
        face.set_diff_surf(False, ind_i=-1, ind_j=-1, ind_u=0, ind_v=0)
	
        for ind_j in xrange(2 + self._num_surf_wing):
            face.set_diff_surf(False, ind_i=0, ind_j=ind_j, ind_u=0)
            face.set_diff_surf(False, ind_i=-1, ind_j=ind_j, ind_u=2)
	for ind_i in [0,-1]:
            face.set_diff_surf(False, ind_i=ind_i, ind_j=0, ind_v=0)
            face.set_diff_surf(False, ind_i=ind_i, ind_j=-1, ind_v=2)
        '''
        loc = self._loc

        fu = self._fcomp._num_cp_list['u']
        fv = self._fcomp._num_cp_list['v']
        fu,fv = self._flip(fu,fv)
        fu1 = sum(fu[:loc['u']])
        fu2 = sum(fu[:loc['u']+2+self._num_u])
        fv1 = sum(fv[:loc['v']])
        fv2 = sum(fv[:loc['v']+2+self._num_v])
        fFace_inds = self._rotate(self._fcomp.vec_inds['cp_bez'])
        fFace_inds = fFace_inds[fu1:fu2+1,fv1:fv2+1]

    def set_hidden_surfaces(self):
        '''
        loc = self._loc
        fInds = self._rotate(self._fcomp._surf_indices)
        for ind_v in xrange(2 + self._num_surf_wing):
            for ind_u in xrange(2):
                isurf = fInds[ind_u + loc['u'], ind_v + loc['v']]
                self._bse.hidden[isurf] = True
        '''

    def compute(self, name):
        loc = self._loc

        # Get lists that contains the number of control points along each
        # surface of the female component
        fu = self._fcomp._num_cp_list['u']
        fv = self._fcomp._num_cp_list['v']
        # Flip according to the junction orientation
        fu,fv = self._flip(fu,fv)
        # Now we can find how many control points are in the female component
        # up to the junction location
        fu1 = sum(fu[:loc['u']]) # Total number of control points at the beginning of the junction
        fu2 = sum(fu[:loc['u']+2+self._num_u]) # Total number of control points at the end of the junction
        fv1 = sum(fv[:loc['v']]) # Total number of control points at the beginning of the junction
        fv2 = sum(fv[:loc['v']+2+self._num_v]) # Total number of control points at the end of the junction
        fFace_inds = self._rotate(self._fcomp.vec_inds['cp_bez'])
        fFace_inds = fFace_inds[fu1:fu2+1,fv1:fv2+1]

        num_u = [fv[loc['u']] + 1, 
                 sum(fv[loc['u']+1:loc['u'] + 1 + self._num_u]) + 1, 
                 fv[loc['u'] + 1 + self._num_u] + 1]
        num_v = [fv[loc['v']] + 1, 
                 sum(fv[loc['v']+1:loc['v'] + 1 + self._num_v]) + 1, 
                 fv[loc['v'] + 1 + self._num_v] + 1]
        nu0 = sum(num_u) - 2
        nv0 = sum(num_v) - 2

        print 'num_u: ',num_u
        print 'num_v: ',num_v


        if name == 'cp_prim': #If we are at the cp_prim step...
            return super(PGMjunction2, self).compute(name) #Call the function that sets up the normal properties
        elif name == 'cp_bez':

            # Now we need to get the last two layers of control points that define
            # the male component outline
            #
            #
            #    +---> v
            #    | +--------+--------------+---------------+-----------+  ---
            #    | |        |              |               |           |   |
            #    V |        |              |               |           |   |
            #    u |        |              |               |           |   |
            #      +--------+--------------+---------------+-----------+   |
            #      |        |   |-->       N               |           |   |
            #      |        |                              |           |   |
            #      |        | ^                          ^ |           |   |
            #      +--------+ |          male            | +-----------+   |
            #      |        |---       component        ---|           |   |  # of male faces + 2
            #      |        |W          outline          E |           |   |
            #      +--------+                              +-----------+   |
            #      |        |                              |           |   |
            #      |        |                              |           |   |
            #      |        |   |-->       S               |           |   |
            #      +--------+--------------+---------------+-----------+   |
            #      |        |              |               |           |   |
            #      |        |              |               |           |   |
            #      |        |              |               |           |   |
            #      |        |              |               |           |   |
            #      +--------+--------------+---------------+-----------+  ---
            #
            #      |---------------------------------------------------|
            #                       # of male faces + 2
        
            # Get reference to the male component
            mcomp = self._mcomp

            # Get the names of the faces of the male component
            keys = mcomp.faces.keys()

            # Now gather the control points according to the figure shown above
            if self._side == 'right':
                N = mcomp.faces[keys[0]].vec_inds['cp_prim'][::-1,:2,:]
                W = mcomp.faces[keys[1]].vec_inds['cp_prim'][::-1,:2,:]
                S = mcomp.faces[keys[2]].vec_inds['cp_prim'][:,:2,:]
                E = mcomp.faces[keys[3]].vec_inds['cp_prim'][:,:2,:]
            elif self._side == 'left':
                N = mcomp.faces[keys[0]].vec_inds['cp_prim'][:,-1:-3:-1,:]
                W = mcomp.faces[keys[3]].vec_inds['cp_prim'][:,-1:-3:-1,:]
                S = mcomp.faces[keys[2]].vec_inds['cp_prim'][::-1,-1:-3:-1,:]
                E = mcomp.faces[keys[1]].vec_inds['cp_prim'][::-1,-1:-3:-1,:]

            print 'N_shape: ',N.shape
            print 'W_shape: ',W.shape
            print 'E_shape: ',E.shape
            print 'S_shape: ',S.shape

            fInds = -numpy.ones((nu0, nv0, 3), dtype=int, order='F')
            fInds[:num_u[0],:] = fFace_inds[:num_u[0],:]
            fInds[-num_u[2]:,:] = fFace_inds[-num_u[2]:,:]

            nD = 3 * 2 * nv0 + 3 * 2 * (nu0 - 2)
            nD += 3 * 2
            if num_u[1] != 1:
                nD += 3 * 2
            nD += 3 * 2 * (num_v[1] - 2)
            if num_u[1] != 1:
                nD += 3 * 2 * (num_u[1] - 2)
            nD += 3 * 4 * 2 * (num_u[0] - 2)
            nD += 3 * 4 * 2 * (num_u[2] - 2)
            nD += 3 * 4 * (num_v[0] - 2)
            nD += 3 * 4 * (num_v[2] - 2)
            if num_u[1] != 1:
                nD += 3 * 4 * (num_v[0] - 2)
                nD += 3 * 4 * (num_v[2] - 2)

            Da, Di, Dj = PGMlib.computejunctionwireframe(nD, nu0, nv0, num_u[0], num_u[1], num_u[2], num_v[0], num_v[1], num_v[2], self._fweight, self._mweight, W, E, N, S, fInds, self.faces[''].vec_inds['cp_bez'])
            Da = Da * (-1 != Dj)
            Dj = numpy.maximum(0, Dj)

            Das, Dis, Djs = super(PGMjunction2, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            return Das + [Da], Dis + [Di], Djs + [Dj]
        elif name == 'cp_coons': #If we are at the cp_coons step...
            nD = 0
            for i in range(3):
                for j in range(3):
                    if (num_u[1] != 1) or (i !=1):
                        nD += 3 * 8 * (num_u[i]-2) * (num_v[j]-2)
        
            Das, Dis, Djs = super(PGMjunction2, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            Da, Di, Dj = PGMlib.computejunctioncoons(nD, nu0, nv0, num_u[0], num_u[1], num_u[2], num_v[0], num_v[1], num_v[2], self.faces[''].vec_inds['cp_coons'])
            return Das + [Da], Dis + [Di], Djs + [Dj]
