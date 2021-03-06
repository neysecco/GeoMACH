"""
GeoMACH wing class
This class has 4 faces (upper skin, leading edge, lower skin, and trailing edge)
John Hwang, July 2014
Ney Secco, March 2016
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse, scipy.sparse.linalg
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM.components.PGMprimitive import PGMprimitive
from GeoMACH.PGM.core.PGMface import PGMface


class PGMwing2(PGMprimitive):
    """ Wing component """

    def __init__(self, num_x=1, num_z=1,
                 left_end='open', right_end='open', round_te=False):
        """
        Parameters
        ----------
        num_x : ``integer``
           Number of surfaces in the chord-wise direction
        num_z : ``integer``
           Number of surfaces in the span-wise direction
        left_end, right_end : ``string``
           ``open`` if open,
           ``tip`` if attached to tip component
           ``junction`` if attached to a junction component
        """
        super(PGMwing2, self).__init__()

        self._num_surf['x'] = num_x
        self._num_surf['y'] = 1
        self._num_surf['z'] = num_z

        self._left_end = left_end
        self._right_end = right_end
        self._round_te = round_te

        self._ax1 = 3
        self._ax2 = 2

        # Define faces
        # We will assume that the LE and TE have only one surface
        self.faces['upp'] = PGMface(num_x, num_z)
        self.faces['le'] = PGMface(1, num_z)
        self.faces['low'] = PGMface(num_x, num_z)
        self.faces['te'] = PGMface(1, num_z)

    def assemble_sizes(self, bse):
        super(PGMwing2, self).assemble_sizes(bse)

        if self._bse is None:
            # Here we assign a wing with rectangular cross section of thickness 10%

            # The shapes array are of the form
            # self._shapes[<face>] = [num_u, num_v, 3]
            # The first two indices indicates the parametric position of the control point
            # while the last index represents the physical coordinate of the control point (x,y, or z)

            # Get references to the control point coordinates
            upp = self._shapes['upp']
            le = self._shapes['le']
            low = self._shapes['low']
            te = self._shapes['te']
            # Initialize all coordinates to zero
            upp[:, :, :] = 0.0
            low[:, :, :] = 0.0
            le[:, :, :] = 0.0
            te[:, :, :] = 0.0

            # We will set coordinates so that we have a rectangular cross section of thickness 10%

            # UPPER AND LOWER SKINS
            # We will add linearly-spaced control points in the x (chord) direction
            num = upp.shape[0]
            for ind in range(num):
                upp[ind, :, 0] = 1 - ind / (num-1)
            num = low.shape[0]
            for ind in range(num):
                low[ind, :, 0] = ind / (num-1)
            # all points in the upper skin constrained to y=0.05
            upp[:, :, 1] = 0.05
            # all points in the lower skin constrained to y=0.05
            low[:, :, 1] = -0.05

            # LEADING AND TRAILING EDGES
            # We will add linearly-spaced control points in the y (vertical) direction
            num = le.shape[0]
            for ind in range(num):
                le[ind, :, 1] = 0.05 - (ind / (num-1))*0.1
            num = te.shape[0]
            for ind in range(num):
                te[ind, :, 1] = -0.05 + (ind / (num-1))*0.1
            # all points in the leading edge constrained to x=0.0
            le[:, :, 0] = 0.0
            # all points in the trailing edge constrained to x=1.0
            te[:, :, 0] = 1.0

            # COMMON STEPS
            # We will also add linearly-spaced control points in the z (span) direction
            num_span = upp.shape[1]
            for ind in range(num_span):
                z_station = ind / (num_span-1)
                upp[:, ind, 2] = z_station
                low[:, ind, 2] = z_station
                le[:, ind, 2] = z_station
                te[:, ind, 2] = z_station

        else:
            self.set_airfoil()

    def compute(self, name):
        '''
        This method will run the "compute" method of the PGMprimitive class.
        It also will update the trailing edge control points if necessary.
        '''

        # Run the "compute" method from the PGMprimitive class
        vals_list, rows_list, cols_list = super(PGMwing2, self).compute(name)

        '''
        # Now we update the trailing edge points
        if name == 'cp_prim':
            face = self.faces['te']

            # Import bezier curve function
            from GeoMACH.PGM.PGMlib import sparsebezier
            # Get number of control points in the trailing edge
            num_u = self.faces['te']._num_cp_total['u']
            for i in range(num_u):
                # Use a Bezier curve to get the weighting coefficients
                C = sparsebezier(i/(num_u-1), -0.5, 0.5)
                # Apply weights to the jacobian for each coordiante (x,y, and z)
                for k in range(3):
                    vals = C[:]
                    

                # Apply weights to find coordinates of the trailing edge points
                #self._shapes['te'][i,-1,:] = C[0]*self._shapes['te'][0,-1,:] + \
                #                             C[1]*self._shapes['te'][0,-2,:] + \
                #                             C[2]*self._shapes['te'][-1,-1,:] + \
                #                             C[3]*self._shapes['te'][-1,-2,:]
        '''

        #RETURNS
        return vals_list, rows_list, cols_list

    def set_diff(self):
        # This method set differentiability options for the control points
        # We need to use C0 continuity if we want sharp corners

        # Get the number of faces
        num_faces = len(self.faces)

        # First set C1 everywhere
        for ind in xrange(num_faces):
            face = self.faces.values()[ind]
            face.set_diff_surf(True)

        # Now we turn off the continuities we want

        # UPPER SKIN
        # u - TE to LE
        # v - right to left (from pilot's view)
        # Get reference to the upper surface
        face = self.faces['upp']
        #C0 left edge if we have and open end
        if self._left_end is 'open':
            face.set_diff_surf(False, ind_j=-1, ind_v=2)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v1', ind_j=-1)
        #C0 right edge if we have an open end
        if self._right_end is 'open':
            face.set_diff_surf(False, ind_j=0, ind_v=0)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v0', ind_j=0)
        #C0 trailing edge if we have it sharp
        if not self._round_te:
            face.set_diff_surf(False, ind_i=0, ind_u=0)
            # Keep continuity along the edge only
            face.set_diff_edge(True, 'u0', ind_i=0)
            # Turn off continuity at TE corners where we have open ends
            if self._left_end is 'open':
                face.set_diff_corner(False, ind_i=0, ind_j=-1)
            if self._right_end is 'open':
                face.set_diff_corner(False, ind_i=0, ind_j=0)

        # LEADING EDGE
        # u - up to down
        # v - right to left (from pilot's view)
        # Get reference to the leading edge
        face = self.faces['le']
        #C0 left edge if we have and open end
        if self._left_end is 'open':
            face.set_diff_surf(False, ind_j=-1, ind_v=2)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v1', ind_j=-1)
        #C0 right edge if we have an open end
        if self._right_end is 'open':
            face.set_diff_surf(False, ind_j=0, ind_v=0)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v0', ind_j=0)

        # LOWER SKIN
        # u - LE to TE
        # v - right to left (from pilot's view)
        # Get reference to the lower surface
        face = self.faces['low']
        #C0 left edge if we have and open end
        if self._left_end is 'open':
            face.set_diff_surf(False, ind_j=-1, ind_v=2)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v1', ind_j=-1)
        #C0 right edge if we have an open end
        if self._right_end is 'open':
            face.set_diff_surf(False, ind_j=0, ind_v=0)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v0', ind_j=0)
        #C0 trailing edge if we have it sharp
        if not self._round_te:
            face.set_diff_surf(False, ind_i=-1, ind_u=2)
            # Keep continuity along the edge only
            face.set_diff_edge(True, 'u1', ind_i=-1)
            # Turn off continuity at TE corners where we have open ends
            if self._left_end is 'open':
                face.set_diff_corner(False, ind_i=-1, ind_j=-1)
            if self._right_end is 'open':
                face.set_diff_corner(False, ind_i=-1, ind_j=0)

        # TRAILING EDGE
        # u - down to up
        # v - right to left (from pilot's view)
        # Get reference to the trailing edge
        face = self.faces['te']
        #C0 left edge if we have and open end
        if self._left_end is 'open':
            face.set_diff_surf(False, ind_j=-1, ind_v=2)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v1', ind_j=-1)
        #C0 right edge if we have an open end
        if self._right_end is 'open':
            face.set_diff_surf(False, ind_j=0, ind_v=0)
            # Still leave edge-edge continuity
            face.set_diff_edge(True, 'v0', ind_j=0)
        #C0 trailing edge if we have it sharp
        if not self._round_te:
            # C0 between lower skin and TE
            face.set_diff_surf(False, ind_i=0, ind_u=0)
            # C0 between upper skin and TE
            face.set_diff_surf(False, ind_i=-1, ind_u=2)
            # But still keep edge-to-edge continuity along the edges between TE and skins
            face.set_diff_edge(True, 'u0', ind_i=0)
            face.set_diff_edge(True, 'u1', ind_i=-1)
            # Turn on discontinuity at corners where we have open ends
            if self._left_end is 'open':
                face.set_diff_corner(False, ind_i=0, ind_j=-1)
                face.set_diff_corner(False, ind_i=-1, ind_j=-1)
            if self._right_end is 'open':
                face.set_diff_corner(False, ind_i=0, ind_j=0)
                face.set_diff_corner(False, ind_i=-1, ind_j=0)

    '''
    def set_airfoil(self, filename=['naca0012','naca0012'], pos_v =[0,1], order_v = 2, LE_pos=0.05, bunch_LE=1.0, bunch_TE=1.0):
        from GeoMACH.PGM.PGMlib import computebspline
        from numpy import prod

        if type(filename) == str:
            filename = [filename, filename]

        numSections = len(filename)

        for name in ['upp', 'le', 'low', 'te']: #Loop for each face
            numCpsPerChordSurf = self.faces[name]._num_cp_list['u'] #Get the number of chord-wise control points
            numPtsPerChordSurf = self.faces[name]._num_pt_list['u'] #Get the number of chord-wise points
            numCpsPerSpanSurf = self.faces[name]._num_cp_list['v'] #Get the number of span-wise control points
            numPtsPerChord = sum(numPtsPerChordSurf) + 1

            #First we need to create an array that contains all the normalized control points that belongs to the specified sections. These control points will be interpolated from the airfoil coordinate points.
            Q_total = numpy.zeros((sum(numCpsPerChordSurf)+1,numSections,3))
            for index in range(len(filename)): #For each specified airfoil
                current_filename = filename[index] #Get the airfoil name
                if current_filename[:4]=='naca' or current_filename[:4]=='NACA':
                    airfoils = self._get_airfoil_naca(current_filename[4:],LE_pos)
                else:
                    airfoils = self._get_airfoil_file(current_filename, LE_pos)

                P = self._get_P(numPtsPerChord, airfoils, name, bunch_LE, bunch_TE)
                P[:, 0] /= numpy.max(P[:, 0])
                P[:, 1] -= numpy.linspace(P[0,1], P[-1,1], P.shape[0])

                
                if name == 'upp':
                    sign = 1.0
                elif name == 'low':
                    sign = -1.0
                t, p, x = blunt_thk, blunt_pos, P[:, 0]
                P[:, 1] += sign * t * (-2*(x/p)**3 + 3*(x/p)**2) * (x < p)
                P[:, 1] += sign * t * numpy.sqrt(1 - (numpy.maximum(x,2*p-1) - p)**2/(1-p)**2) * (x >= p)
                

                Q = self._get_Q(numCpsPerChordSurf, numPtsPerChordSurf, P) #These are the normalized control points for this section
                Q_total[:,index,:] = Q #Store them into the full matrix

            #Interpolate using B-splines
            #ATTENTION: in this part of code "control points" refers to the control points of the sections where the airfoil is specified, while "points" refers to the control points of the sections where the airfoil is not specified. The transformation will be the same for every chord-wise points, so we just need to find the span-wise interpolation. That's why we use only 2 control points in u-direction for this.
            num_pt = {'u':0, 'v':0}
            order = {'u':0, 'v':0}
            num_cp = {'u':0, 'v':0}
            pos = {'u':0, 'v':0}

            num_pt['u'] = 2
            num_pt['v'] = sum(numCpsPerSpanSurf)+1
            order['u'] = 2
            order['v'] = order_v
            num_cp['u'] = 2
            num_cp['v'] = len(pos_v)
            nnz = num_pt['u'] * num_pt['v'] * order['u'] * order['v']
            cp_indices = numpy.array(range(num_cp['u']*num_cp['v'])).reshape((num_cp['u'],num_cp['v']))
            pt_indices = numpy.array(range(num_pt['u']*num_pt['v'])).reshape((num_pt['u'],num_pt['v']))
            pos['u'] = numpy.array([0,1])
            pos['v'] = pos_v

            vals, rows, cols \
                = computebspline(nnz, 
                                 order['u'], order['v'],
                                 num_cp['u'], num_cp['v'],
                                 num_pt['u'], num_pt['v'],
                                 cp_indices, pt_indices,
                                 pos['u'], pos['v'])

            #Reshape transformation matrix
            dQbar_dQ = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(num_pt['u']*num_pt['v'], num_cp['u']*num_cp['v']))

            #TRANFORMATION
            #Preallocate the interpolated control points
            Qbar = numpy.zeros((sum(numCpsPerChordSurf)+1,sum(numCpsPerSpanSurf)+1,3))
            vec_a = numpy.zeros((2, num_cp['v'], 3))
            vec_b = numpy.zeros((2, num_pt['v'], 3))
            for chord_index in range(sum(numCpsPerChordSurf)+1):
                vec_a[0, :, :] = Q_total[chord_index,:,:]
                outvec = dQbar_dQ.dot(vec_a.reshape((2 * num_cp['v'], 3)))
                vec_b[:, :, :] = outvec.reshape((2, num_pt['v'], 3))
                Qbar[chord_index,:,:] = vec_b[0, :, :]

            #Save the interpolated control points
            for j in range(self.faces[name]._num_cp_total['v']):
                self._shapes[name][:,j,:] = Qbar[:,j,:]
                self._shapes[name][:,j,2] = 0.0
    '''

    def set_airfoil(self, filename='naca0012', blunt_thk=0.0, blunt_pos=0.95, LE_pos=0.05, bunch_LE=1.0, bunch_TE=1.0):
        if filename[:4]=='naca' or filename[:4]=='NACA':
            airfoils = self._get_airfoil_naca(filename[4:],LE_pos)
        else:
            airfoils = self._get_airfoil_file(filename,LE_pos)

        import matplotlib.pyplot as plt
        fig = plt.figure()

        for name in ['upp', 'le', 'low', 'te']:
            ms = self.faces[name]._num_cp_list['u']
            ns = self.faces[name]._num_pt_list['u']
            nP = sum(ns) + 1

            P = self._get_P(nP, airfoils, name, bunch_LE, bunch_TE)
            #P[:, 1] -= numpy.linspace(P[0,1], P[-1,1], P.shape[0])

            Q = self._get_Q(ms, ns, P)
            num_v = self.faces[name]._num_cp_total['v']
            for j in range(num_v):
                if name in ['upp','low']:
                    # Make sure the last point is at x=1.0
                    Q[:, 0] /= numpy.max(Q[:, 0])
                elif name == 'te':
                    # Set TE to x=1.0
                    Q[:, 0] = Q[0, 0]
                self._shapes[name][:,j,:] = Q[:,:]
                self._shapes[name][:,j,2] = j/(num_v-1)
            
            '''
            plt.plot(self._shapes[name][:,0,0],self._shapes[name][:,0,1])
            plt.show()
            '''

        '''
        # If one of the sides is closed with a wingtip, we need to do a rounded trailing edge
        # This will be shaped as a B-spline using the upper and lower sides of the trailing edge
        if not self._left_closed:
            # Import bezier curve function
            from GeoMACH.PGM.PGMlib import sparsebezier
            # Get number of control points in the trailing edge
            num_u = self.faces['te']._num_cp_total['u']
            for i in range(num_u):
                # Use a Bezier curve to get the weighting coefficients
                C = sparsebezier(i/(num_u-1), -0.5, 0.5)
                # Apply weights to find coordinates of the trailing edge points
                self._shapes['te'][i,-1,:] = C[0]*self._shapes['te'][0,-1,:] + \
                                             C[1]*self._shapes['te'][0,-2,:] + \
                                             C[2]*self._shapes['te'][-1,-1,:] + \
                                             C[3]*self._shapes['te'][-1,-2,:]
            self._shapes['te'][:,-1,2] = self._shapes['te'][:,-1,2] + (1.0 - max(self._shapes['te'][:,-1,2]))
            print self._shapes['te'][:,-1,2]
                
        if not self._right_closed:
            # Import bezier curve function
            from GeoMACH.PGM.PGMlib import sparsebezier
            # Get number of control points in the trailing edge
            num_u = self.faces['te']._num_cp_total['u']
            for i in range(num_u):
                # Use a Bezier curve to get the weighting coefficients
                C = sparsebezier(i/(num_u-1), -0.5, 0.5)
                # Apply weights to find coordinates of the trailing edge points
                self._shapes['te'][i,0,:] = C[0]*self._shapes['te'][0,0,:] + \
                                             C[1]*self._shapes['te'][0,1,:] + \
                                             C[2]*self._shapes['te'][-1,0,:] + \
                                             C[3]*self._shapes['te'][-1,1,:]
        '''

    def _get_Q(self, ms, ns, P0):
        nsurf = ns.shape[0]

        Ps = []
        for i in xrange(nsurf):
            P = numpy.zeros((ns[i]+1,4,3), order='F')
            for j in xrange(4):
                P[:,j,:] = P0[sum(ns[:i]):sum(ns[:i+1])+1,:]
                P[:,j,2] = j
            Ps.append(P)

        bse = BSEmodel(Ps)
        for i in xrange(nsurf):
            bse.set_bspline_option('num_cp', i, 'u', ms[i]+1)
            bse.set_bspline_option('num_pt', i, 'u', ns[i]+1)
            bse.set_bspline_option('num_pt', i, 'v', 4)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array
        jac = bse.jac['d(pt_str)/d(cp_str)']

        for i in xrange(nsurf):
            bse.vec['pt_str'].surfs[i][:,:,:] = Ps[i]
        mtx = jac.T.dot(jac)
        rhs = jac.T.dot(pt)
        for dim in xrange(3):
            cp[:, dim] = scipy.sparse.linalg.cg(mtx, rhs[:, dim])[0]

        Q = numpy.zeros((sum(ms) + 1,3),order='F')
        for i in xrange(nsurf):
            Q[sum(ms[:i]):sum(ms[:i+1])+1,:] = \
                bse.vec['cp_str'].surfs[i][:,0,:]
        return Q

    def _get_P(self, nP, airfoils, name, bunch_LE, bunch_TE):
        def bunch_start(ind, a):
            ind[:] = ind**a
        def bunch_end(ind, a):
            ind[:] = 1 - (1-ind)**a

        airfoil = airfoils[name]

        P = numpy.zeros((airfoil.shape[0],4,3),order='F')
        for j in range(4):
            P[:,j,:2] = airfoil[:,:]
            P[:,j,2] = j
        bse = BSEmodel([P])

        bse.set_bspline_option('num_pt', 0, 'u', 
                               airfoil.shape[0])
        bse.set_bspline_option('num_cp', 0, 'u', 
                               int(airfoil.shape[0]/2))
        bse.set_bspline_option('num_pt', 0, 'v', 4)
        bse.set_bspline_option('num_cp', 0, 'v', 4)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array
        jac = bse.jac['d(pt_str)/d(cp_str)']

        bse.vec['pt_str'].surfs[0][:, :, :] = P
        mtx = jac.T.dot(jac)
        rhs = jac.T.dot(pt)
        for dim in xrange(3):
            cp[:, dim] = scipy.sparse.linalg.cg(mtx, rhs[:, dim])[0]
        fit = numpy.array(cp)

        bse.set_bspline_option('num_pt', 0, 'u', nP)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array

        cp[:, :] = fit[:, :]
        surfs = numpy.zeros(nP).astype(int)
        ind_u = numpy.linspace(0, 1, nP)

        if name == 'upp':
            bunch_start(ind_u, bunch_TE)
            bunch_end(ind_u, bunch_LE)
        elif name == 'low':
            bunch_start(ind_u, bunch_LE)
            bunch_end(ind_u, bunch_TE)
        elif name == 'le':
            bunch_start(ind_u, bunch_TE)
            bunch_end(ind_u, bunch_LE)

        ind_v = numpy.zeros(nP)
        bse.add_jacobian('new', surfs, ind_u, ind_v, 3)
        bse.apply_jacobian('new', 'd(new)/d(cp_str)', 'cp_str')

        return bse.vec['new'].array[:, :]

    def _get_airfoil_naca(self, naca, LE_pos):

        # This function generates NACA coordinates and assigns them
        # to each face.

        # Define number of points
        num = 100 # Total number of points from Leading Edge to Trailing Edge

        # Read airfoil data from the given numbers
        max_cmb = int(naca[0]) / 100.0
        pos_cmb = int(naca[1]) / 10.0
        thickness = int(naca[2:4]) / 100.0

        # Generate cossine-spaced x coordinates
        x = 0.5 * (1 - numpy.cos(numpy.pi * 
                                 numpy.linspace(0, 1, num)))
        # Compute thickness distribution
        y_sym = thickness/0.2 * (0.2969*x**0.5 - 0.1260*x - 
                                 0.3516*x**2 + 0.2843*x**3 - 
                                 0.1015*x**4)
        # Compute camber line
        y_cmb = numpy.zeros(num)
        if pos_cmb != 0:
            y1 = max_cmb * x / pos_cmb**2 * \
                 (2 * pos_cmb - x)
            y2 = max_cmb * (1-x) / (1-pos_cmb)**2 * \
                 (1 + x - 2 * pos_cmb)
            for i in range(num):
                if x[i] < p:
                    y_cmb[i] = y1[i]
                else:
                    y_cmb[i] = y2[i]

        # Find how many points belong to the LE using the given LE_pos.
        # LE_pos should be between 0 and 1. Points with x < LE_pos will be assigned
        # to the LE face, while the others will be assigned to the skins
        
        # Find the index that crosses the LE_pos limit
        LE_limit = next(enum[0] for enum in enumerate(x) if enum[1] > LE_pos)

        # Initialize arrays with coordinates on each surface
        upper = numpy.zeros((num-LE_limit+1, 2), order='F')
        lower = numpy.zeros((num-LE_limit+1, 2), order='F')
        le = numpy.zeros((2*LE_limit, 2), order='F')
        te = numpy.zeros((12, 2), order='F') # We'll add twelve vertical points on the TE

        # Add upper skin coordinates (from TE to LE)
        upper[:,0] = x[:LE_limit-2:-1]
        upper[:,1] = y_sym[:LE_limit-2:-1] + y_cmb[:LE_limit-2:-1]

        # Add lower skin coordinates (from LE to TE)
        lower[:,0] = x[LE_limit-1:]
        lower[:,1] = -y_sym[LE_limit-1:] + y_cmb[LE_limit-1:]

        # Add LE coordinates
        # Upper part
        le[:LE_limit,0] = x[LE_limit-1::-1]
        le[:LE_limit,1] = y_sym[LE_limit-1::-1] + y_cmb[LE_limit-1::-1]
        # Lower part
        le[LE_limit:,0] = x[:LE_limit]
        le[LE_limit:,1] = -y_sym[:LE_limit] + y_cmb[:LE_limit]

        # Add TE coordinates (just a straigth line at x=1)
        te[:,0] = 1.0
        te[:,1] = numpy.linspace(lower[-1,1],upper[0,1],12) # Connection between TE and lower skin

        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.plot(upper[:,0],upper[:,1])
        plt.plot(lower[:,0],lower[:,1])
        plt.plot(le[:,0],le[:,1])
        plt.plot(te[:,0],te[:,1])
        plt.show()
        '''

        return {'upp': upper, 'le': le, 'low': lower, 'te': te}

    def _get_airfoil_file(self, filename, LE_pos):

        path = __import__(__name__).__file__
        index_slash = path[::-1].index('/')
        path = path[:-index_slash]
        data = numpy.genfromtxt(path+'PGM/airfoils/'+filename)

        if data[0,0] > data[1,0]:
            mark = numpy.argmin(data,0)[0]
            upper = data[mark+1::-1,:]
            lower = data[mark:,:]
        else:
            for i in range(data.shape[0]-1):
                if abs(data[i+1,0]-data[i,0]) > 0.8:
                    mark = i
            upper = data[:mark,:]
            lower = data[mark+1:,:]

        # Rescale things so that chord ends at 1.0
        upper = upper/numpy.max(upper[:,0])
        lower = lower/numpy.max(lower[:,0])

        # Now we need to crop the upper and lower skins to give room for the LE
        # First for the upper skin
        # Find the index that crosses the LE_pos limit
        LE_limit = next(enum[0] for enum in enumerate(upper[:,0]) if enum[1] > LE_pos)
        print LE_limit
        # Add upper skin coordinates (from TE to LE)
        upper_crop = upper[:LE_limit-2:-1,:]
        # Add upper LE coordinates
        le_crop = upper[LE_limit-1::-1,:]

        # Now the lower skin
        # Find the index that crosses the LE_pos limit
        LE_limit = next(enum[0] for enum in enumerate(lower[:,0]) if enum[1] > LE_pos)
        # Add lower skin coordinates (from LE to TE)
        lower_crop = lower[LE_limit-1:,:]
        # Add lower LE coordinates
        le_crop = numpy.vstack([le_crop, lower[:LE_limit,:]])

        # Add TE coordinates (just a straigth line at x=1)
        te_crop = numpy.zeros((12, 2), order='F') # We'll add twelve vertical points on the TE
        te_crop[:,0] = 1.0
        te_crop[:,1] = numpy.linspace(lower_crop[-1,1],upper_crop[0,1],12) # Connection between TE and lower skin

        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        #plt.plot(upper[:,0],upper[:,1])
        #plt.plot(lower[:,0],lower[:,1])
        plt.plot(upper_crop[:,0],upper_crop[:,1],label='upp')
        plt.plot(lower_crop[:,0],lower_crop[:,1],label='low')
        plt.plot(le_crop[:,0],le_crop[:,1],label='le')
        plt.plot(te_crop[:,0],te_crop[:,1],label='te')
        plt.legend()
        plt.show()
        '''

        return {'upp': upper_crop, 'le': le_crop, 'low': lower_crop, 'te': te_crop}

    def add_thk_con(self, name, urange, vrange, factor):
	self.funcs[name] = WingThicknessFunction(self, name, urange, vrange, factor)



class WingFunction(object):

    def __init__(self, comp, name, urange, vrange, factor):
	self.bse = comp._bse
	self.comp = comp
	self.name = name
	self.num_u, self.num_v = len(urange), len(vrange)
	self.factor = factor

	num_u, num_v = self.num_u, self.num_v
	self.pt = numpy.zeros(2*num_u*num_v*3)
	self.pt_array = self.pt.reshape((2,num_u,num_v,3), order='F')

	ni, nj = self.num_u, self.num_v
	locations = {'u':numpy.zeros((ni,nj)), 'v':numpy.zeros((ni,nj))}
	for i in range(ni):
	    for j in range(nj):
		locations['u'][i,j] = urange[i]
		locations['v'][i,j] = vrange[j]

	face = self.comp.faces['upp']

        increments = {'upp': {'u': None, 'v': None}, 'low': {'u': None, 'v': None}}
	for f in ['upp', 'low']:
	    for d in ['u', 'v']:
		n = face._num_surf[d]
		increments[f][d] = numpy.zeros(n+1)
		for i in xrange(n+1):
		    increments[f][d][i] = sum(face._num_cp_list[d][:i]) / sum(face._num_cp_list[d])
	    if f=='upp':
		increments[f]['u'][:] = 1 - increments[f]['u'][::-1]

	surf = {'upp':numpy.zeros((ni,nj), dtype=int, order='F'), 'low':numpy.zeros((ni,nj), dtype=int, order='F')}
	locs = {'upp':{'u':numpy.zeros((ni,nj), order='F'), 'v':numpy.zeros((ni,nj), order='F')}, 'low':{'u':numpy.zeros((ni,nj), order='F'), 'v':numpy.zeros((ni,nj), order='F')}}


	loc_face = {'u': None, 'v': None}
	for f in ['upp', 'low']:
	    for i in range(ni):
		for j in range(nj):
		    if f=='upp':
			loc_face['u'] = 1-locations['u'][i,j]
			loc_face['v'] = locations['v'][i,j]
		    else:
			loc_face['u'] = locations['u'][i,j]
			loc_face['v'] = locations['v'][i,j]
		    loc_surf = {'u':-1, 'v':-1}
		    for d in ['u', 'v']:
			n = face._num_surf[d]
			for k in range(n):
			    if increments[f][d][k] <= loc_face[d] <= increments[f][d][k+1]:
				locs[f][d][i,j] = (loc_face[d] - increments[f][d][k]) / \
				    (increments[f][d][k+1] - increments[f][d][k])
				loc_surf[d] = k
		    if loc_surf['u'] == -1 or loc_surf['v'] == -1:
			raise Exception('Invalid thickness constraint locations')
		    ind = 0 if f=='upp' else 1
		    surf[f][i,j] = self.comp.faces.values()[ind]._surf_indices[loc_surf['u'], loc_surf['v']]

	surf_flat = numpy.zeros(2*ni*nj)
	locs_u_flat = numpy.zeros(2*ni*nj)
	locs_v_flat = numpy.zeros(2*ni*nj)

	a = surf['upp'].flatten(order='F')
	b = surf['low'].flatten(order='F')
	for k in range(ni*nj):
	    surf_flat[2*k] = a[k]
	    surf_flat[2*k + 1] = b[k]

	c = locs['upp']['u'].flatten(order='F')
	d = locs['low']['u'].flatten(order='F')
	for k in range(ni*nj):
	    locs_u_flat[2*k] = c[k]
	    locs_u_flat[2*k + 1] = d[k]

	e = locs['upp']['v'].flatten(order='F')
	f = locs['low']['v'].flatten(order='F')
	for k in range(ni*nj):
	    locs_v_flat[2*k] = e[k]
	    locs_v_flat[2*k + 1] = f[k]

	self.bse.add_jacobian('constr_pts', surf_flat, locs_u_flat, locs_v_flat, ndim=3)
	J = self.bse.jac['d(constr_pts)/d(cp_str)']

#	self.bse.apply_jacobian('constr_pts', 'd(constr_pts)/d(cp_str)', 'cp_str')
#       self.bse.vec['constr_pts'].export_tec_scatter()
#	exit()

	self.dpt_dcp = scipy.sparse.bmat(
            [
                [J, None, None],
                [None, J, None],
                [None, None, J]
            ],
            format = 'csc')
	

    def initialize(self):
	self.compute_all()
	self.func0[:] = self.func[:]

    def get_func(self):
	self.compute_all()
	return self.factor**2 * self.func0[:] - self.func[:]



class WingThicknessFunction(WingFunction):

    def __init__(self, comp, name, urange, vrange, factor):
	super(WingThicknessFunction, self).__init__(comp, name, urange, vrange, factor)
	self.size = self.num_u * self.num_v
	num_u, num_v = self.num_u, self.num_v

	self.func = numpy.zeros(num_u*num_v)
	self.func_array = self.func.reshape((num_u, num_v), order='F')

	self.func0 = numpy.array(self.func)

    def compute_all(self):
	self.pt[:] = self.dpt_dcp.dot(self.bse.vec['cp_str'].array.reshape(3*self.bse.vec['cp_str'].size, order='F'))
	num_u, num_v = self.num_u, self.num_v

	self.func_array[:,:] = 0.0
	for k in xrange(3):
	    self.func_array[:,:] += (self.pt_array[0,:,:,k] - self.pt_array[1,:,:,k])**2
#	self.func_array[:,:] += self.pt_array[0,:,:,1]


    def get_jacobian(self):
	self.compute_all()
	num_u, num_v = self.num_u, self.num_v

	nD = 6 * num_u * num_v
	Da = numpy.zeros((2,num_u,num_v,3), order = 'F')
	Di = numpy.zeros((2,num_u,num_v,3), dtype=int, order='F')
	Dj = numpy.zeros((2,num_u,num_v,3), dtype=int, order='F')

	pt_indices = numpy.array(numpy.linspace(0, 2*num_u*num_v*3-1, 2*num_u*num_v*3), int).reshape((2,num_u,num_v,3), order='F')

#	for f in [0]:
#	    for k in [1]:
#		Da[f,:,:,k] = -1
#		Di[f,:,:,k] = numpy.linspace(0, num_u*num_v-1, num_u*num_v).reshape((num_u,num_v), order='F')
#		Dj[f,:,:,k] = pt_indices[f,:,:,k]
	for f in xrange(2):
	    for k in xrange(3):
		Da[f,:,:,k] = -2 * (self.pt_array[f,:,:,k] - self.pt_array[1-f,:,:,k])
		Di[f,:,:,k] = numpy.linspace(0, num_u*num_v-1, num_u*num_v).reshape((num_u,num_v), order='F')
		Dj[f,:,:,k] = pt_indices[f,:,:,k]

	Da = Da.reshape(2*num_u*num_v*3, order='F')
	Di = Di.reshape(2*num_u*num_v*3, order='F')
	Dj = Dj.reshape(2*num_u*num_v*3, order='F')
	df_dpt = scipy.sparse.csr_matrix((Da, (Di, Dj)), shape=(num_u*num_v, 2*num_u*num_v*3))
        df_dcp = df_dpt * self.dpt_dcp
	return df_dcp
