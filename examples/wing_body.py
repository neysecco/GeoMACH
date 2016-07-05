# John Hwang, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone
from GeoMACH.PGM.components import PGMjunction2, PGMtip2, PGMwing2
from GeoMACH.PSM import Airframe
import numpy


class Conventional(PGMconfiguration):

    def _define_comps(self):
        self.comps['fuse'] = PGMbody(num_x=12, num_y=4, num_z=2) #num_x is the number of surfacees
        self.comps['lwing'] = PGMwing2(num_x=4, num_z=4, left_end='tip')
        #self.comps['rwing'] = PGMwing2(num_x=4, num_z=4, right_end='junction')

        self.comps['fuse_f'] = PGMcone(self, 'fuse', 'front', 2)
        self.comps['fuse_r'] = PGMcone(self, 'fuse', 'rear', 2)
        self.comps['lwing_t'] = PGMtip2(self, 'lwing', 'left', 1)
        #self.comps['rwing_t'] = PGMtip2(self, 'rwing', 'right', 1)

        self.comps['lwing_fuse'] = PGMjunction2(self, 'fuse', 'lft', 'E', [2,1], 'lwing', 'right')
        #self.comps['rwing_fuse'] = PGMjunction2(self, 'fuse', 'rgt', 'W', [2,5], 'rwing', 'left')

    def _define_params(self):
        fuse = self.comps['fuse'].props #1st parameter: streamwise
        fuse['pos'].params[''] = PGMparameter(2, 3) #Specify the cylinder front and end
        fuse['nor'].params[''] = PGMparameter(1, 1) #Normal of each section
        fuse['scl'].params[''] = PGMparameter(1, 1)
        fuse['flt'].params[''] = PGMparameter(2, 4, pos_u=[0.28,0.53])
        fuse['pos'].params['nose'] = PGMparameter(3, 3, pos_u=[0,0.065,0.13], order_u=3)
        fuse['scl'].params['nose'] = PGMparameter(3, 1, pos_u=[0,0.07,0.14], order_u=3)
        fuse['scl'].params['tail'] = PGMparameter(2, 1, pos_u=[0.7,1.0])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.35,1.0])
        lwing['pos'].params['lin'] = PGMparameter(2, 3)
	lwing['shY','upp'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)
	lwing['shY','low'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)

        '''
        rwing = self.comps['rwing'].props
        rwing['pos'].params[''] = PGMparameter(1, 3)
        rwing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.65,1.0])
        rwing['pos'].params['lin'] = PGMparameter(2, 3)
	rwing['shY','upp'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)
	rwing['shY','low'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)
        '''

        #lwing_fuse = self.comps['lwing_fuse'].props
        #lwing_fuse['shN',''].params[''] = PGMparameter(7,7)

    def _define_dvs(self):
        dvs = self.dvs
        dvs['span'] = PGMdv((1), 23.3).set_identity_param('lwing', 'pos', 'lin', (1,2))
        dvs['mid_chord'] = PGMdv((1), 4.5).set_identity_param('lwing', 'scl', '', (1,0))
        dvs['tip_chord'] = PGMdv((1), 1.2).set_identity_param('lwing', 'scl', '', (2,0))
	dvs['shape_wing_upp'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'upp'), '')
	dvs['shape_wing_low'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'low'), '')
	#dvs['lwing_fuse_normal'] = PGMdv((7,7)).set_identity_param('lwing_fuse', ('shN', ''), '')

    def _compute_params(self):

        fuse = self.comps['fuse'].props
#	fuse['pos'].params[''].data[:,:] = [[0,0,0],[50,0,0]]
        fuse['pos'].params[''].val([[0,0,0],[50,0,0]])
        fuse['nor'].params[''].val([1.0])
        fuse['scl'].params[''].val([2.6])
        fuse['flt'].params[''].val([[0,0,0.5,0.5],[0,0,0.5,0.5]])
        fuse['pos'].params['nose'].val([[0,-1.1,0],[0,0,0],[0,0,0]])
        fuse['scl'].params['nose'].val([-2, 0, 0])
        fuse['scl'].params['tail'].val([0, -2.3])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([16,-1,2.6])
        lwing['scl'].params[''].val([10,4.5,1.2])
        lwing['pos'].params['lin'].val([[0,0,0],[16.5,4.4,23.3]])

        '''
        rwing = self.comps['rwing'].props
        rwing['pos'].params[''].val([16,-1,-2.6])
        rwing['scl'].params[''].val([1.2,4.5,10])
        rwing['pos'].params['lin'].val([[16.5,4.4,-23.3],[0,0,0]])
        '''

        #lwing_fuse = self.comps['lwing_fuse'].props
        #lwing_fuse['shN',''].params[''].data[:,:] = 0.0 #No normal perturbation initially

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['fuse'].faces['rgt'].set_option('num_cp', 'u', [4,4,4,4])
        comps['fuse'].faces['rgt'].set_option('num_cp', 'v', [18,4,4,4,4,8,4,15,4,4,10,4])
        comps['fuse'].faces['rgt'].set_option('num_pt', 'v', [40,16,16,16,16,60,16,60,16,16,70,16], both=False)
        comps['fuse'].faces['top'].set_option('num_cp', 'u', [8,8])
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', [6,4,4,20])
        comps['lwing'].faces['te'].set_option('num_cp', 'u', [9])


#=============================================#

if __name__ == '__main__':

    pgm = Conventional()
    bse = pgm.initialize()

    #pgm.comps['lwing'].set_airfoil('rae2822.dat')
    pgm.dvs['shape_wing_upp'].data[2,2] = 0.0
    pgm.dvs['span'].data[0] = 23.3
    pgm.dvs['tip_chord'].data[0] = 1.2
    #pgm.dvs['lwing_fuse_normal'].data[:,:] = 0.0 #Select normal perturbation for wing-fuselage junction
    pgm.compute_all()
    pgm.compute_normals()
    pgm.compute_all()

    #bse.vec['pt_str']._hidden[:] = False
    bse.vec['pt_str'].export_tec_str()
    bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
    bse.vec['cp_str'].export_plot3d()
