# John Hwang, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing2, PGMtip2
from GeoMACH.PSM import Airframe
import numpy


class Wing(PGMconfiguration):

    def _define_comps(self):
        self.comps['wing'] = PGMwing2(num_x=2, num_z=4, left_closed=True)
        self.comps['tip'] = PGMtip2(self, 'wing', 'left', 0.1)

    def _define_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''] = PGMparameter(3, 3, pos_u=[0,0.37,1.0])
        wing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.37,1.0])

    def _compute_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''].data[0, :] = [904.294, 174.126, 0.0]
        wing['pos'].params[''].data[1, :] = [1225.82, 181.071, 427.999]
        wing['pos'].params[''].data[2, :] = [1780.737, 263.827, 1156.753]
        wing['scl'].params[''].data[:, 0] = [536.181, 285.782, 107.4]
        return [], [], []

    def _set_bspline_options(self):
        wing = self.comps['wing'].faces
        wing['upp'].set_option('num_cp', 'u', [10])
        wing['upp'].set_option('num_cp', 'v', [10])
        wing['low'].set_option('num_cp', 'u', [10])
        wing['low'].set_option('num_cp', 'v', [10])
        wing['te'].set_option('num_cp', 'u', [8])
        wing['te'].set_option('num_cp', 'v', [8])
        wing['le'].set_option('num_cp', 'u', [8])
        wing['le'].set_option('num_cp', 'v', [8])
        pass


# EXECUTION

if __name__ == '__main__':

    pgm = Wing()
    bse = pgm.initialize()
    #pgm.comps['wing'].set_airfoil(filename=['naca0019','rae2822.dat'])
    pgm.compute_all()
    bse.vec['pt_str'].export_tec_str()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
    exit()
