""" Harmonic potential """

import sys
from .dummy import Dummy_driver
import numpy as np
from ipi.utils import units
from .lammps import lammps_driver
from ctypes import c_double
__DRIVER_NAME__ = "lammps_ghost"
__DRIVER_CLASS__ = "lammps_ghost_driver"

ERROR_MSG = """
Lammps  driver requires specification of the inputfile
Note: the lammps input file should contain all the information as usual
except the 'fix ipi' and 'run' commands

Example: python driver.py -m lammps -u -o <inputfile>,<displacement_z>
"""
class lammps_ghost_driver(lammps_driver):

    def __init__(self, args=None, verbose=False):
        super(lammps_ghost_driver, self).__init__(args, verbose=verbose, error_msg=ERROR_MSG)

    def initialize_lammps(self):
        super(lammps_ghost_driver, self).initialize_lammps()
        if(np.absolute(self._displacement) >0.0001):
         self.natoms3_real  = self.natoms3//2
         self.natoms3_ghost = self.natoms3//2
         self.natoms_real  = self.natoms//2
         self.natoms_ghost = self.natoms//2
        else:
         self.natoms3_real  = self.natoms3
         self.natoms3_ghost = 0
         self.natoms_real  = self.natoms
         self.natoms_ghost = 0

       


    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""

        if len(self.args) == 2:
            self._inputfile = self.args[0]
            self._displacement = float(self.args[1])
            self.disp_v = np.array([0,0,self._displacement])
        else:
            sys.exit(self.error_msg)

    def send_info(self, cell, pos):
        """Send information to Lammps"""
        x = (self.natoms3 * c_double)()
        p = pos.flatten()

        #Real system
        for i in range(self.natoms3_real):
            x[i] = p[i]

        #Ghost system
        for i in range(self.natoms_ghost):
         for ii in range(3):
            print(i,ii,self.natoms3_real+3*i+ii,self.natoms3)
            x[self.natoms3_real+3*i+ii] = p[3*i+ii] + self.disp_v[ii] 

        self.lmp.scatter_atoms("x", 1, 3, x)
        # boxlo = [-0.5*cell[0,0],-0.5*cell[1,1],-0.5*cell[2,2]]
        # boxhi = [0.5*cell[0,0],0.5*cell[1,1],0.5*cell[2,2]]
        boxlo = [0.0, 0.0, 0.0]
        boxhi = [cell[0, 0], cell[1, 1], cell[2, 2]]
        xy = cell[0, 1]
        xz = cell[0, 2]
        yz = cell[1, 2]
        if (
            np.absolute(xy) < 0.001
            and np.absolute(xz) < 0.001
            and np.absolute(yz) < 0.001
        ):
            self.lmp.command(
                "change_box all x final {} {} y final {} {} z final {} {}".format(
                    boxlo[0], boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]
                )
            )
        else:
            raise NotImplementedError("Only orthorhombic boxes are implemented.")

    def get_info(self):
        """Get information from Lammps"""
        # pot = self.lmp.extract_variable("eng", None, LMP_VAR_EQUAL)
        pot = self.lmp.extract_compute("thermo_pe", 0, 0)

        f = self.lmp.gather_atoms("f", 1, 3)
        force = np.zeros((self.natoms_real * 3))
        for i in range(self.natoms3_real):
            force[i] = f[i]
            print(i,force[i])
        current_vol = self.lmp.get_thermo("vol")
        vir = np.eye(3) * 0.0

        return pot, force, vir

