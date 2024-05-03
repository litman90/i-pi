""" Harmonic potential """

import sys
from .dummy import Dummy_driver
import numpy as np

__DRIVER_NAME__ = "lammps"
__DRIVER_CLASS__ = "lammps_driver"

ERROR_MSG = """
lammps  driver requires specification of the inputfile
Example: python driver.py -m lammps -u -o <inputfile>
"""
from lammps import lammps, LMP_VAR_EQUAL

# from laamps import  LAMMPS_INT, LMP_STYLE_GLOBAL, LMP_VAR_ATOM
from ctypes import c_double


class lammps_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        super(lammps_driver, self).__init__(args, error_msg=ERROR_MSG)

        self.initialize_lammps()

    def initialize_lammps(self):
        """Initialize the lammps driver"""
        self.lmp = lammps()
        self.lmp.file(self._inputfile)
        self.natoms = self.lmp.get_natoms()

    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""

        if len(self.args) == 1:
            self._inputfile = self.args[0]
        else:
            sys.exit(self.error_msg)

    def send_info(self, cell, pos):
        """Send information to Lammps"""
        n3 = self.natoms * 3
        x = (n3 * c_double)()
        for i in range(n3):
            x[i] = pos.flatten()[i]
        self.lmp.scatter_atoms("x", 1, 3, x)

    def get_info(self):
        """Get information from Lammps"""
        pot = self.lmp.extract_variable("eng", None, LMP_VAR_EQUAL)
        f = self.lmp.gather_atoms("f", 1, 3)
        pot = self.lmp.extract_compute("thermo_pe", 0, 0)
        vir = np.eye(3) * 0.0  # self.lmp.extract_compute("thermo_press", 0, 0)
        force = np.zeros((self.natoms * 3))
        for i in range(3 * self.natoms):
            force[i] = f[i]

        return pot, force, vir

    def __call__(self, cell, pos):
        """Lammps potential"""

        print("\n")
        print("\n")
        print("here0", self.natoms, pos.shape)
        print("\n")
        print("\n")
        self.send_info(cell, pos)
        self.lmp.command("run 0")
        pot, force, vir = self.get_info()
        print("here", self.natoms, pos.shape, force.shape)
        # vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "nada"
        return pot, force, vir, extras
