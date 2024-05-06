""" Harmonic potential """

import sys
from .dummy import Dummy_driver
import numpy as np
from ipi.utils import units

__DRIVER_NAME__ = "lammps"
__DRIVER_CLASS__ = "lammps_driver"

ERROR_MSG = """
Lammps  driver requires specification of the inputfile
Note: the lammps input file should contain all the information as usual
except the 'fix ipi' and 'run' commands

Example: python driver.py -m lammps -u -o <inputfile>
"""
from lammps import lammps, LMP_VAR_EQUAL
from lammps.constants import (
    LMP_STYLE_ATOM,
    LMP_TYPE_VECTOR,
    LMP_STYLE_GLOBAL,
    LMP_TYPE_ARRAY,
)

# from laamps import  LAMMPS_INT, LMP_STYLE_GLOBAL, LMP_VAR_ATOM
from ctypes import c_double

units_lammps = {}
units_lammps["electron"] = {
    "boltz" == 3.16681534e-6,
    "nktv2p" == 2.94210108e13,
    "angstrom" == 1.88972612,
}
units_lammps["metal"] = {
    "boltz" == 8.617343e-5,
    "nktv2p" == 1.6021765e6,
    "angstrom" == 1.0,
}
A2au = units.unit_to_internal("length", "angstrom", 1.0)
eV2au = units.unit_to_internal("energy", "electronvolt", 1.0)


class lammps_driver(Dummy_driver):

    def __init__(self, args=None, verbose=False,error_msg=ERROR_MSG):
        super(lammps_driver, self).__init__(args, error_msg=error_msg)
        self.posconv = 0.52917721
        self.potconv = 3.1668152e-06
        self.initialize_lammps()

    def initialize_lammps(self):
        """Initialize the lammps driver"""
        self.lmp = lammps()
        self.lmp.file(self._inputfile)
        self.units = {}
        self.units["name"] = self.lmp.extract_global("units")
        assert self.units["name"] in [
            "electron",
            "metal",
        ], 'This driver is only implemented for Lammps "electron" and "metal" units'
        if self.units["name"] == "electron":
            self.units["pos"] = 1.0
            self.units["energy"] = 1.0
            self.units["force"] = 1.0
            self.units["virial"] = 1.0

        elif self.units["name"] == "metal":
            self.units["pos"] = 1.0 / A2au
            self.units["energy"] = 1.0 / eV2au
            self.units["force"] = A2au / eV2au
            self.units["virial"] = 1.0

        self._initial_box = self.lmp.extract_box()

        self.natoms = self.lmp.get_natoms()
        self.natoms3 = self.natoms * 3

        self.lmp.command("velocity all set 0.0 0.0 0.0")
        # self.lmp.command("compute IPI_TEMP all temp")
        # self.lmp.command("compute IPI_PRESS all pressure IPI_TEMP virial")

    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""

        if len(self.args) == 1:
            self._inputfile = self.args[0]
        else:
            sys.exit(self.error_msg)

    def send_info(self, cell, pos):
        """Send information to Lammps"""
        x = (self.natoms3 * c_double)()
        p = pos.flatten()
        for i in range(self.natoms3):
            x[i] = p[i]

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
            raise NotImplementedError("This hasnt been tested.")
            self.lmp.command(
                "change_box all x final {} {} y final {} {} z final {} {} xy final {} xz final {} yz final {}".format(
                    boxlo[0],
                    boxhi[0],
                    boxlo[1],
                    boxhi[1],
                    boxlo[2],
                    boxhi[2],
                    xy,
                    xy,
                    yz,
                )
            )

    def get_info(self):
        """Get information from Lammps"""
        # pot = self.lmp.extract_variable("eng", None, LMP_VAR_EQUAL)
        pot = self.lmp.extract_compute("thermo_pe", 0, 0)

        f = self.lmp.gather_atoms("f", 1, 3)
        force = np.zeros((self.natoms * 3))
        for i in range(3 * self.natoms):
            force[i] = f[i]

        current_vol = self.lmp.get_thermo("vol")
        vir = np.eye(3) * 0.0

        # aux_vir = self.lmp.numpy.extract_compute("IPI_PRESS",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR)
        # vir[0,0] = aux_vir[0]*current_vol
        # vir[1,1] = aux_vir[1]*current_vol
        # vir[2,2] = aux_vir[2]*current_vol
        # vir[0,1] = aux_vir[3]*current_vol
        # vir[0,2] = aux_vir[4]*current_vol
        # vir[1,2] = aux_vir[5]*current_vol

        # vir[1,0] = vir[0,1]
        # vir[2,0] = vir[0,2]
        # vir[2,1] = vir[1,2]

        return pot, force, vir

    def __call__(self, cell, pos):
        """Lammps potential"""
        # Check box is not changing since virial is not implemented,after that this check should be deleted
        self.box = self.lmp.extract_box()
        assert (
            self.box == self._initial_box
        ), "Sorry the communication of the  virial is not implemented yet, until that happens the box must not change during the simulation"

        # Send information
        self.send_info(cell * self.units["pos"], pos * self.units["pos"])

        # Compute quantities this configuration
        self.lmp.command("run 0")

        # Gather quantities this configuration
        pot, force, vir = self.get_info()
        pot /= self.units["energy"]
        force /= self.units["force"]
        vir /= self.units["virial"]
        extras = "nada"

        return pot, force, vir, extras

