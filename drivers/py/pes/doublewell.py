""" Harmonic potential """

import sys
from .dummy import Dummy_driver
import numpy as np
from ipi.utils import units
import json

# np.set_printoptions(precision=14, suppress=True,threshold='nan',linewidth=1000)

invcm2au = units.unit_to_internal("frequency", "inversecm", 1.0)
A2au = units.unit_to_internal("length", "angstrom", 1.0)

# ------------DOUBLE WELL POTENTIAL-----------------------------
#
#                                                 m^2*w_b^4
# V(x,w_b,v0) =  - 0.5 *m*w_b^2*(x-delta)^2 +    ---------- x^4
#                                                  16V0
#
# ---------------------------------------------------------


class DoubleWell_driver(Dummy_driver):
    def __init__(self, args=None):

        self.error_msg = """\nDW driver accepts 0 or 4 arguments.\nExample: python driver.py -m DoubleWell -o omega_b (cm^-1) V0 (cm^-1) mass(a.u) delta(angs) \n
        python driver.py -m DoubleWell -o 500,2085,1837,0.00 \n"""
        super(DoubleWell_driver, self).__init__(args)

    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""
        self.k = 1836 * (3800.0 / 219323.0) ** 2
        if self.args == "":
            # We used Craig's values (J. Chem. Phys. 122, 084106, 2005)
            w_b = 500 * invcm2au  # Tc = 115K
            v0 = 2085 * invcm2au
            m = 1837.36223469
            self.delta = 00
        else:
            try:
                arglist = self.args.split(",")
                param = list(map(float, arglist))
                assert len(param) == 4
                w_b = param[0] * invcm2au
                v0 = param[1] * invcm2au
                m = param[2]
                self.delta = param[3] * A2au
            except:
                sys.exit(self.error_msg)

        self.A = -0.5 * m * (w_b) ** 2
        self.B = ((m ** 2) * (w_b) ** 4) / (16 * v0)

    def __call__(self, cell, pos):
        """DoubleWell potential l"""
        pot = 0
        pos3 = pos.reshape(-1, 3)
        force3 = np.zeros(pos.shape)

        # DW
        pot += self.A * (pos3[:, 0] - self.delta) ** 2 + self.B * (pos3[:, 0] ** 4)
        force3[:, 0] = -2.0 * self.A * (pos3[:, 0] - self.delta) - 4.0 * self.B * (
            pos3[:, 0] ** 3
        )

        # Harmonic
        pot += 0.5 * self.k * (pos3[:, 1] ** 2).sum()
        pot += 0.5 * self.k * (pos3[:, 2] ** 2).sum()
        force3[:, 1] = -self.k * pos3[:, 1]
        force3[:, 2] = -self.k * pos3[:, 2]

        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "nada"
        pos = pos3.reshape(pos.shape)
        force = force3.reshape(pos.shape)

        return pot, force, vir, extras


class DoubleWell_with_friction_driver(DoubleWell_driver):
    """Adds to the double well potential the calculation of the friction tensor.

    friction(q) = eta0 [\partial sd(q) \partial q ]^2
    with
    q = position, and
    sd(q) = [1+eps1 exp( (q-0)^2 / (2deltaQ^2) ) ] + eps2 tanh(q/deltaQ)
    """

    def __init__(self, args=None):

        self.error_msg = """\nDW+fric driver excepts 8 arguments.\n
        Example: python driver.py -m DoubleWell_with_fric -o omega_c (cm^-1) V0 (cm^-1) mass delta(\AA) eta0  eps1 eps2  deltaQ      \n
        python driver.py -m DoubleWell -o 500,2085,1837,0.00,1,0,0,1\n"""
        self.args = args
        self.check_arguments()

    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""

        self.k = 1836 * (3800.0 / 219323.0) ** 2
        try:
            arglist = self.args.split(",")
            param = list(map(float, arglist))
            assert len(param) == 8
            w_b = param[0] * invcm2au
            v0 = param[1] * invcm2au
            m = param[2]
            self.delta = param[3] * A2au
            self.eta0 = param[4]
            self.eps1 = param[5]
            self.eps2 = param[6]
            self.deltaQ = param[7]
        except:
            sys.exit(self.error_msg)

        self.A = -0.5 * m * (w_b) ** 2
        self.B = ((m ** 2) * (w_b) ** 4) / (16 * v0)

    def check_dimensions(self, pos):
        """Functions that checks dimensions of the received position"""
        assert pos.shape == (1, 3)

    def SD(self, q):
        """Auxiliary function to compute friction tensor"""
        dx = q / self.deltaQ
        SD = 1.0 + self.eps1 * np.exp(-0.5 * (dx ** 2)) + self.eps2 * np.tanh(dx)
        return SD

    def dSD_dq(self, q):
        """Auxiliary function to compute friction tensor"""
        dx = q / self.deltaQ
        dsddq1 = self.eps1 * np.exp(-0.5 * (dx ** 2)) * (-dx / self.deltaQ)
        dsddq2 = self.eps2 * (1 - np.tanh(dx) ** 2) / self.deltaQ
        dSD_dq = q * (dsddq1 + dsddq2) + self.SD(q)

        return dSD_dq

    def get_friction_tensor(self, pos):
        """Function that computes spatially dependent friction tensor"""

        self.check_dimensions(pos)
        x = pos[0, 0]
        friction_tensor = np.zeros((3, 3))

        friction_tensor[0, 0] = self.eta0 * self.dSD_dq(x) ** 2
        return friction_tensor

    def __call__(self, cell, pos):
        """DoubleWell potential l"""

        pot, force, vir, extras = super(DoubleWell_with_friction_driver, self).__call__(
            cell, pos
        )

        friction_tensor = self.get_friction_tensor(pos)
        extras = json.dumps({"friction": friction_tensor.tolist()})

        return pot, force, vir, extras
