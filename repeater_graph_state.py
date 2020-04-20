from scipy.interpolate import interp2d
from utils import template_plot
import os
import sys
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import cycler
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import physical_constants

from fiber import Fiber
from detector import Detector
from tree_cluster import TreeCluster
dir_path = os.path.dirname(os.path.realpath(__file__))


class RepeaterGraphState(object):
    """Define a Repeater graph state (RGS).

    With:
    - 2 * m arms  (m=3 by default)
    - Collection efficiency of the generated photon : p  (p=1 by default)
    - 1st-leaf qubits encoded with an error-corrected tree cluster state (single qubit by default)

    The operation times (expressed in s!!!) for generating it are:
    - T_ph for the photon emission time (default 1ns)
    - T_CZ for the CZ gate time (default 1ns)
    - T_H for the Hadamard gate time (default 1ns)
    - T_Me for the emitter qubit measurement time (default 1ns)

    The coherence time of the emitter and ancillary qubits are:
    - T_2e for the emitter (default 1ns)
    - T_2 for the ancilla (default 1ns)

    The gate fidelities are:
    - F_ph for the photon emission (default 1)
    - F_CZ for the CZ gate (default 1)
    - F_H for the Hadamard gate (default 1)
    - F_Me for the emitter qubit measurement (default 1)
    They are not used for the moment in the program.
    """

    def __init__(self, m=3, p=1, T_ph=0,
                 T_CZ=10**-6, T_H=0, T_Me=0,
                 T_2='infinity', T_2e='infinity',
                 F_ph=1, F_CZ=1, F_H=1, F_Me=1,
                 tree=TreeCluster()):
        super(RepeaterGraphState, self).__init__()
        self.m = m
        self.p = p
        self.T_ph = T_ph
        self.T_CZ = T_CZ
        self.T_H = T_H
        self.T_Me = T_Me
        self.T_2 = T_2
        self.T_2e = T_2e
        self.F_ph = F_ph
        self.F_CZ = F_CZ
        self.F_H = F_H
        self.F_Me = F_Me
        self.tree = tree

    @property
    def T_arm(self):
        """Time to generate one arm of the RGS.

        Works for tree of depth 2."""
        time_m = (self.tree.N_M + 2) * self.T_Me
        time_ph = (self.tree.N_Ph + 1) * self.T_ph
        time_cz = (self.tree.N_CZ + 2) * self.T_CZ
        time_h = self.tree.N_H * self.T_H
        return time_m + time_h + time_cz + time_ph

    @property
    def T(self):
        """Generation time of an RGS."""
        return self.T_arm * (2 * self.m) + self.T_Me

    @property
    def N_ph(self):
        """Number of photon in the RGS."""
        return 2 * self.m * (1 + self.tree.N_qubits)

    @property
    def N_ss(self):
        """Number of matter qubits required for the generation of the RGS."""
        return 2 + self.tree.N_ss
