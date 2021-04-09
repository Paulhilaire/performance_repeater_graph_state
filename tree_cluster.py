import numpy as np
import matplotlib.pyplot as plt
from fiber import Fiber
from detector import Detector

class TreeCluster(object):
    """Tree cluster state.

    Example of tree
                     O
                   /   \ b_0 = 2
                  O     O
                 /|\   /|\  b_1 = 3
                O O O O O O
                depth=2
    It is determined by its depth n and its number of arm per level given by the branching parameters:
    branches=[b_0, ..., b_{n-1}]
    The total number of photon qubits in this tree is N_qubits
    The error correction code is given for an intrinsic error epsilon_0  (see Varnava et al. (2006)).
    The number of solid state qubits required for the realization of this tree is equal to N_ss of the tree (see Buterakos et al. (2017))
    """

    def __init__(self, branches=[]):
        super(TreeCluster, self).__init__()
        self.depth = len(branches)
        self._branches = branches

    @property
    def branches(self):
        """Branching parameter of the tree graph state."""
        return self._branches

    @branches.setter
    def branches(self, value):
        self._branches = value
        self.depth = len(value)

    def add_branch(self, arm):
        """Add one level with "arm" number of arms at the bottom of the cluster tree."""
        self.depth += 1
        self.branches.append(arm)

    @property
    def N_qubits(self):
        """Number of photonic qubits"""
        N_qubits = 1
        N_prev = 1
        for i in range(self.depth):
            N_prev = N_prev * self.branches[i]
            N_qubits += N_prev
        return N_qubits

    @property
    def N_ss(self):
        """Number of solid-state qubits required to generate it."""
        if self.depth == 0:
            return 0
        else:
            return self.depth - 1

    @property
    def N_CZ(self):
        """Number of CZ operation in the general case."""
        return self.func_op(self.depth - 1)

    @property
    def N_M(self):
        """Number of single matter qubit measurement operation in the general case. It might be possible to do much less to be discussed!!!"""
        return self.func_op(self.depth - 1)

    @property
    def N_H(self):
        """Number of Hadamard operation in the general case. It might be possible to do much less to be discussed!!!"""
        return self.func_op(self.depth - 1)

    @property
    def N_Ph(self):
        """Number of Ph operation in the general case. It might be possible to do much less to be discussed!!!"""
        return 1 + self.func_op(self.depth)

    def func_op(self, n):
        """Function to count the number of operations. (f(\vec{b}, n) in the article)"""
        N = 0
        N_prev = 1
        for i in range(n):
            N_prev = N_prev * self.branches[i]
            N += N_prev
        return N
