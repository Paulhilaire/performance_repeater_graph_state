import numpy as np
from rgs_network import RGSNetwork
from scipy.special import factorial, binom


class ErrorModel(object):
    """Model for error correction of the tree graph state."""

    def __init__(self):
        super(ErrorModel, self).__init__()

    @property
    def qber(self):
        """Single qubit error rate."""
        return 1 - self.fidelity

    @property
    def binary_entropy(self):
        """"Binary entropy for a single qubit error rate self.qber"""
        return self.binary_entropy_func(self.qber)

    def binary_entropy_func(self, Q):
        """"Binary entropy for a single qubit error rate Q.
        Used for the binary entropy definition of the protocol."""
        if Q < 0:
            print("Q is negative in binary entropy: Q = %s" % (Q))
            print("Check manually if the error is numerical (close to 0) or more serious")
            return 0
        elif Q > 1.0:
            print("Q is above 1 in binary entropy: Q = %s" % (Q))
            print("Check manually if the error is numerical (close to 1) or more serious")
            return 1
        elif Q == 0:
            return 0
        else:
            return - Q * np.log2(Q) - (1 - Q) * np.log2(1 - Q)

    @property
    def loss_ph(self):
        """Single photon loss probability"""
        return self.network.loss_ph

    @property
    def branches(self):
        """Branches of the tree graph state used for error correction"""
        return self.network.RGS.tree.branches

    @property
    def depth(self):
        """Depth of the tree graph state used for error correction"""
        return len(self.branches)

    @property
    def b(self):
        """Adding two zeros terms to the branches to help calculations"""
        b = [x for x in self.branches]
        b.append(0)
        b.append(0)
        return b

    def R_k(self):
        """Probability of successful indirect measurement of a photon at level k in the tree."""
        R = np.zeros((len(self.b),))
        for k in range(self.depth - 1, -1, -1):
            R[k] = 1 - (1 - (1 - self.loss_ph) * (1 - self.loss_ph +
                                                  self.loss_ph * R[k + 2]) ** self.b[k + 1]) ** self.b[k]
        return R

    def S_k(self):
        """Probability of individual successful indirect measurement of a photon at level k in the tree.
        (using only one "child qubit")"""
        S = np.zeros((len(self.b),))
        for k in range(self.depth - 1, -1, -1):
            S[k] = (1 - self.loss_ph) * (1 - self.loss_ph +
                                         self.loss_ph * self.R[k + 2]) ** self.b[k + 1]
        return S

    def e_Is(self):
        """
        e_I[k]: Error of indirect measurement of a photon at level k in the tree.
        """
        self.e_I = np.zeros((len(self.b) + 1,))
        for k in range(self.depth - 1, -1, -1):
            value = 0
            for m_k in range(1, self.b[k] + 1):
                value += self.p_k(k, m_k) * self.e_Ik_mk(k, m_k)
            self.e_I[k] = value / self.R[k]
        return self.e_I

    def e_Ik_B(self, k):
        """
        e_I_B[k]: Error of individual indirect measurement of a photon at level k in the tree.
        (using only one "child qubit")"""

        value = 0

        # n_k: number of successful measurement made only directly (not indirectly)
        for n_k in range(self.b[k + 1] + 1):
            a = self.binom_proba(
                p=self.R[k + 2] / (1 - self.loss_ph + self.loss_ph * self.R[k + 2]),
                n=self.b[k+1],
                m=n_k)

            b = 0
            for i in range(n_k + 1 + 1):
                c = self.binom_proba(
                    p=self.epsilon,
                    n=n_k + 1,
                    m=i,
                )
                for j in range(self.b[k+1] - n_k + 1):
                    if i + j % 2 == 0:  # two parity errors cancel each other
                        continue
                    else:
                        d = self.binom_proba(
                            p=self.e_I[k+2],
                            n=self.b[k+1] - n_k,
                            m=j,
                        )
                        b += c * d
            value += a * b
        return value

    def binom_proba(self, p, n, m):
        """Binomial probability.
        n trials, m success with success probability p.
        """
        return binom(n, m) * p ** m * (1 - p) ** (n-m)

    def p_k(self, k, m_k):
        """probability that exactly m_k indirect measurement succeeds (perhaps with errors)."""
        return binom(self.b[k], m_k) * self.S[k] ** m_k * (1 - self.S[k]) ** (self.b[k] - m_k)

    def e_Ik_mk(self, k, m_k):
        """Error in the case of m_k indirect measurement that succeeds.
        We use a majority vote to reduce this error."""
        m_k2 = int(m_k / 2)
        value = 0
        if m_k % 2 == 1:
            # m_k is odd
            for j in range(m_k2 + 1, m_k + 1):
                value += binom(m_k, j) * self.e_Ik_B(k) ** j * (1 - self.e_Ik_B(k)) ** (m_k - j)

        else:
            # m_k is even
            for j in range(m_k2 + 1, m_k):
                value += binom(m_k - 1, j) * self.e_Ik_B(k) ** j * \
                    (1 - self.e_Ik_B(k)) ** (m_k - j - 1)
        return value

    def error_X(self):
        """Return the error corrected residual error e_z.
        """
        return self.e_I[0]

    def error_Z(self):
        """Return the error corrected residual error e_z.
        """
        error_z = 0
        for n in range(self.b[0] + 1):
            term_1 = binom(self.b[0], n) * (1 - self.R[1] / (1 - self.loss_ph +
                                                             self.loss_ph * self.R[1])) ** n
            term_2 = (self.R[1] / (1 - self.loss_ph + self.loss_ph *
                                   self.R[1])) ** (self.b[0] - n)
            e_n = 0
            for i in range(n + 1):
                a = binom(n, i) * self.epsilon ** i * (1 - self.epsilon) ** (n - i)
                b = 0
                for j in range(self.b[0] - n + 1):
                    if (i + j) % 2 == 0:
                        continue
                    b += binom(self.b[0] - n, j) * self.e_I[1] ** j * \
                        (1 - self.e_I[1]) ** (self.b[0] - n - j)
                e_n += a * b

            error_z += self.binom_proba((1 - self.R[1] / (1 - self.loss_ph + self.loss_ph * self.R[1])),
                                        self.b[0],
                                        n) * e_n

        return error_z

    def proba_X(self):
        """Calculation of the success probability of a X measurement with error correction.
        The intrinsic single qubit loss error is epsilon_0."""
        # If the tree is empty, there is only one single qubit and thus no error correction.
        if len(self.branches) == 0:
            return 1 - self.loss_ph
        # Else, the error correction is given by:
        else:
            return self.R[0]

    def proba_Z(self):
        """Calculation of the success probability of a Z measurement with error correction.
        The intrinsic single qubit loss error is epsilon_0."""
        # If the tree is empty, there is only one single qubit and thus no error correction.
        if len(self.branches) == 0:
            return 1 - self.loss_ph
        # Else, the error correction is given by:
        else:
            return (1 - self.loss_ph + self.loss_ph * self.R[1]) ** self.b[0]

    def proba_A(self):
        """Calculation of the success probability of a measurement with error correction.
        The intrinsic single qubit loss error is loss_ph."""
        # If the tree is empty, there is only one single qubit and thus no error correction.
        if len(self.branches) == 0:
            return 1 - self.loss_ph
        # Else, the error correction is given by:
        else:
            return ((1 - self.loss_ph + self.loss_ph * self.R[1]) ** self.b[0] - (self.loss_ph * self.R[1]) ** self.b[0]) * (1 - self.loss_ph + self.loss_ph * self.R[2]) ** self.b[1]

    def calculate_proba_err(self):
        """Calculations of both the errors and the measurement success proba."""
        return self.proba_X(), self.proba_Z(), self.proba_A(), self.error_X(), self.error_Z()

    def calculate_proba(self):
        """Calculations of only measurement success proba."""
        return self.proba_X(), self.proba_Z(), self.proba_A()

    def calculate_all(self):
        """Global function to calculate all that is necessary.

        If also_error is true, means that there are also single qubit errors (we do not consider perfect qubits).
        """
        self.R = self.R_k()
        if self.also_error:
            self.S = self.S_k()
            self.e_I = self.e_Is()
            self.P_X, self.P_Z, self.P_A, self.X_error, self.Z_error = self.calculate_proba_err()
        else:
            self.P_X, self.P_Z, self.P_A = self.calculate_proba()


class SingleQubitError(ErrorModel):
    """Calculations of the global protocol error using a single qubit error model."""

    def __init__(self, network=RGSNetwork(), epsilon=0):
        super(SingleQubitError, self).__init__()
        self.network = network
        self.epsilon = epsilon
        self.also_error = True

    @property
    def loss_ph(self):
        """Single photon loss probability."""
        return self.network.loss_ph

    @property
    def error_XZ(self):
        """Global X and Z error rates of the protocol for the photon pair shared between Alice and Bob."""
        return 1 / 4. - 1 / 4. * ((1 - 2 * self.epsilon) ** (2 * self.network.N_RGS + 2) * (1 - 2 * self.X_error) ** (2 * self.network.N_RGS))

    @property
    def error_Y(self):
        """Global Y error rate of the protocol for the photon pair shared between Alice and Bob."""
        a = 1 / 4. + 1 / 4. * ((1 - 2 * self.epsilon) ** (2 * self.network.N_RGS + 2)
                               * (1 - 2 * self.X_error) ** (2 * self.network.N_RGS))
        b = 1 / 2. * (1 - 2 * self.epsilon) ** (2 * self.network.N_RGS + 2) * (1 - 2 * self.X_error) ** self.network.N_RGS * \
            (1 - 2 * self.Z_error) ** ((2 * self.network.RGS.m - 2) * self.network.N_RGS)
        return a - b

    @property
    def fidelity(self):
        """Fidelity of the photon pair shared between Alice and Bob."""
        self.also_error = True
        self.calculate_all()
        # print(self.network.RGS.tree.error_tree.X_error)
        return 1 - (2 * self.error_XZ + self.error_Y)
