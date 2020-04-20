from repeater_graph_state import RepeaterGraphState
from fiber import Fiber
from detector import Detector
from source import SinglePhotonSource


class RGSNetwork(object):
    """Create a Network based on repeater graph states (RGS).

    With:
    - N_RGS  RGS sources (N_RGS = 3 by default)
    - The RGS created have the properties of RGS and are always equally separated on the network distance.
      (properties of the RGS are default properties of RepeaterGraphState())
    - A fiber connect each RGS. (properties of the fiber are default properties of Fiber())
    - The detection efficiency is given by eta_d = 1 by default
    - The success probability of a Bell measurement is P_b (default P_b=0.5)
    - The time of measurement is T_mes (in s, default is 100ns)

    The total transfer time of the protocol is given by:
    time_T / P_AB
    with:
    - time_T = T_creation + T_propagation + 2 * T_measurement + T_classical, the total time of the protocol, with:
        T_creation: the repeater creation time
        T_propagation the propagation time of the photons to the quantum node.
        T_measurement: the time to read the photon detection outcome.
        T_classical: the communication time in the classical channel of all the protocols.
    - P_AB the probability of a success link between Alice and Bob given by:
        P_AB = P_2RGS ** self.N_QR
        with P_2RGS the probability of a success link between two adjacent RGS given by:
            (1 - (1 - P_bell) ** m) * P_X ** 2 * P_Z ** (2 * (m - 1))
            with:
                P_bell = P_b * (eta_t eta_d p) ** 2 the probability of a success Bell measurement
                P_X,Z = eta_t eta_d p the probability of a successful single photon measurement in the X,Z basis
    """

    def __init__(self, L_0=1, distance=1, P_b=0.5, T_mes=0,
                 RGS=RepeaterGraphState(), fiber=Fiber(), detector=Detector(), source=SinglePhotonSource()):
        super(RGSNetwork, self).__init__()
        self._L_0 = L_0
        self._distance = distance
        # self._N_RGS = N_RGS
        self._N_RGS = distance / L_0 - 1
        self.RGS = RGS
        self.fiber = fiber
        self.detector = detector
        self.P_b = P_b
        self.T_measurement = T_mes
        self.source = source

    @property
    def distance(self):
        return self._distance

    @property
    def L_0(self):
        return self._L_0

    @property
    def N_RGS(self):
        return self._N_RGS

    @distance.setter
    def distance(self, value):
        self._distance = value
        self._N_RGS = self._distance / self.L_0 - 1

    @L_0.setter
    def L_0(self, value):
        self._L_0 = value
        self._N_RGS = self._distance / self._L_0 - 1

    @property
    def eta_detector(self):
        return self.detector.det

    @property
    def loss_ph(self):
        return 1 - self.eta_detector * self.source.proba * self.fiber.transmission(self.L_0 / 2) * self.RGS.p
