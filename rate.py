import numpy as np
from rgs_network import RGSNetwork
from error_model import SingleQubitError


class Rate(object):
    """Calculate the protocol rates for a given RGS network.

    A network and an error model should be specified"""

    def __init__(self, network=RGSNetwork(), error=SingleQubitError()):
        # super(Rate, self).__init__()
        self.network = network
        self.error = error

    def get_rates(self):
        """Calculate the probability of generating an entanglement connection "proba"
        and the rate and secret key rate (skr) for a given network and a given error model."""
        self.error.calculate_all()
        self.proba = self.P_AB()
        self.rate = self.proba / self.time_T
        self.skr = max(0, self.rate * (1 - 2 * self.error.binary_entropy))

    def P_AB(self):
        """Probability to generate a link between Alice and Bob.
        It is the probability that an entanglement connection is realized between all RGS at intermediary nodes."""
        return self.P_2RGS() ** (self.network.N_RGS + 1)

    def P_2RGS(self):
        """ Probability of generating an entanglement connection between two adjacent RGS."""
        bell_measurements = (1 - (1 - self.P_bell()) ** self.network.RGS.m)
        x_measurements = self.error.P_X ** 2
        z_measurements = self.error.P_Z ** (2 * (self.network.RGS.m - 1))
        return bell_measurements * x_measurements * z_measurements

    def P_bell(self):
        """Success probability of a photonic Bell state measurement."""
        return (1 - self.network.loss_ph) ** 2 * self.network.P_b

    @property
    def time_T(self):
        """Repetition time of the protocol.
        It is only limited by the generation time of an RGS
        """
        return self.network.RGS.T

    @property
    def proba_ch(self):
        """Probability per channel use."""
        return self.proba / (2 * self.network.RGS.N_ph)

    @property
    def rate_ch(self):
        """Rate per channel use."""
        return self.rate / (2 * self.network.RGS.N_ph)

    @property
    def skr_ch(self):
        """Secret key rate per channel use."""
        return self.skr / (2 * self.network.RGS.N_ph)

    @property
    def proba_ph(self):
        """Probability per photon used."""
        return self.proba / (self.network.RGS.N_ph * self.network.N_RGS)

    @property
    def rate_ph(self):
        """Rate per photon used."""
        return self.rate / (self.network.RGS.N_ph * self.network.N_RGS)

    @property
    def skr_ph(self):
        """Secret key rate per photon used."""
        return self.skr / (self.network.RGS.N_ph * self.network.N_RGS)

    @property
    def proba_ss(self):
        """Probability per matter qubit used."""
        return self.proba / (self.network.RGS.N_ss * self.network.N_RGS)

    @property
    def rate_ss(self):
        """Rate per matter qubit used."""
        return self.rate / (self.network.RGS.N_ss * self.network.N_RGS)

    @property
    def skr_ss(self):
        """Secret key rate per matter qubit used."""
        return self.skr / (self.network.RGS.N_ss * self.network.N_RGS)

    @property
    def proba_qn(self):
        """Probability per quantum node."""
        return self.proba / self.network.N_RGS

    @property
    def rate_qn(self):
        """Rate per quantum node."""
        return self.rate / self.network.N_RGS

    @property
    def skr_qn(self):
        """Secret key rate per quantum node."""
        return self.skr / self.network.N_RGS
