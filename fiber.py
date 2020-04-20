import numpy as np
from scipy.constants import physical_constants


class Fiber(object):

    """
    Define a fiber to calculate its photon loss.
    -> set or get db_km/att_distance are used to define the properties of
    the fiber attenuation in distance or in dB/km
    -> Transmission calculate the fiber transmission
    """

    def __init__(self, att=0.2):
        """
        Define a fiber and sets its default attenuation fiber value
        to 0.2dB/km (typical fiber loss in telecommunication fiber)
        """
        super(Fiber, self).__init__()
        self._attenuation_fiber_dB_km = att
        self._attenuation_distance = self.convert(
            self._attenuation_fiber_dB_km)
        # Speed of light in a fiber
        self.light_speed = 2 / 3. * physical_constants["speed of light in vacuum"][0] / 10 ** 3

    def convert(self, f):
        return 10 / (np.log(10) * f)

    def transmission(self, distance):
        """Fiber transmission for a given distance"""
        return np.exp(- distance / self.attenuation_distance)

    @property
    def attenuation_fiber_dB_km(self):
        return self._attenuation_fiber_dB_km

    @property
    def attenuation_distance(self):
        return self._attenuation_distance

    @attenuation_fiber_dB_km.setter
    def attenuation_fiber_dB_km(self, value):
        self._attenuation_fiber_dB_km = value
        self._attenuation_distance = self.convert(value)

    @attenuation_distance.setter
    def attenuation_distance(self, value):
        self._attenuation_distance = value
        self._attenuation_fiber_dB_km = self.convert(value)

    @property
    def t_att(self):
        return self.attenuation_distance / self.light_speed
