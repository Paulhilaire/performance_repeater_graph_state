import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler
from matplotlib import cm
import sys
import os
from scipy.constants import physical_constants


class Detector(object):
    """Detector with given detection efficiency.

    Possibility of improvement:
    dark count noise...
    """

    def __init__(self, det=1):
        super(Detector, self).__init__()
        if det < 0 or det > 1:
            raise "Error Detector.det should be between 0 and 1"
        self._det = det

    @property
    def det(self):
        return self._det

    @det.setter
    def det(self, value):
        if value < 0 or value > 1:
            raise "Error Detector.det should be between 0 and 1"
        self._det = value
