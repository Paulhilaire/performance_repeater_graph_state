import numpy as np


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
        """Single photon detection efficiency."""
        return self._det

    @det.setter
    def det(self, value):
        """Single photon detection efficiency."""
        if value < 0 or value > 1:
            raise "Detector.det should be between 0 and 1"
        self._det = value
