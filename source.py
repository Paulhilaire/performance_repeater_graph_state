import numpy as np

class SinglePhotonSource(object):
    """Single photon source

    It has:
    - a repetition rate r (default 1MHz)
    - a photon emission probability (default 100%).
    """

    def __init__(self, r=10 ** 6, proba=1):
        super(SinglePhotonSource, self).__init__()
        self.r = r
        self.p = proba
