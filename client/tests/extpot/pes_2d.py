import gauss_peak
import numpy as np

class PES_2D():
    def __init__(self, x=None, y=None):
        self._E = None
        self._fx = None
        self._fy = None

        landscape = self._build_landscape()

        if x is not None:
            r = np.array([x,y])
            self._calculate_landscape(r, landscape)


    def get_E(self):
        return self._E

    def get_fx(self):
        return self._fx

    def get_fy(self):
        return self._fy

    ########################## THE LANDSCAPE #####################

    def _calculate_landscape(self, r, landscape):
        self._E = 0.
        self._fx = 0.
        self._fy = 0.
        for component in landscape:
            self._E += component.getValue(r)
            f = -component.getFirstDerivative(r)
            self._fx += f[0]
            self._fy += f[1]


    def _build_landscape(self):
        landscape = []
        landscape.append(gauss_peak.GaussPeak([ 0.,     -1.], [1.,2.], -12.))
        landscape.append(gauss_peak.GaussPeak([ 0.,      1.], [1.,2.], -12.))
        landscape.append(gauss_peak.GaussPeak([ 0.,      .0], [1.,2.], 13.))

        return landscape
