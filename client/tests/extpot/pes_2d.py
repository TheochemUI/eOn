import gauss_peak
import numpy as np
#import pylab

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

#    def map_out_pes(self, interval_x, interval_y, res):
#
#        landscape = self._build_landscape()
#        dx = (interval_x[1] - interval_x[0]) / float(res)
#        dy = (interval_y[1] - interval_y[0]) / float(res)
#
#        data = []
#        for i in range(res):
#            data_line = []
#            for j in range(res):
#                x = interval_x[0] + i * dx
#                y = interval_y[0] + j * dy
#                r = np.array([x,y])
#                self._calculate_landscape(r, landscape)
#                data_line.append(self.get_E())
#            data.append(data_line)
#
#            x_list = []
#            y_list = []
#        for j in range(res):
#            x_list.append(interval_x[0] + j * dx)
#            y_list.append(interval_y[0] + j * dy)
#        x_data, y_data = pylab.meshgrid(np.array(x_list), np.array(y_list))
#
#        return np.array(data), x_data, y_data
#
#    def plot(self, interval_x, interval_y, res):
#        data = self.map_out_pes(interval_x, interval_y, res)
#        pylab.pcolor(data[1], data[2], data[0].T)
#        pylab.xlim([interval_x[0],interval_x[1]])
#        pylab.ylim([interval_y[0],interval_y[1]])
#        pylab.axis('equal')

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
