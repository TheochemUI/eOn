import scipy

class GaussPeak():
    def __init__(self, r0, sigma, height):
        self.setParameters(r0, sigma, height)
        return

    def setParameters(self, r0, sigma, height):
        self._r0 = r0
        self._height = height
        sigma = scipy.array(sigma)
        if scipy.size(sigma) == 1:
            self._sigmaX2 = sigma**2
            self._sigmaY2 = sigma**2
            self._sigmaX4 = sigma**4
            self._sigmaY4 = sigma**4
        else:
            self._sigmaX2 = sigma[0]**2
            self._sigmaY2 = sigma[1]**2
            self._sigmaX4 = sigma[0]**4
            self._sigmaY4 = sigma[1]**4
        self._sigmaX2sigmaY2 = self._sigmaX2 * self._sigmaY2
        return


    def _getValue(self, r):

        x = r[0] - self._r0[0]
        y = r[1] - self._r0[1]

        v = self._height * scipy.exp(-((x**2 / self._sigmaX2) + (y**2 / self._sigmaY2)))
        return (v, [x, y])


    def getValue(self, r):
        return self._getValue(r)[0]


    def getFirstDerivative(self, r):
        v, [x, y] = self._getValue(r)

        dfdx = -2. * x / self._sigmaX2 * v
        dfdy = -2. * y / self._sigmaY2 * v

        return scipy.array([dfdx, dfdy])


    def getSecondDerivative(self, r):
        v, [x, y] = self._getValue(r)

        d2fd2x = (-2. / self._sigmaX2 + 4. * x**2 / self._sigmaX4) * v
        d2fd2y = (-2. / self._sigmaY2 + 4. * y**2 / self._sigmaY4) * v
        d2fdxdy = 4. * x * y / self._sigmaX2sigmaY2 * v

        return scipy.array([[ d2fd2x, d2fdxdy],
                            [d2fdxdy, d2fd2y ]])
