from matplotlib.pylab import *
from numpy import loadtxt
from scipy.interpolate import interp1d

data = loadtxt('neb.dat').T
energy = data[2]
rc = data[1]

cs_int_f = interp1d(rc, energy, kind='cubic')


figure(figsize = (3.2,2.5), dpi = 200)

#plot(data[0], data[1])

#plot(data[0],data[2])

rc_fine = linspace(rc.min(), rc.max(), num=100)


plot(rc_fine, cs_int_f(rc_fine))

plot(rc, energy, linestyle = '', marker = 'o')

xlabel('Reaction Coordinate, ($\AA$)')
ylabel('Energy, ($eV$)')
minorticks_on()

ymax = gca().get_ylim()[1]
gca().set_ylim(0,ymax)

xmax = gca().get_xlim()[1]
gca().set_xlim(0,xmax)

tight_layout(pad = 0.2)
subplots_adjust(left= 0.20)

savefig('neb_energy_path.pdf', transparent = True)
show()
