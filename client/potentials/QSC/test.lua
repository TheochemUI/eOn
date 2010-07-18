r = {.0,.0,.0,1.,.0,.0}
z = {46,46}
box = {100.,100.,100.}
set_r(r)
set_z(z)
set_box(box)

a = 1.75
b = 5.0
range = b-a
stepsize = .07
npts = range/stepsize
for i=0,npts do
    r[3] = i*stepsize+a
    set_r(r)
    energy, force = get_force()
    io.write(string.format("%10.4f %10.4f\n", r[3], energy))
    for j=1,#force do
        s=string.format("%10.4f\n", force[j])
        io.write(s)
    end
end
