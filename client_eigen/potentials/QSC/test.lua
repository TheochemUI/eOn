function print_forces(forces)
    fs = ""
    for i,force in pairs(forces) do
        fs = fs..string.format("%10.4f", force)
        if (i%3 == 0) then
           io.write(fs.."\n")
           fs = ""
        else
           fs = fs.." "
        end
    end
end

z = {46,46}
r = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0}
box = {100.,100.,100.}
s = get_new_sim(2)
set_z(s,z)
set_box(s, box)
set_r(s,r)

a = 2.2
b = 3.5
range = b-a
stepsize = .05
npts = range/stepsize
energies = {}
for i=0,npts do
    r_ij = i*stepsize+a
    r = {0.0, 0.0, 0.0, r_ij, 0.0, 0.0}
    set_r(s, r)
    energy, force = get_force(s)
    energies[i+1] = {r_ij, energy}
    io.write(string.format("%10.4f %10.4f %10.4f\n", r_ij, energy, force[1]))
end
