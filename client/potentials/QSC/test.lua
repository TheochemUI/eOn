function get_r(r_ij) 
    r = {0.0,0.0,0.0,
         r_ij,0.0,0.0,
         0.0,r_ij,0.0,
         0.0,0.0,r_ij}
    return r
end
z = {46,46}
box = {100.,100.,100.}
set_r(get_r(1.0))
set_z(z)
set_box(box)

a = 2.0
b = 5.0
range = b-a
stepsize = .001
npts = range/stepsize
energies = {}
for i=0,npts do
    r_ij = i*stepsize+a
    set_r(get_r(r_ij))
    energy, force = get_force()
    energies[i+1] = {r_ij, energy}
    --io.write(string.format("%10.4f %10.4f %10.4f\n", r_ij, energy, force[1]))
end

mini = 1
for k,v in ipairs(energies) do
    if (v[2]<energies[mini][2]) then
        mini = k
    end
end
set_r(get_r(energies[mini][1]))
energy, force = get_force()
print("Minimum Energy at: ", energies[mini][1], energies[mini][2], force[1])
--print("Forces:")
--for i=1,#force do
--    s=string.format("%10.4f\n", force[i])
--    io.write(s)
--end
