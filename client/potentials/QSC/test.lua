function get_r(r_ij) 
    r = {0.0,0.0,0.0,
         r_ij,0.0,0.0,
         0.0,r_ij,0.0,
         0.0,0.0,r_ij,
         -r_ij,0.0,0.0,
         0.0,-r_ij,0.0,
         0.0,0.0,-r_ij}
    return r
end

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
z = {46,46,46}--,46,46,46,46}
box = {100.,100.,100.}
set_r(get_r(1.0))
set_z(z)
set_box(box)

a = 2.3
b = 2.5
range = b-a
stepsize = .05
npts = range/stepsize
energies = {}
for i=0,npts do
    r_ij = i*stepsize+a
    set_r(get_r(r_ij))
    energy, force = get_force()
    energies[i+1] = {r_ij, energy}
    --io.write(string.format("%10.4f %10.4f %10.4f\n", r_ij, energy, force[1]))
    io.write(string.format("%10.4f %10.4f\n", r_ij, energy))
    print_forces(force)
    io.write("\n")
end

--mini = 1
--for k,v in ipairs(energies) do
--    if (v[2]<energies[mini][2]) then
--        mini = k
--    end
--end
--set_r(get_r(energies[mini][1]))
--energy, force = get_force()
--io.write(string.format("Min at: %10.4f %10.4f %10.4g\n", energies[mini][1], energies[mini][2], force[1]))
