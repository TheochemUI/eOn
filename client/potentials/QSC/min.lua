function print_vector(forces)
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

function magnitude(v)
    m = 0
    for i=1,#v do
        m = m + v[i]*v[i]
    end
    m = math.sqrt(m)
    return m
end

function steepest_descent(gamma, tol, maxiters)
    --r = r0
    --set_r(r)
    U,F = get_force()
    itercount = 0
    repeat
        for i=1,#r do
            r[i] = r[i] + gamma*F[i]
        end
        set_r(r)
        U,F = get_force()
        mag_force = magnitude(F)
        itercount = itercount + 1
    until mag_force < tol or itercount > maxiters
    io.write(string.format("mag_force = %10.4f\n", mag_force))
    return r
end

function random_system(N,boxsize)
    math.randomseed(os.time())
    r = {}
    for i=1,3*N do
        r[i] = math.random()*boxsize
    end
    set_r(r)

    z = {}
    for i=1,N do z[i]=46 end
    set_z(z)

    set_box({boxsize,boxsize,boxsize})
    return r
end

--do
--    random_system(2,100)
--    rmin = steepest_descent(.01, 1e-4, 5000)
--    print("\nMinimum at:")
--    print_vector(rmin)
--end

do
    --N = 2
    --box = 2
    --r = random_system(N,box)
    --rmin = steepest_descent(.01, 1e-4, 1000)
    set_z({46,46})
    --set_box({100.,100.,100.})
    set_box({4,4,4})
    r = {0,0,0,
         2,0,0}
    set_r(r)
    --print(string.format("box: %i x %i x %i", box,box,box))
    print("r:")
    print_vector(r)
    e1,f1 = get_force()
    print("F:")
    print_vector(f1)
    dr = {}
    mag_force = magnitude(f1)

    epsilon = 1e-5
    for i=1,#r do 
        dr[i] = epsilon*(f1[i]/mag_force)
    end

    for i=1,#r do
        r[i] = r[i]+dr[i]
    end

    --print("r+dr:")
    --print_vector(r)

    set_r(r)
    e2,f2 = get_force()
    cmp1 = e1-e2
    cmp2 = 0
    for i=1,#f1 do cmp2 = cmp2 + f1[i]*dr[i] end
    print(string.format("%12s = %10.4g","U(r)",e1))
    print(string.format("%12s = %10.4g","U(r+dr)",e2))
    print(string.format("%12s = %10.4g","U(r)-U(r+dr)",cmp1))
    print(string.format("%12s = %10.4g","F.dr",cmp2))
    print(string.format("%12s = %10.4g","diff",cmp2-cmp1))
end
