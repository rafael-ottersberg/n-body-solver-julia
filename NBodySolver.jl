using CSV
using DataFrames

function readdata(filename)
    df = CSV.read(filename, DataFrame)
    au = 1.5e11
    m_sol = 2e30
    day = 24.0 * 60.0 * 60.0
    
    scaledist(vector) = vector * au
    scalespeed(vector) = vector * au / day

    x = scaledist(df.x)
    y = scaledist(df.y)
    z = scaledist(df.z)
    vx = scalespeed(df.vx)
    vy = scalespeed(df.vy)
    vz = scalespeed(df.vz)

    m = df.m * m_sol

    return x, y, z, vx, vy, vz, m
end

function calculateacceleration(x, y, z, m, g)
    n = length(x)
    ax = Vector{Float64}(undef, n)
    ay = Vector{Float64}(undef, n)
    az = Vector{Float64}(undef, n)

    for i in 1:n
        axi, ayi, azi = 0., 0., 0.
        for j in 1:n
            if i != j
                r_x = x[j] - x[i]
                r_y = y[j] - y[i]
                r_z = z[j] - z[i]
                rcube = sqrt(r_x^2 + r_y^2 + r_z^2)^3
                a = g * m[j] / rcube
                axi += a * r_x
                ayi += a * r_y
                azi += a * r_z
            end
        end
        ax[i] = axi
        ay[i] = ayi
        az[i] = azi
    end
    return ax, ay, az
end

function leapfrog!(x, y, z, vx, vy, vz, m, g, dt, nsteps)
    dt_ = 0.5 * dt

    for i in 1:nsteps
        x = dt_ * vx
        y = dt_ * vy
        z = dt_ * vz
        ax, ay, az = calculateacceleration(x, y, z, m, g)
        vx += ax * dt
        vy += ay * dt
        vz += az * dt
        x += dt_ * vx
        y += dt_ * vy
        z += dt_ * vz
    end
end