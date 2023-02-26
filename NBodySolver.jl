using CSV
using DataFrames

function readdata(filename)
    df = CSV.read(filename, DataFrame)
    au = 1.5e11
    m_sol = 2e30
    day = 24.0 * 60.0 * 60.0

    x = df.x * au
    y = df.y * au
    z = df.z * au
    vx = df.vx * au / day
    vy = df.vy * au / day
    vz = df.vz * au / day
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
    for i in 1:nsteps
        x += 0.5 * dt * vx
        y += 0.5 * dt * vy
        z += 0.5 * dt * vz
        ax, ay, az = calculateacceleration(x, y, z, m, g)
        vx += ax * dt
        vy += ay * dt
        vz += az * dt
        x += 0.5 * dt * vx
        y += 0.5 * dt * vy
        z += 0.5 * dt * vz
    end
end