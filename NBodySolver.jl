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
    ax = zeros(length(x))
    ay = zeros(length(x))
    az = zeros(length(x))

    for i in eachindex(x)
        axi, ayi, azi = 0., 0., 0.
        for j in eachindex(x)
            if i != j
                r_x = x[j] - x[i]
                r_y = y[j] - y[i]
                r_z = z[j] - z[i]
                r = sqrt(r_x^2 + r_y^2 + r_z^2)
                if r > 0
                    axi += g * m[j] * r_x / r^3
                    ayi += g * m[j] * r_y / r^3
                    azi += g * m[j] * r_z / r^3
                end
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
        x += 0.5 * vx * dt
        y += 0.5 * vy * dt
        z += 0.5 * vz * dt
        ax, ay, az = calculateacceleration(x, y, z, m, g)
        vx += ax * dt
        vy += ay * dt
        vz += az * dt
        x += 0.5 * vx * dt
        y += 0.5 * vy * dt
        z += 0.5 * vz * dt
    end
end