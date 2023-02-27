using CSV
using DataFrames

mutable struct Body
    x::Float64
    y::Float64
    z::Float64

    vx::Float64
    vy::Float64
    vz::Float64

    ax::Float64
    ay::Float64
    az::Float64

    m::Float64
end

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

    bodies = Vector{Body}(undef,length(x))
    for i in eachindex(x)
        bodies[i] = Body(x[i], y[i], z[i], vx[i], vy[i], vz[i], 0., 0., 0., m[i])
    end

    return bodies
end

function calculateacceleration!(bodies, g)
    Threads.@threads for i in eachindex(bodies)
        axi, ayi, azi = 0., 0., 0.
        bi = bodies[i]
        for j in eachindex(bodies)
            if i != j
                bj = bodies[j]
                r_x = bj.x - bi.x
                r_y = bj.y - bi.y
                r_z = bj.z - bi.z
                rcube = sqrt(r_x^2 + r_y^2 + r_z^2)^3
                a = g * bj.m / rcube
                axi += a * r_x
                ayi += a * r_y
                azi += a * r_z
            end
        end
        bi.ax = axi
        bi.ay = ayi
        bi.az = azi
        bodies[i] = bi
    end
end

function leapfrog!(bodies, g, dt, nsteps)
    dt_ = 0.5 * dt

    for i in 1:nsteps
        for i in eachindex(bodies)
            body = bodies[i]
            body.x += dt_ * body.vx
            body.y += dt_ * body.vy
            body.z += dt_ * body.vz
        end

        calculateacceleration!(bodies, g)

        for i in eachindex(bodies)
            body = bodies[i]
            body.vx += body.ax * dt
            body.vy += body.ay * dt
            body.vz += body.az * dt
            body.x += dt_ * body.vx
            body.y += dt_ * body.vy
            body.z += dt_ * body.vz
            bodies[i] = body
        end
    end
end