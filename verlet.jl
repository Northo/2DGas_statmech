function initial_values(numberOfParticles, radius)
    pos = rand((-radius:radius), 2, numberOfParticles)
    angles = rand((-pi:pi), numberOfParticles)
    return pos, angles
end


function billiard(
    numberOfParticles,
    numberOfIterations,
    radius,
    dt,
    KK;
    initial_value=nothing
)

    if (isnothing(initial_value))
        initial_value = initial_values(numberOfParticles, radius)
    end
    pos = zeros(2, numberOfParticles, numberOfIterations+1)
    velocity = Array{Float64, 3}(undef, 2, numberOfParticles, numberOfIterations+1)

    pos[:, :, 1], angles = initial_value

    vx = cos.(angles)
    vy = sin.(angles)
    velocity[:, :, 1] = transpose([vx vy])

    for i in 1:numberOfIterations, j in 1:numberOfParticles
        distance = sqrt(pos[1,j,i]^2 + pos[2,j,i]^2)
        if distance > radius
            acceleration = -KK * ((distance - radius)/distance) * pos[:, j, i]
        else
            acceleration = [0, 0]
        end

        pos[:, j, i+1] = pos[:, j, i] + velocity[:, j, i]*dt + acceleration*dt^2 / 2

        distance = sqrt(pos[1, j, i+1]^2 + pos[2, j, i+1]^2)
        if distance > radius
            acceleration2 = -KK * ((distance - radius)/distance) * pos[:, j, i+1]
        else
            acceleration2 = [0, 0]
        end

        velocity[:, j, i+1] = velocity[:, j, i] + (acceleration + acceleration2)/2 * dt
    end

    return pos, velocity
end

function potential(x, y, radius, KK)
    distance = sqrt(x^2 + y^2)
    if distance <= radius
        return 0
    end
    return (distance - radius)^2 * KK / 2
end

function energy(x, y, vel_x, vel_y, radius, KK)
    pot = potential(x, y, radius, KK)
    total_vel_square = vel_x^2 + vel_y^2

    return 0.5*total_vel_square + pot
end
