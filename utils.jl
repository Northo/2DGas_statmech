function initial_values(numberOfParticles, radius)
    pos = rand((-radius:radius), 2, numberOfParticles)
    angles = rand((-pi:pi), numberOfParticles)
    return pos, angles
end


function lennard_jones(pos1_x, pos1_y, pos2_x, pos2_y, epsilon)
    """Calculates the lennard_jones force"""

    dx = pos2_x - pos1_x
    dy = pos2_y - pos1_y

    distance = sqrt(
        (dx)^2 + (dy)^2
    )
    if distance > 1
        return [0, 0]
    else
        strength = 12*epsilon * (-1/distance^13 + 1/distance^7)
        return [dx, dy] * strength
    end
end

function lennard_jones_potential(pos1_x, pos1_y, pos2_x, pos2_y, epsilon)
    dx = pos2_x - pos1_x
    dy = pos2_y - pos1_y

    distance = sqrt(
        (dx)^2 + (dy)^2
    )
    if distance > 1
        return 0
    else
        return epsilon * (1/distance^12 - 2/distance^6 + 1)
    end
end

function billiard(
    numberOfParticles,
    numberOfIterations,
    radius,
    dt,
    KK;
    initial_value=nothing,
    epsilon=1,
)

    if (isnothing(initial_value))
        initial_value = initial_values(numberOfParticles, radius)
    end
    pos = zeros(2, numberOfParticles, numberOfIterations+1)
    velocity = Array{Float64, 3}(undef, 2, numberOfParticles, numberOfIterations+1)

    pos[:, :, 1], angles = initial_value

    # Normalized velocities (average kinetic energy = 1)
    vx = cos.(angles) * sqrt(2)
    vy = sin.(angles) * sqrt(2)
    velocity[:, :, 1] = transpose([vx vy])

    for i in 1:numberOfIterations, j in 1:numberOfParticles
        distance = sqrt(pos[1,j,i]^2 + pos[2,j,i]^2)
        if distance > radius
            acceleration = -KK * ((distance - radius)/distance) * pos[:, j, i]
        else
            acceleration = [0, 0]
        end

        for k in 1:numberOfParticles
            if k==j
                continue
            end
            lennard_jones_force = lennard_jones(
                pos[1, j, i],
                pos[2, j, i],
                pos[1, k, i],
                pos[2, k, i],
                epsilon,
            )
            acceleration += lennard_jones_force
        end

        pos[:, j, i+1] = pos[:, j, i] + velocity[:, j, i]*dt + acceleration*dt^2 / 2

        distance = sqrt(pos[1, j, i+1]^2 + pos[2, j, i+1]^2)
        if distance > radius
            acceleration2 = -KK * ((distance - radius)/distance) * pos[:, j, i+1]
        else
            acceleration2 = [0, 0]
        end

        for k in 1:numberOfParticles
            if k==j
                continue
            end
            lennard_jones_force = lennard_jones(
                pos[1, j, i+1],
                pos[2, j, i+1],
                pos[1, k, i],
                pos[2, k, i],
                epsilon,
            )
            acceleration2 += lennard_jones_force
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

function estimate_distribution(collection)
    """Given a collection, estimate the distribution of the quantity it containts.
    The colletion must have the shape (numberOfParticles, numberOfIterations)"""

    return mean(collection)
end
