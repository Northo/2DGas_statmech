function initial_values(numberOfParticles, radius)
    radii = rand((0.0:radius), numberOfParticles)
    radii = complex.(0, radii)
    angles = rand((-pi:pi), numberOfParticles)
    pos = radii .* exp.(complex.(0, angles))

    angles_vel = rand((-pi:pi), numberOfParticles)
    return pos, angles_vel
end


function lennard_jones(pos1, pos2, epsilon)
    """Calculates the lennard_jones force"""
    distance = abs(pos2 - pos1)

    if distance > 1
        return 0
    else
        strength = 12*epsilon * (-1/distance^13 + 1/distance^7)
        return (pos2 - pos1) * strength
    end
end

function lennard_jones_potential(pos1, pos2, epsilon)
    distance = abs(pos2 - pos1)
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
    pos = zeros(Complex, numberOfParticles, numberOfIterations+1)
    velocity = Array{Complex, 2}(undef, numberOfParticles, numberOfIterations+1)

    pos[:, 1], angles = initial_value

    # Normalized velocities (average kinetic energy = 1)
    velocity[:, 1] = exp.(complex.(0, angles)) * sqrt(2)

    for i in 1:numberOfIterations, j in 1:numberOfParticles
        distance = abs(pos[j, i])
        if distance > radius
            acceleration = -KK * ((distance - radius)/distance) * pos[j, i]
        else
            acceleration = 0
        end

        for k in 1:numberOfParticles
            if k==j
                continue
            end
            lennard_jones_force = lennard_jones(
                pos[j, i],
                pos[k, i],
                epsilon,
            )
            acceleration += lennard_jones_force
        end

        pos[j, i+1] = pos[j, i] + velocity[j, i]*dt + acceleration*dt^2 / 2

        distance = abs(pos[j, i+1])
        if distance > radius
            acceleration2 = -KK * ((distance - radius)/distance) * pos[j, i+1]
        else
            acceleration2 = 0
        end

        for k in 1:numberOfParticles
            if k==j
                continue
            end
            lennard_jones_force = lennard_jones(
                pos[j, i+1],
                pos[k, i],
                epsilon,
            )
            acceleration2 += lennard_jones_force
        end

        velocity[j, i+1] = velocity[j, i] + (acceleration + acceleration2)/2 * dt
    end

    return pos, velocity
end

function potential(pos, radius, KK)
    distance = abs(pos)
    if distance <= radius
        return 0
    end
    return (distance - radius)^2 * KK / 2
end

function energy(pos, vel, radius, KK)
    pot = potential(pos, radius, KK)
    return 0.5*abs2(vel) + pot
end

function relative_error(a::Vector)
    """Returns a list of the relative error with respect to the first element"""
    return (a .- a[1])/a[1]
end

function get_interaction_potential(pos, epsilon)
    """Given list of positions, with shape pos[num_particles, num_iterations],
    returns list with interaction potential for each particle, with the same shape"""

    interaction_energy = zero(pos)
    num_particles, num_iterations = size(pos)

    for i in 1:num_particles, j in 1:num_particles
        if i==j
            continue
        end
        interaction_energy[i, :] += lennard_jones_potential.(
            pos[i, :],
            pos[j, :],
            epsilon,
        )
    end

    return interaction_energy
end
