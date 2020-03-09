function generate_initial_positions(numberOfParticles, radius)
    radii = rand((0.0:radius), numberOfParticles) .+ 0im
    angles = rand((-pi:pi), numberOfParticles)*im
    pos = radii .* exp.(angles)

    return pos
end


function generate_initial_velocities(numberOfParticles)
    angles = rand((-pi:pi), numberOfParticles)*im
    # Normalized velocities (average kinetic energy = 1)
    return exp.(angles) * sqrt(2)
 end


function safe_initial_values(numberOfParticles, radius)
    """Deprecated. Generate 'safe' initial values, where noe particles starts too close"""
    angles = range(0, 2*pi, length=numberOfParticles)*im
    radii = range(0.1, 0.9, length=numberOfParticles) .* radius .+ 0im
    pos = radii .* exp.(angles)

    angles_vel = rand((-pi:pi), numberOfParticles)*im
    velocities = exp.(angles_vel)

    return pos, velocities
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
    initial_positions=nothing,
    initial_velocities=nothing,
    epsilon=1,
)

    if (isnothing(initial_positions))
        initial_positions = generate_initial_positions(numberOfParticles, radius)
    end

    if (isnothing(initial_velocities))
        initial_velocities = generate_initial_velocities(numberOfParticles)
    end

    pos = zeros(Complex, numberOfParticles, numberOfIterations+1)
    velocity = Array{Complex, 2}(undef, numberOfParticles, numberOfIterations+1)

    pos[:, 1] = initial_positions
    velocity[:, 1] = initial_velocities

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

function plot_velocity_distributions(vel)
    vel_distribution_fig, vel_distribution_ax = plt.subplots()

    vel_x = real.(vel)
    v_min, v_max = minimum(vel_x), maximum(vel_x)
    v_list = range(v_min, v_max, length=100)

    # Boltzman dist.
    kbT = mean(x->x^2, vel_x)  # From equipartition thm.
    boltzman_std = sqrt(kbT)
    boltzman_dist = Normal(0, boltzman_std)

    plt.hist(collect(Iterators.flatten(vel_x)), bins=40, density=true, label="Simulated distribution")
    plt.plot(v_list, Distributions.pdf.(boltzman_dist, v_list), label="Boltzman distribution")
    if FIT
        plt.plot(v_list, Distributions.pdf.(fit_velocity_distribution, v_list), label="Fitted curve")
    end
    plt.legend()
    plt.show()

end
