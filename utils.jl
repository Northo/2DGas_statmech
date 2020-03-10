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


function safe_initial_positions(numberOfParticles, radius)
    """Deprecated. Generate 'safe' initial values, where noe particles starts too close"""
    radii = range(0.1, 0.9, length=numberOfParticles) .* radius .+ 0im
    angles = range(0, 2*pi, length=numberOfParticles)*im
    pos = radii .* exp.(angles)
    return pos
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

    interaction_energy = zeros(Float64, size(pos))
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

function plot_velocity_distributions(vel, fit_distribution=nothing)
    vel_distribution_fig, vel_distribution_ax = plt.subplots()

    vel_x = real.(vel)
    v_lim = maximum(abs.(vel_x))
    v_list = range(-v_lim, v_lim, length=100)

    # Boltzman dist.
    boltzman_dist = Normal(0, 1)

    plt.hist(collect(Iterators.flatten(vel_x)), bins=40, density=true, label="Simulated distribution")
    plt.plot(
        v_list,
        Distributions.pdf.(boltzman_dist, v_list),
        label="Boltzman distribution",
    )
    if !isnothing(fit_distribution)
        plt.plot(v_list, fit_distribution(v_list), label="Fitted curve")
    end

    # Add some info

    plt.text(
        -7, 0.2,
        @sprintf("Mean \$v_x\$: %.3f\nRMS \$v\$: %.3f", mean(vel_x), sqrt(mean(abs2, vel))),
        bbox = Dict("boxstyle"=>"round", "ec"=>(1., 0.5, 0.5), "fc"=>(1., 0.8, 0.8))
    )

    plt.legend()
    plt.xlabel("\$v_x\$")
    plt.ylabel("\$P(v_x)\$")
end


function plot_trajectories(pos, scatter=false, figname=nothing)
    num_particles, num_iterations = size(pos)

    trajectory_fig, trajectory_ax = plt.subplots()
    for i in 1:num_particles
        if scatter
            trajectory_ax.scatter(real(pos[i, 1:100:end]), imag(pos[i, 1:100:end]), label=string("Particle", i), s=1, c="#aaaaaa")
        else
            trajectory_ax.plot(real(pos[i, 1:100:end]), imag(pos[i, 1:100:end]), label=string("Particle", i))
        end
    end

    circ=plt.Circle((0, 0), radius=radius, fill=false)
    plt.gca().add_artist(circ)
    plt.gca().set_aspect("equal")
    #trajectory_ax.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(string(num_particles, " particles, T = ", total_time))
    if !isnothing(figname)
        plt.savefig(figname)
    end
    trajectory_fig.show()
end


function plot_energy_distribution(energies, figname=nothing)
    plt.subplots()
    num_particles, num_iterations = size(energies)
    for i in 1:num_particles
        plt.plot(0:dt:total_time, energies[i, :])
    end
    plt.xlabel("Time")
    plt.ylabel("Energy")
    if !isnothing(figname)
        plt.savefig(figname)
    end
    plt.show()
end
