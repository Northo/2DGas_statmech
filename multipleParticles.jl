"""This program investigates the problems conserning one particle"""

include("utils.jl")
using Statistics  # Used for mean
using Distributions  # fit
#########
# Setup #
#########

num_particles = 10
radius = 10
KK = 15
total_time = 200
dt = 0.001
num_iterations = floor(Int, total_time / dt)

epsilon = 0.5
initial_vals = initial_values(num_particles, radius)
# initial_vals = [
#     [5 -7;
#      5 -7],
#     [0, -pi]
# ]

println("Inital_values: ", initial_vals)

################
# Calculations #
################

println("Calculating trajectories...")
pos, vel = billiard(
    num_particles,
    num_iterations,
    radius,
    dt,
    KK,
    initial_value=initial_vals,
    epsilon=epsilon,
)

# Numerical validation
println("Calculating energy for validation...")
x, y = pos[1, :, :], pos[2, :, :]  # x and y still has dimensions for particles and iterations!
vel_x, vel_y = vel[1, :, :], pos[2, :, :]

engy = energy.(x, y, vel_x, vel_y, radius, KK)  # Not including interactions between particles
interaction_energy = zero(engy)
for i in 1:num_particles, j in 1:num_particles
    if i==j
        continue
    end
    interaction_energy[i, :] += lennard_jones_potential.(
        x[i, :],
        y[i, :],
        x[j, :],
        y[j, :],
        epsilon,
    )
end

total_energy = sum(engy + interaction_energy, dims=2)
relative_energy_error = (total_energy .- total_energy[1]) ./ total_energy[1]


## Estimating velocity distribution ##
velocity_distribution = collect(Iterators.flatten(vel_x))
println(typeof(velocity_distribution))
fit_velocity_distribution = fit(Normal, velocity_distribution)
μ, σ = fit_velocity_distribution.μ, fit_velocity_distribution.σ
#println("Fitted: ", μ, σ)
println("Mean of vel_x: ", mean(vel_x))
println("Mean of vel_x^2: ", mean(x->x^2, vel_x))
println("Mean of x: ", mean(x))
println("RMS of x: ", sqrt(mean(x->x^2, x)))

###################
# Writing results #
###################

############
# Plotting #
############
println("Loading PyPlot...")
using PyPlot
println("PyPlot loaded.")


## Trajectory ##
trajectory_fig, trajectory_ax = plt.subplots()
for i in 1:num_particles
    trajectory_ax.plot(pos[1, i, 1:100:end], pos[2, i, 1:100:end], label=string("Particle", i))
end

circ=plt.Circle((0, 0), radius=radius, fill=false)
plt.gca().add_artist(circ)
#plt.gca().set_aspect("equal")
trajectory_ax.legend()
trajectory_fig.show()

## Energy validation ##
engy_fig, engy_ax = plt.subplots()
engy_ax.plot(relative_energy_error)
engy_fig.show()

vel_distribution_fig, vel_distribution_ax = plt.subplots()
v_min, v_max = minimum(vel_x), maximum(vel_x)
v_list = range(v_min, v_max, length=100)

# Boltzman dist.
kbT = mean(x->x^2, vel_x)  # From equipartition thm.
boltzman_std = sqrt(kbT)
boltzman_dist = Normal(0, boltzman_std)

plt.hist(collect(Iterators.flatten(vel_x)), bins=40, density=true, label="Simulated distribution")
plt.plot(v_list, Distributions.pdf.(fit_velocity_distribution, v_list), label="Fitted curve")
plt.plot(v_list, Distributions.pdf.(boltzman_dist, v_list), label="Boltzman distribution")
plt.legend()
plt.show()
