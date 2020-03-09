"""This program investigates the problems conserning one particle"""

include("utils.jl")
#########
# Setup #
#########

num_particles = 2
radius = 10
KK = 5
total_time = 100
dt = 0.0001
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
