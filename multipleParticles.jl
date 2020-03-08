"""This program investigates the problems conserning one particle"""

include("verlet.jl")
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

###################
# Writing results #
###################

############
# Plotting #
############
println("Loading PyPlot...")
using PyPlot
println("PyPlot loaded.")

for i in 1:num_particles
    plt.plot(pos[1, i, 1:100:end], pos[2, i, 1:100:end], label=string("Particle", i))
end

circ=plt.Circle((0, 0), radius=radius, fill=false)
plt.gca().add_artist(circ)
#plt.gca().set_aspect("equal")
plt.legend()
plt.show()
