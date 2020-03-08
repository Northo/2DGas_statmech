"""This program investigates the problems conserning one particle"""

include("verlet.jl")
#########
# Setup #
#########

num_particles = 1
radius = 10
KK = 5
total_time = 200

#dt = 0.005
#num_iterations = Int(1E5)


################
# Calculations #
################
dt_list = [0.005, 0.01, 0.05, 0.1, 0.6]

positions = []
velocities = []

initial_vals = initial_values(num_particles, radius)

for (i,dt) in enumerate(dt_list)
    num_iterations = floor(Int, total_time / dt)
    pos, vel = billiard(
        num_particles,
        num_iterations,
        radius,
        dt,
        KK,
        initial_value=initial_vals
    )

    push!(positions, pos)
    push!(velocities, vel)
end

###################
# Writing results #
###################

############
# Plotting #
############
println("Loading PyPlot...")
using PyPlot
println("PyPlot loaded.")

for (i, dt) in enumerate(dt_list)
    t = 0:dt:total_time
    plt.plot(positions[i][1, 1, :], positions[i][2, 1, :], label=dt)
end

circ=plt.Circle((0, 0), radius=radius, fill=false)
gca().add_artist(circ)
gca().set_aspect("equal")
legend()
#show()

pos = positions[1]
vel = velocities[1]

engy = energy.(pos[1, 1, :], pos[2, 1, :], vel[1, 1, :], vel[2, 1, :], radius, KK)
pot = potential.(pos[1, 1, :], pos[2, 1, :], radius, KK)

plt.figure()
plt.plot(engy)
plt.show()

# for (i, dt) in enumerate(dt_list)
#     t = 0:dt:total_time
#     pos = positions[:, 1:length(t), i]
#     distance = sqrt.(pos[1,:].^2 + pos[2,:].^2)
#     vs = diff(pos, dims=2)
#     t = t[1:end-1]  # Because of diff
#     E = sqrt.(vs[1, :].^2 + vs[2, :].^2)
#     potential = KK * ifelse.(distance .> radius, 1.0, 0.0).*(distance .- radius).^2/2
#     println(size(potential))
#     E += potential[1:end-1]
#     relative_E = E./E[1]
#     plt.plot(t, relative_E, label=dt)
# end
# plt.legend()
# plt.show()
