"""This program investigates the problems conserning one particle"""

include("verlet.jl")
#########
# Setup #
#########

num_particles = 1
radius = 10
KK = 5
total_time = 100

#dt = 0.005
#num_iterations = Int(1E5)


################
# Calculations #
################
dt = 0.01
num_iterations = floor(Int, total_time/ dt)
@time pos1 = billiard(num_particles, num_iterations, radius, dt, KK)
@time pos2 = billiard2(num_particles, num_iterations, radius, dt, KK)


# dt_list = [0.005, 0.01, 0.05, 0.1]
# positions = zeros(
#     2,
#     Int(total_time/minimum(dt_list))+1,
#     length(dt_list)
# )

# for (i,dt) in enumerate(dt_list)
#     num_iterations = floor(Int, total_time / dt)
#     pos = billiard(
#         num_particles,
#         num_iterations,
#         radius,
#         dt,
#         KK,
#     )
#     num = size(pos)[3]
#     positions[:, 1:num, i] = pos[:, 1, :]
# end

###################
# Writing results #
###################

############
# Plotting #
############
# println("Loading PyPlot...")
# using PyPlot
# println("PyPlot loaded.")
# for (i, dt) in enumerate(dt_list)
#     t = 0:dt:total_time
#     plt.plot(positions[2, 1:length(t), i], positions[1, 1:length(t), i], label=dt)
# end
# plt.legend()
# plt.show()

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
