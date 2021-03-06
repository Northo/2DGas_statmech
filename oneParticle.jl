"""This program investigates the problems conserning one particle"""

include("utils.jl")
#########
# Setup #
#########

num_particles = 1
radius = 10
KK = 10
total_time = 1000

#dt = 0.005
#num_iterations = Int(1E5)


################
# Calculations #
################
dt_list = [
#    0.005,
#    0.01,
    0.1,
#    0.2,
#    0.3,
#    0.5,
]

positions = []
velocities = []

initial_positions = generate_initial_positions(num_particles, radius)
initial_velocities = generate_initial_velocities(num_particles)

for (i,dt) in enumerate(dt_list)
    num_iterations = floor(Int, total_time / dt)
    pos, vel = billiard(
        num_particles,
        num_iterations,
        radius,
        dt,
        KK,
        initial_positions=initial_positions,
        initial_velocities=initial_velocities
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

plt.style.use("bmh")

for (i, dt) in enumerate(dt_list)
    t = 0:dt:total_time
    plt.plot(real(positions[i][1, :]), imag(positions[i][1, :]), label=dt)
end

circ=plt.Circle((0, 0), radius=radius, fill=false)
gca().add_artist(circ)
#gca().set_aspect("equal")
#plt.title("Trajectories")
#legend()
#plt.savefig("media/trajectory_one_particle.pdf")
show()

#energies = [energy.(positions[i][1, :], velocities[i][1, :], radius, KK) for i in eachindex(dt_list)]

plt.figure()
markers = ["-", "--", ":", "-.", "."]
for (i, dt) in enumerate(dt_list)
    engy = energy.(positions[i][1, :], velocities[i][1, :], radius, KK)
    pot = potential.(positions[i][1, :], radius, KK)
    kin = 0.5*abs2.(velocities[i][1, :])
    plt.plot(0:dt:total_time, relative_error(engy), label=dt)
end
plt.legend()
plt.ylim(-0.04, 0.04)
plt.title("Energy")
plt.ylabel("Relative error in energy")
plt.xlabel("Time")
#plt.savefig("media/errors_one_particle.pdf", dpi=300)
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
