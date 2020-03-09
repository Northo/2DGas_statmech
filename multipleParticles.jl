"""This program investigates the problems conserning one particle"""

include("utils.jl")
using Statistics  # Used for mean
using Distributions  # fit
using DelimitedFiles
#########
# Setup #
#########
SAVE_DATA = false
DATA_DIR = "datadir/"
FIG_DIR = "media/"
FIT = false

num_particles = 2
num_particles_word = "five"
radius = 10
KK = 15
total_time = 400
dt = 0.005
num_iterations = floor(Int, total_time / dt)

epsilon = 0.5

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
    epsilon=epsilon,
)

# Numerical validation
println("Calculating energy for validation...")

engy = energy.(pos, vel, radius, KK)  # Not including interactions between particles
interaction_energy = get_interaction_potential(pos, epsilon)

total_energy = sum(engy + interaction_energy, dims=2)
relative_energy_error = (total_energy .- total_energy[1]) ./ total_energy[1]

x = real.(pos)
vel_x = real.(vel)

## Estimating velocity distribution ##
if FIT
    velocity_distribution = collect(Iterators.flatten(vel_x))
    try
        fit_velocity_distribution = fit(Normal, velocity_distribution)
    catch e
        println("Was unable to fit curve, error: ", e)
        fit_veloctiy_distribution = Normal(0, 1)
    end
    μ, σ = fit_velocity_distribution.μ, fit_velocity_distribution.σ
end

println("Mean of vel_x: ", mean(vel_x))
println("Mean of vel_x^2: ", mean(x->x^2, vel_x))
println("Mean of x: ", mean(x))
println("RMS of x: ", sqrt(mean(x->x^2, x)))

###################
# Writing results #
###################
if SAVE_DATA
    filename = string(DATA_DIR, "pos_vel_num_particles_", num_particles, "_num_iterations_", num_iterations, ".txt")
    println("Writing data to file...")
    open(filename, "w") do file
        writedlm(file, [pos, vel])
    end
end

############
# Plotting #
############
println("Loading PyPlot...")
using PyPlot
println("PyPlot loaded.")


## Trajectory ##
trajectory_fig, trajectory_ax = plt.subplots()
for i in 1:num_particles
    #trajectory_ax.scatter(real(pos[i, 1:100:end]), imag(pos[i, 1:100:end]), label=string("Particle", i), s=1, c="#aaaaaa")
    trajectory_ax.plot(real(pos[i, 1:100:end]), imag(pos[i, 1:100:end]), label=string("Particle", i))
end

circ=plt.Circle((0, 0), radius=radius, fill=false)
plt.gca().add_artist(circ)
plt.gca().set_aspect("equal")
#trajectory_ax.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title(string(num_particles, " particles, T = ", total_time))
plt.savefig(string(FIG_DIR, "trajectory_", num_particles_word ,"_particles.pdf"))
trajectory_fig.show()

## Energy distribution ##
plt.subplots()
for i in 1:num_particles
    plt.plot(engy[i, :])
end
plt.show()

# ## Energy validation ##
# engy_fig, engy_ax = plt.subplots()
# engy_ax.plot(relative_energy_error)
# engy_fig.show()

#plot_velocity_distributions(vel)
