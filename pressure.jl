using Statistics  # Used for mean
using Printf
using PyPlot

include("utils.jl")

#########
# Setup #
#########
SAVE_DATA = false
DATA_DIR = "datadir/"
FIG_DIR = "media/"

VALIDATION = false

KK = 15
total_time = 200
dt = 0.003
num_iterations = floor(Int, total_time / dt)

epsilon = 0.5
#initial_velocities = zeros(num_particles)
#initial_velocities[1] = sqrt(num_particles*2)

num_particles_list = [5, 10, 15, 5]
radius_list = [6, 10, 13, 10]


## Helper function ##
function my_abs(pos::Array{Complex}, radius)
    """Calculates distances where distance > radius, ie
    (r_i - R) * step(r_i - R)"""

    absolutes = zeros(Float64, size(pos))
    for i in eachindex(pos)
        if abs(pos[i]) <= radius
            absolutes[i] = 0
        else
            absolutes[i] = abs(pos[i]) - radius
        end
    end
    return absolutes
end


################
# Calculations #
################
energy_errors = []
for (i, num_particles, radius) in Iterators.zip(1:length(radius_list), num_particles_list, radius_list)
    @printf(" - N = %i, R = %.2f\n", num_particles, radius)

    println("Calculating trajectories...")
    pos, vel = billiard(
        num_particles,
        num_iterations,
        radius,
        dt,
        KK,
        epsilon=epsilon,
        initial_positions=safe_initial_positions(num_particles, radius),
        #    initial_velocities=initial_velocities
    )

    # Validation #
    if VALIDATION
        engy = energy.(pos, vel, radius, KK)  # Not including interactions between particles
        interaction_energy = get_interaction_potential(pos, epsilon)
        total_energy = sum(engy + interaction_energy, dims=1)
        relative_energy_error = (total_energy .- total_energy[1]) ./ total_energy[1]
        push!(energy_errors, relative_energy_error)
    end

    println("Calculating pressure...")
    force = KK * sum(my_abs(pos, radius), dims=1)
    pressure = force/ (2*pi*radius)

    # pAN is the value p*A/N, should be constant = k_B T
    pAN = mean(pressure) * pi * radius^2 / num_particles
    println("pAN: ", pAN)

    kbT = 1
    pressure_theory = num_particles * kbT / (pi*radius^2)
    #    plt.figure()
    plt.subplot(2, 2, i)
    plt.title(string("N= ", num_particles, ", R = ", radius))
    plt.plot(pressure[:])
    plt.axhline(pressure_theory)
    # plt.annotate(
    #     @sprintf("\$E[P] A / N\$ from simulation : %.3f", pAN),
    #     (0.25, 0.8),
    #     xycoords="figure fraction",
    #     bbox = Dict("boxstyle"=>"round", "ec"=>(1., 0.5, 0.5), "fc"=>(1., 0.8, 0.8))
    # )
    #plt.show()

    println()
end
plt.tight_layout()
#plt.savefig("media/pressures.pdf")

if VALIDATION
    println("plotting errors")
    using PyPlot

    for (i, energy_error) in enumerate(energy_errors)
        plt.plot(energy_error[1:20:end], label=i)
    end
    plt.legend()
    plt.show()
end
