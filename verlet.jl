function billiard(numberOfParticles, numberOfIterations, radius, dt, KK)
    pos = zeros(2, numberOfParticles, numberOfIterations+1)
    pos[1,:,1] = rand((-radius:radius), numberOfParticles)
    pos[2,:,1] = rand((-radius:radius), numberOfParticles)

    angles = rand((-pi:pi), numberOfParticles)
    vx = cos.(angles)
    vy = sin.(angles)

    for i in 1:numberOfIterations
        distance = sqrt.(pos[1,:,i].^2 + pos[2,:,i].^2)
        acceleration = zeros((2, numberOfParticles))

        acceleration[1, distance .> radius] = -KK * (distance[distance .> radius] .- radius) .* pos[1, distance .> radius, i] ./ distance[distance .> radius]
        acceleration[2, distance .> radius] = -KK * (distance[distance .> radius] .- radius) .* pos[2, distance .> radius, i] ./ distance[distance .> radius]

        pos[1, :, i+1] = pos[1, :, i] .+ vx*dt + (acceleration[1, :]*dt^2)/2
        pos[2, :, i+1] = pos[2, :, i] .+ vy*dt + (acceleration[2, :]*dt^2)/2

        distance = sqrt.(pos[1, :, i+1].^2 + pos[2, :, i+1].^2)
        acceleration2 = zeros((2, numberOfParticles))
        acceleration2[1, distance .> radius] = -KK * (distance[distance .> radius] .- radius) .* pos[1, distance .> radius, i+1] ./ distance[distance .> radius]
        acceleration2[2, distance .> radius] = -KK * (distance[distance .> radius] .- radius) .* pos[2, distance .> radius, i+1] ./ distance[distance .> radius]

        vx += (acceleration[1,:] + acceleration2[1,:])/2 * dt
        vy += (acceleration[2,:] + acceleration2[2,:])/2 * dt
    end

    return pos
end


function billiard2(numberOfParticles, numberOfIterations, radius, dt, KK)
    pos = zeros(2, numberOfParticles, numberOfIterations+1)
    pos[1,:,1] = rand((-radius:radius), numberOfParticles)
    pos[2,:,1] = rand((-radius:radius), numberOfParticles)
    velocity = Array{Float64, 3}(undef, 2, numberOfParticles, numberOfIterations+1)

    angles = rand((-pi:pi), numberOfParticles)
    vx = cos.(angles)
    vy = sin.(angles)
    velocity[:, :, 1] = transpose([vx vy])

    for i in 1:numberOfIterations, j in 1:numberOfParticles
        distance = sqrt(pos[1,j,i]^2 + pos[2,j,i]^2)
        if distance > radius
            acceleration = -KK * ((distance - radius)/distance) * pos[:, j, i]
        else
            acceleration = [0, 0]
        end

        pos[:, j, i+1] = pos[:, j, i] + velocity[:, j, i]*dt + acceleration*dt^2 / 2

        distance = sqrt(pos[1, j, i+1]^2 + pos[2, j, i+1]^2)
        if distance > radius
            acceleration2 = -KK * ((distance - radius)/distance) * pos[:, j, i+1]
        else
            acceleration2 = [0, 0]
        end

        velocity[:, j, i+1] = velocity[:, j, i] + (acceleration + acceleration2)/2 * dt
    end

    return pos, velocity
end

function potential(x, y, radius, KK)
    distance = sqrt(x^2 + y^2)
    if distance <= radius
        return 0
    end
    return (distance - radius)^2 * KK / 2
end

function energy(x, y, vel_x, vel_y, radius, KK)
    pot = potential(x, y, radius, KK)
    total_vel_square = vel_x^2 + vel_y^2

    return 0.5*total_vel_square + pot
end
