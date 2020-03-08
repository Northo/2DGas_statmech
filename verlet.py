import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def billiard(numberOfParticles, numberOfIterations, radius, dt, KK):
    pos = np.zeros((2, numberOfParticles, numberOfIterations+1))
    pos[0,:,0] = np.random.uniform(-radius, radius, numberOfParticles)
    pos[1,:,0] = np.random.uniform(-radius, radius, numberOfParticles)

    angles = np.random.uniform(-np.pi, np.pi)
    vx = np.cos(angles)
    vy = np.sin(angles)

    for i in range(numberOfIterations):
        distance = np.sqrt(pos[0,:,i]**2 + pos[1,:,i]**2)
        acceleration = np.zeros((2, numberOfParticles))

        #test if we are outside the circle and calculate force from wall
        acceleration[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius, i] / distance[distance > radius]
        acceleration[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius, i] / distance[distance > radius]

        #update positions
        pos[0, :, i+1] =  pos[0, :, i] + vx * dt + (acceleration[0, :] * dt**2)/2
        pos[1, :, i+1] =  pos[1, :, i] + vy * dt + (acceleration[1, :] * dt**2)/2

        distance = np.sqrt(pos[0,:,i+1]**2 + pos[1,:,i+1]**2)
        acceleration2 = np.zeros((2, numberOfParticles))
        #test if we are outside the circle and calculate force from wall
        acceleration2[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius, i+1] / distance[distance > radius]
        acceleration2[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius, i+1] / distance[distance > radius]

        #update velocities
        vx += (acceleration[0,:] + acceleration2[0,:])/2 * dt
        vy += (acceleration[1,:] + acceleration2[1,:])/2 * dt

    return pos


pos = billiard(1, 100000, 10, 0.005, 5)

plt.plot(pos[0,0,0:-1:100], pos[1,0,0:-1:100])
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.show()
