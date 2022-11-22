import numpy as np
import matplotlib.pyplot as plt
"""
Assume a system of N point particles, indexed by i = 1,...,N
Each particle has:
mass
position [x,y,z]
velocity [x,y,z]
Each particle experiences gravitational attraction to each other
"""

def getAcceleration(pos, mass, G, softening):
	"""
	Calculate the acceleration of the particle according to Newton's Law
	pos is an N x 3 matrix position
	mass is an N x 1 vector of masses
	G is a constant
	softening is the softening length
	a is N x 3 matrix of acceleration
	"""

	#positions r = [x,y,z] for all particles
	x = pos[:,0:1] #Get column 0
	y = pos[:,1:2] #Get column 1
	z = pos[:,2:3] #Get column 2

	#Matrix to store pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	#Matrix to store 1/r^3 for all particles
	inv_r3 = (dx**2 + dy**2 + softening**2)
	inv_r3[inv_r3 > 0] = inv_r3[inv_r3 > 0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass

	#Pack together acceleration components
	a = np.hstack((ax,ay, az))

	return a


def getEnergy(pos, vel, mass, G):
	"""
	Get the kinetic energy and potential of simulation
	pos is N x 3 matrix of positions
	vel is N x 3 matrix of velocities
	G is a constant
	KE is the kinetic energy of the system
	PE is the potential energy of the system
	"""

	#Kinetic Energy:
	KE = 0.5 * np.sum(np.sum(mass * vel ** 2))

	#Potential energy 

	#Positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	#Matrix that stores all the pairwise separations
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	#Matrix that stores 1/r for all particles
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r > 0] = 1.0/inv_r[inv_r > 0]

	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r, 1)))

	return KE, PE

def main():
	"""
	N-body simulation
	"""
	N = 100 #Number of particles
	t = 0 #Current time of simulation
	tEnd = 10.0 #Time at which the simulation ends
	dt = 0.01 #Time step
	softening = 0.1 #Softening length
	G = 1.0 #Newton's Gravitational Constant
	plotRealTime = True

	#Generate initial conditions
	np.random.seed(21)

	mass = 20.0 * np.ones((N,1))/N #Total mass of particles
	pos = np.random.randn(N,3)
	vel = np.random.randn(N,3)




















