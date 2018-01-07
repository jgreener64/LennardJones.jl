using LennardJones
using Base.Test

n_atoms = 5
box_size = 50.0
starting_velocity = 0.0
timestep = 0.01
n_steps = 1e6

s = Simulation(n_atoms, box_size, starting_velocity, timestep, n_steps)

simulate!(s)
