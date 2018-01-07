# LennardJones

Molecular dynamics simulation of a simple Lennard Jones system.

Status: work in progress.

## Example usage

```julia
using LennardJones
using Plots

# Reduced units
n_atoms = 5
box_size = 50.0
starting_velocity = 0.0
timestep = 0.01
n_steps = 1e6

s = Simulation(n_atoms, box_size, starting_velocity, timestep, n_steps)

coords_start = coord_array(s.universe)
println("Starting PE: ", potential_energy(s.universe))
println("Starting KE: ", kinetic_energy(s.universe.velocities))
println("Starting T:  ", temperature(kinetic_energy(s.universe.velocities), n_atoms))

simulate!(s)

coords_end = coord_array(s.universe)
println("Ending PE:   ", potential_energy(s.universe))
println("Ending KE:   ", kinetic_energy(s.universe.velocities))
println("Ending T:    ", temperature(kinetic_energy(s.universe.velocities), n_atoms))

# View starting and ending coordinates
scatter([coords_start[:,1] coords_end[:,1]], [coords_start[:,2] coords_end[:,2]])
```
