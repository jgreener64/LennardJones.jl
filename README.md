# LennardJones.jl

[![Build Status](https://travis-ci.org/jgreener64/LennardJones.jl.svg?branch=master)](https://travis-ci.org/jgreener64/LennardJones.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/fwf42unlh61ctt68?svg=true)](https://ci.appveyor.com/project/jgreener64/lennardjones-jl)
[![codecov.io](http://codecov.io/github/jgreener64/LennardJones.jl/coverage.svg?branch=master)](http://codecov.io/github/jgreener64/LennardJones.jl?branch=master)

Molecular dynamics simulation of a simple Lennard Jones system.

Status: work in progress.

## Example usage

```julia
using LennardJones
using Plots

# Reduced units
n_atoms = 10
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
plot_coords(coords_start, coords_end, labels=["start", "end"])
```
