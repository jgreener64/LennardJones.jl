# Molecular dynamics simulation of a simple Lennard Jones system
# See https://www.saylor.org/site/wp-content/uploads/2011/06/MA221-6.1.pdf for
#   integration algorithm - used shorter second version
# See http://phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf for equations

module LennardJones

using ProgressMeter
using Plots

export
    Simulation,
    kinetic_energy,
    potential_energy,
    energy,
    temperature,
    coord_array,
    simulate!,
    plot_energy,
    plot_coords

const sqdist_cutoff = 20.0^2

mutable struct Coordinates
    x::Float64
    y::Float64
    z::Float64
end

struct Atom
    id::Int32
    coords::Coordinates
end

mutable struct Velocity
    x::Float64
    y::Float64
    z::Float64
end

struct Universe
    atoms::Vector{Atom}
    velocities::Vector{Velocity}
    box_size::Float64
end

mutable struct Simulation
    universe::Universe
    timestep::Float64
    n_steps::Int32
    steps_made::Int32
    pes::Vector{Float64}
    kes::Vector{Float64}
    energies::Vector{Float64}
    temps::Vector{Float64}
end

function Atom(id, box_size::Real)
    return Atom(id, Coordinates(rand()*box_size, rand()*box_size, rand()*box_size))
end

function Velocity(starting_velocity)
    θ = rand()*π
    ϕ = rand()*2π
    r = rand()*starting_velocity
    return Velocity(r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ))
end

function Simulation(n_atoms, box_size, starting_velocity, timestep, n_steps)
    a = [Atom(i, box_size) for i in 1:n_atoms]
    v = [Velocity(starting_velocity) for i in 1:n_atoms]
    u = Universe(a, v, box_size)
    return Simulation(u, timestep, n_steps, 0, [], [], [], [])
end

function Base.show(io::IO, s::Simulation)
    print(io, "Simulation of ", length(s.universe.atoms),
            " atoms in a box of size ", Int(round(s.universe.box_size, 0)),
            ": ", s.steps_made, "/", s.n_steps, " steps made (",
            round(100*s.steps_made/s.n_steps, 1), " %)")
end

mutable struct Acceleration
    x::Float64
    y::Float64
    z::Float64
end

function update_coordinates!(atoms, velocities, accels, timestep, box_size)
    for (i, a) in enumerate(atoms)
        a.coords.x += velocities[i].x*timestep + 0.5*accels[i].x*timestep^2
        a.coords.y += velocities[i].y*timestep + 0.5*accels[i].y*timestep^2
        a.coords.z += velocities[i].z*timestep + 0.5*accels[i].z*timestep^2
        a.coords.x = bound!(a.coords.x, box_size)
        a.coords.y = bound!(a.coords.y, box_size)
        a.coords.z = bound!(a.coords.z, box_size)
    end
    return atoms
end

function bound!(coord, box_size)
    if coord <= 0.0
        coord += box_size
        return bound!(coord, box_size)
    elseif coord > box_size
        coord -= box_size
        return bound!(coord, box_size)
    end
    return coord
end

function update_velocities!(velocities, accels_t, accels_t_dt, timestep)
    for (i, v) in enumerate(velocities)
        v.x += 0.5*(accels_t[i].x+accels_t_dt[i].x)*timestep
        v.y += 0.5*(accels_t[i].y+accels_t_dt[i].y)*timestep
        v.z += 0.5*(accels_t[i].z+accels_t_dt[i].z)*timestep
    end
    return velocities
end

function sqdist(a1, a2, box_size)
    return bounddist(a1.coords.x, a2.coords.x, box_size)^2 +
        bounddist(a1.coords.y, a2.coords.y, box_size)^2 +
        bounddist(a1.coords.z, a2.coords.z, box_size)^2
end

function bounddist(x1, x2, box_size)
    if x1 < x2
        return min(x2-x1, x1-x2+box_size)
    else
        return min(x1-x2, x2-x1+box_size)
    end
end

function lennardjones_potential(a1, a2, box_size)
    r2 = sqdist(a1, a2, box_size)
    return 4.0*(r2^-6 - r2^-3)
end

function force(a1, a2, box_size)
    r2 = sqdist(a1, a2, box_size)
    if r2 < sqdist_cutoff
        f = 48*(r2^-7 - 0.5*r2^-4)
        return f*(a1.coords.x-a2.coords.x), f*(a1.coords.y-a2.coords.y), f*(a1.coords.z-a2.coords.z)
    else
        return 0.0, 0.0, 0.0
    end
end

function update_accelerations!(accels, atoms, box_size)
    update_x = 0.0
    update_y = 0.0
    update_z = 0.0
    for (i, a1) in enumerate(atoms)
        accel_x = 0.0
        accel_y = 0.0
        accel_z = 0.0
        for a2 in atoms
            if a1.id != a2.id
                update_x, update_y, update_z = force(a1, a2, box_size)
                accel_x += update_x
                accel_y += update_y
                accel_z += update_z
            end
        end
        accels[i].x = accel_x
        accels[i].y = accel_y
        accels[i].z = accel_z
    end
    return accels
end

function kinetic_energy(velocities)
    ke = 0.0
    for v in velocities
        ke += v.x^2 + v.y^2 + v.z^2
    end
    return 0.5*ke
end

function potential_energy(universe)
    pe = 0.0
    for (i, a1) in enumerate(universe.atoms)
        for j in 1:i-1
            pe += lennardjones_potential(a1, universe.atoms[j], universe.box_size)
        end
    end
    return pe
end

energy(universe) = kinetic_energy(universe.velocities) + potential_energy(universe)

temperature(ke, n_atoms) = (2*ke)/(3*n_atoms)

coord_array(universe) = vcat([[a.coords.x a.coords.y a.coords.z] for a in universe.atoms]...)

empty_accelerations(n_atoms) = [Acceleration(0.0, 0.0, 0.0) for i in 1:n_atoms]

simulate!(s) = simulate!(s, s.n_steps-s.steps_made)

function simulate!(s, n_steps)
    n_atoms = length(s.universe.atoms)
    a_t = update_accelerations!(empty_accelerations(n_atoms), s.universe.atoms, s.universe.box_size)
    a_t_dt = empty_accelerations(n_atoms)
    @showprogress for i in 1:n_steps
        update_coordinates!(s.universe.atoms, s.universe.velocities, a_t, s.timestep, s.universe.box_size)
        update_accelerations!(a_t_dt, s.universe.atoms, s.universe.box_size)
        update_velocities!(s.universe.velocities, a_t, a_t_dt, s.timestep)
        if i % 100 == 0
            pe = potential_energy(s.universe)
            ke = kinetic_energy(s.universe.velocities)
            push!(s.pes, pe)
            push!(s.kes, ke)
            push!(s.energies, pe+ke)
            push!(s.temps, temperature(ke, n_atoms))
        end
        a_t = a_t_dt
        s.steps_made += 1
        #i%10000==0 && println(s.universe.atoms[1].coords, s.universe.velocities[1], a_t[1])
    end
    return s
end

function plot_energy(simulation; kwargs...)
    v = simulation.energies[1:10:end]
    return plot(collect(1:length(v)), v; kwargs...)
end

function plot_coords(coords...; kwargs...)
    xs = hcat((i[:,1] for i in coords)...)
    ys = hcat((i[:,2] for i in coords)...)
    zs = hcat((i[:,3] for i in coords)...)
    return scatter(xs, ys, zs; kwargs...)
end

end # LennardJones
