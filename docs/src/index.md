# AttenuatedTotalReflectance.jl

A Julia package for simulating Attenuated Total Reflection (ATR) through the Transfer Matrix Method (TMM) in multilayer thin films. This module provides as a convient, flexible and fast method of preforming attenuated total reflection simulations.

## Quick Start
```julia
using AttenuatedTotalReflectance
using Plots

# 1. Setup Simulation Parameters
λ_nm = 633.0
λ_m = λ_nm * 1e-9

# Define materials:
mats = [
    (material="Glass", n=1.51, k=0.0, thickness=1e-5),
    (material="Silver", n=0.056206, k=4.2776, thickness=47e-9),
    (material="Air", n=1.0, k=0.0, thickness=5e-6),
]

stack = material_stack(mats)
println("System to Simulate: $(stack)")
# 2. Run Angular Sweep (35° to 55°)
angles_deg = collect(35.0:0.001:55.0)
angles_rad = deg2rad.(angles_deg)

# Calculate Reflection, Transmittance, and Field Enhancement 100 nm from the Silver Surface
R, T, fx, fy, fz, fp = angular_ATR(stack, angles_rad, λ_m, 100e-9, 2)

# 3. Visualize the Results
p1 = plot(angles_deg, [R, T],
    label=["Reflectivity" "Transmittance"],
    ylabel="Intensity [au]", lc=[:red :green],grid=true, legend=:bottomright, lw = 2, title="Reflection & Transmission Curve\n at $λ_nm nm")

p2 = plot(angles_deg, fp,
    label="Field Enhancement",
    ylabel="|E/E₀|²", lc=:purple, lw=2,grid=true, title="Electric Field Enhancement\nat 100 nm from $(stack[2].material) Surface")

plot(p1, p2, layout=(1, 2), size=(900, 400), xlabel="Angle (deg)",gridstyle = :dash, gridlinewidth = 3, margin=5Plots.mm)
```
![Example simulation of silver film on a glass prism. Characteristic Surface plasmon response is observed](/assets/Quickstart.png)