# AttenuatedTotalReflectance.jl

A Julia package for simulating Attenuated Total Reflection (ATR) and the Transfer Matrix Method (TMM) in multilayer thin films. This module provides as a convient, flexible and fast method of preforming attenuatted total internal reflection simulations.

## Quick Start

```julia
using AttenuatedTotalReflectance
using Plots

# Define a stack representing a silver thin film on a sapphire prism.
mats = [
    (material = "Sapphire", n = 1.7659, k = 0.0, thickness = 5e-6),
    (material="Silver", n=0.056206, k=4.2776, thickness=45e-9),
    (material="Air", n=1.0, k=0.0, thickness=5e-6),
]
# Compile the simulation stack
stack = material_stack(mats)

# Run an angular sweep
angles_deg = collect(30.0:0.001:50.0)
angles = deg2rad.(angles_deg)
R, T, fx, fy, fz, fp = angular_ATR(stack, angles, 633.0e-9, 100e-9, 2)

# Plot the result

fontsize = 18
p1 = plot(angles_deg, 
    [R, T],
    xlabel = "Input Angle [deg]",
    ylabel = "Intensity [au]",
    title = "Surface Plamson Resonance of Silver Thin Film at $wavelength_nm nm",
    label=["Reflected" "Transmited"],
    lc=[:red :green],
    lw=3,
    xtickfontsize=fontsize, ytickfontsize=fontsize, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontsize=fontsize, titlefontsize = fontsize + 6,
    dpi = 300)

p2 = plot(angles_deg,
    [fx fz fp],
    xlabel = "Input Angle [deg]",
    ylabel = "Electric Field Intensity [au]",
    title = "Surface Plamson Resonance Electric Field Enhancement \n at $100 nm from the Silver Surface",
    label=["x-component" "z-component" "total"],
    lc= [:blue :red :purple],
    lw=3,
    xtickfontsize=fontsize, ytickfontsize=fontsize, xguidefontsize=fontsize, yguidefontsize=fontsize, legendfontsize=fontsize, titlefontsize = fontsize + 6,
    dpi = 300)

p = plot(p1, p2,  layout=(1,2),
                size = (3000,1500),
                left_margin = 20Plots.mm, 
                right_margin = 5Plots.mm,
                top_margin=15Plots.mm,
                bottom_margin = 20Plots.mm,
                )

display(p)
```