# API Reference

This page provides a detailed list of all functions and types exported by `AttenuatedTotalReflectance.jl`.

## [Index](@id main_index)
```@index

[Core Simulation](@id core_sim)

These are the high-level functions used to generate reflectivity and field enhancement data.

angular_ATR
wavelength_ATR
compute_transfer_coefficents

[System Configuration](@id system_config)

Use these to define your material layers and the geometry of your simulation.

layer
material_stack
target_layer

[Fresnel Coefficients](@id fresnel_coeffs)

The underlying equations for reflection and transmission at single interfaces.

refl_coeff_P
refl_coeff_S
trans_coeff_P
trans_coeff_S

[Optical Utilities](@id optical_utils)

Fundamental calculations for refractive indices, angles, and wavevectors.

complex_n
snells_law
epsilon_to_nk
nk_to_epsilon