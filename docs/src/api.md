# API Reference

This page provides a detailed list of all functions and types exported by `AttenuatedTotalReflectance.jl`.

## [Index](@id main_index)
```@index

## Core Simulation
These are the high-level functions used to generate reflectivity and field enhancement data.

angular_ATR
wavelength_ATR
compute_transfer_coefficients

## Simulation Configuration

Use these to define your material layers and the geometry of your simulation.

layer
material_stack
target_layer

## Fresnel Coefficients

The underlying equations for reflection and transmission at single interfaces.

refl_coeff_P
refl_coeff_S
trans_coeff_P
trans_coeff_S

## Optical Utilities

Fundamental calculations for refractive indices, angles, and wavevectors.

complex_n
snells_law
epsilon_to_nk
nk_to_epsilon