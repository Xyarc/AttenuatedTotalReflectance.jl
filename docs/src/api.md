 API Reference

This page provides a detailed list of all functions and types exported by `AttenuatedTotalReflectance.jl`.

## [Index](@id main_index)

```@index
```

---

## Core Simulation
These are the high-level functions used to generate reflectivity and field enhancement data.

```@docs
AttenuatedTotalReflectance.compute_transfer_coefficients
AttenuatedTotalReflectance.angular_ATR
AttenuatedTotalReflectance.wavelength_ATR
```

## Simulation Configuration
Use these to define your material layers and the geometry of your simulation.

```@docs
AttenuatedTotalReflectance.layer
AttenuatedTotalReflectance.material_stack
AttenuatedTotalReflectance.target_layer
```

## Fresnel Coefficients
The underlying equations for reflection and transmission at single interfaces.

```@docs
AttenuatedTotalReflectance.refl_coeff_P
AttenuatedTotalReflectance.refl_coeff_S
AttenuatedTotalReflectance.trans_coeff_P
AttenuatedTotalReflectance.trans_coeff_S
```

## Optical Utilities
Helper functions for Snell's law and transforming optical constants

```@docs
AttenuatedTotalReflectance.complex_n
AttenuatedTotalReflectance.snells_law 
AttenuatedTotalReflectance.epsilon_to_nk
AttenuatedTotalReflectance.nk_to_epsilon
```