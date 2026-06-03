# AttenuatedTotalReflectance

[![Build Status](https://github.com/Xyarc/AttenuatedTotalReflectance.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Xyarc/AttenuatedTotalReflectance.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/Xyarc/AttenuatedTotalReflectance.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Xyarc/AttenuatedTotalReflectance.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Xyarc.github.io/AttenuatedTotalReflectance.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Xyarc.github.io/AttenuatedTotalReflectance.jl/dev/)


`AttenuatedTotalReflectance.jl` is a fast, flexible, and convenient Julia package for simulating **Attenuated Total Reflection (ATR)** using the **Transfer Matrix Method (TMM)** in multilayer thin films. Whether you are modeling surface plasmon resonance (SPR) in a Kretschmann configuration or designing complex optical coatings, this package provides a high-performance backend to calculate reflectance and transmission.

---

## Features

* **High Performance:** Uses Julia's native speed and type-stability for rapid calculations over large parameter spaces.
* **Arbitrary Multilayers:** Simulate a large number of thin-film layers with complex refractive indices ($n + ik$) to account for absorption.
* **Full Polarisation Support:** Supports both transverse magnetic (**p-polarisation**) and transverse electric (**s-polarisation**) modes.
* **Flexible Sweeps:** Easily simulate across ranges of incident angles, wavelengths, or layer thicknesses.

---
