# Theoretical Introduction: The Transfer Matrix Method (TMM)

To model Attenuated Total Reflection (ATR) and light propagation in multilayered thin films, **AttenuatedTotalReflectance.jl** implements the Transfer Matrix Method (TMM) based on the highly optimized matrix formalism described by Koji Ohta and Hatsuo Ishida (1990). 

This page provides a physical introduction to the mathematical framework utilized under the hood.

---

## 1. Multilayer System Configuration

We consider a stratified, isotropic, and parallel multilayer system consisting of $m$ thin-film layers. This film stack is bounded by an semi-infinite initial incident phase (designated as layer $0$, typically a high-index prism in ATR setups) and a final semi-infinite substrate phase (designated as layer $m+1$, such as air or a cladding medium).
Incident Light (E₀⁺)
            \      ^
             \    /   Reflected Light (E₀⁻)
  Phase 0     \  /

=================/================= Boundary 1
Layer 1

----------------------------------- Boundary 2
Layer 2

:                   :
----------------------------------- Boundary j
Layer j

=================================== Boundary m+1
Phase m+1           ---> Transmitted Light (Eₘ₊₁⁺)

Each layer $j$ is physically characterized by:
* Its physical thickness $h_j$.
* Its complex refractive index $\hat{n}_j$, defined as:
  $$\hat{n}_j = n_j + ik_j$$
  where $n_j$ is the real refractive index and $k_j$ is the extinction coefficient representing material absorption.

When light is incident on the first boundary at an angle $\theta_0$, the complex refractive angles $\theta_j$ in successive layers are strictly governed by Snell's Law:
$$\hat{n}_{j-1} \sin\theta_{j-1} = \hat{n}_j \sin\theta_j \quad (j = 1, 2, \dots, m+1)$$

---

## 2. Abeles's Matrix Formalism

The TMM treats the electromagnetic fields within each layer as a linear superposition of a forward-propagating wave ($E^+$) traveling in the $+z$ direction and a backward-propagating wave ($E^-$) traveling in the $-z$ direction. 

Abeles demonstrated that the field amplitudes at the front interface ($E_0^+, E_0^-$) and the final substrate interface ($E_{m+1}^+, E_{m+1}^-$) can be mapped using a chain of interface propagation matrices:

$$\begin{pmatrix} E_0^+ \\ E_0^- \end{pmatrix} = \frac{C_1 C_2 \dots C_{m+1}}{t_1 t_2 \dots t_{m+1}} \begin{pmatrix} E_{m+1}^+ \\ E_{m+1}^- \end{pmatrix}$$

Where:
* $t_j$ is the local Fresnel transmission coefficient at the interface between layer $j-1$ and layer $j$.
* $C_j$ represents the characteristic propagation matrix for the layer, structured as:

$$C_j = \begin{pmatrix} \exp(-i\delta_{j-1}) & r_j \exp(-i\delta_{j-1}) \\ r_j \exp(i\delta_{j-1}) & \exp(i\delta_{j-1}) \end{pmatrix}$$

Here, $r_j$ is the local Fresnel reflection coefficient. The parameter $\delta_{j-1}$ describes the layer’s optical phase thickness (the phase change experienced as light traverses the layer), formulated as:
$$\delta_0 = 0$$
$$\delta_{j-1} = 2\pi\nu\hat{n}_{j-1}\cos\theta_{j-1}h_{j-1}$$
where $\nu$ represents the wavenumber of the incident light in a vacuum.

### Fresnel Coefficients
The values of $r_j$ and $t_j$ change depending on the polarization of the light wave relative to the plane of incidence:

* **Transverse Magnetic / Parallel ($p$-polarization):**
  $$r_{jp} = \frac{\hat{n}_{j-1}\cos\theta_j - \hat{n}_j\cos\theta_{j-1}}{\hat{n}_{j-1}\cos\theta_j + \hat{n}_j\cos\theta_{j-1}}$$
  $$t_{jp} = \frac{2\hat{n}_{j-1}\cos\theta_{j-1}}{\hat{n}_{j-1}\cos\theta_j + \hat{n}_j\cos\theta_{j-1}}$$

* **Transverse Electric / Perpendicular ($s$-polarization):**
  $$r_{js} = \frac{\hat{n}_{j-1}\cos\theta_{j-1} - \hat{n}_j\cos\theta_j}{\hat{n}_{j-1}\cos\theta_{j-1} + \hat{n}_j\cos\theta_j}$$
  $$t_{js} = \frac{2\hat{n}_{j-1}\cos\theta_{j-1}}{\hat{n}_{j-1}\cos\theta_{j-1} + \hat{n}_j\cos\theta_j}$$

---

## 3. Calculating Total Reflectance ($R$)

Because we assume that the final substrate layer ($m+1$) is semi-infinite, no backward-propagating wave exists deep within it ($E_{m+1}^- = 0$). By defining the total system matrix product as:

$$C_1 C_2 \dots C_{m+1} = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$$

We can easily extract the global amplitude reflection coefficient $r$ and transmission coefficient $t$ for the entire thin-film assembly:
$$r = \frac{E_0^-}{E_0^+} = \frac{c}{a}$$
$$t = \frac{E_{m+1}^+}{E_0^+} = \frac{t_1 t_2 \dots t_{m+1}}{a}$$

The net measurable power **reflectance ($R$)** is the squared magnitude of the complex reflection coefficient:
$$R = |r|^2$$

---

## 4. Internal Electric Field Profiles (Ohta & Ishida's Method)

In Attenuated Total Reflection (ATR) spectroscopy, tracking total reflectance is only half the story. The macroscopic absorption behavior is directly governed by the penetration of the evanescent wave and local electric field intensities ($F(z)$) at specific depths $z$ within the film stack.

Historically, calculating depth-dependent fields required executing computationally intensive matrix inversions (such as Hansen's method). Ohta and Ishida introduced an elegant algorithm that avoids matrix inversions entirely by processing successive matrix products backwards from the substrate side.

By defining a partial backproduct matrix $D_j$ as:
$$D_j = C_{j+1} C_{j+2} \dots C_{m+1} = \begin{pmatrix} a_j & b_j \\ c_j & d_j \end{pmatrix}$$

The forward and backward propagating field amplitudes immediately below any internal interface $j$ are directly related to the input incident field $E_0^+$ via:
$$E_j^+ = t_1 t_2 \dots t_j \frac{a_j}{a} E_0^+$$
$$E_j^- = t_1 t_2 \dots t_j \frac{c_j}{a} E_0^+$$

Once these boundary amplitudes are computed, the field profile at any arbitrary depth $z$ residing within layer $j$ can be continuously evaluated:
$$E^+(z) = E_j^+ \exp(i K_{zj}\Delta z)$$
$$E^-(z) = E_j^- \exp(-i K_{zj}\Delta z)$$

where $K_{zj} = 2\pi\nu\hat{n}_j\cos\theta_j$ is the wavevector component normal to the interfaces, and $\Delta z$ represents the localized distance from the $j$-th boundary boundary. 

This non-inverting architecture allows **AttenuatedTotalReflectance.jl** to solve internal electric fields and point-by-point local absorption characteristics in nearly the exact same processing time required to calculate basic total reflectance.

---

## References
* **Ohta, K., & Ishida, H. (1990).** *Matrix formalism for calculation of electric field intensity of light in stratified multilayered films.* Applied Optics, 29(13), 1952-1959. https://doi.org/10.1364/AO.29.001952