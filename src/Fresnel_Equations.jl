# Load Dependancies
using StaticArrays
using LinearAlgebra


# ---- Internal Calculation Functions -----

function phase_change(wavelength::Float64, ni::ComplexF64, thickness::Float64, theta_i::ComplexF64)
    # Computes the phase change through a layer of given thickness between the upper and lower boundry of the material 
    v = 1/wavelength
    
    cos_term = cos(theta_i)

    phase_change = 2pi*v*ni*thickness*cos_term

    return phase_change
end

# Fresnel Reflection Coefficients for P & S polarisation

function refl_coeff_S(theta_i, ni::ComplexF64, theta_j, nj::ComplexF64)
    #Computes the Fresnel Reflection coefficient between the i th (j-1) and the j th layers for S-polarisation
    return (ni*cos(theta_i) - nj*cos(theta_j))/(ni*cos(theta_i) + nj*cos(theta_j))
end

function refl_coeff_P(theta_i, ni::ComplexF64, theta_j, nj::ComplexF64)
    #Computes the Fresnel Reflection coefficient between the i th (j-1) and the j th layers for P-polarisation
    return (ni*cos(theta_j) - nj*cos(theta_i))/(ni*cos(theta_j) + nj*cos(theta_i))
end

# Fresnel Transmition Coefficients for P & S polarisation

function trans_coeff_S(theta_i, n_i::ComplexF64, theta_j, n_j::ComplexF64)
    #Computes the Fresnel Transmission coefficient between the i th (j-1) and the j th layers for S-polarisation
    return (2*n_i*cos(theta_i))/(n_i*cos(theta_i) + n_j*cos(theta_j))
end

function trans_coeff_P(theta_i, n_i::ComplexF64, theta_j, n_j::ComplexF64)
    #Computes the Fresnel Transmission coefficient between the i th (j-1) and the j th layers for S-polarisation
    return (2*n_i*cos(theta_i))/(n_i*cos(theta_j) + n_j*cos(theta_i))
end

function wave_vector_Z(wavelength::Float64, n::ComplexF64, theta)
    #Computes the z-component of the wavevector for a layer j
    return ((2*pi)/wavelength)*(n*cos(theta))
end

function wave_vector_X(wavelength::Float64, n::ComplexF64, theta)
    #Computes the x-component of the wavevector for a layer j
    return ((2*pi)/wavelength)*(n*sin(theta))
end


#-------- Core Fresnel Computation --------------
"""
compute_transfer_coefficients(stack, theta, wavelength, layer_j=0; S=false)

Calculate the optical transfer matrices and transmission coefficients for a multilayer 
thin-film stack using the Transfer Matrix Method (TMM).

# Arguments
- `stack::Vector`: A vector of layer objects. Each element must have fields `.n` (refractive index), 
  `.k` (extinction coefficient), and `.thickness`. See `material_stack` function in utils.jl
- `theta::Number`: Incident angle in radians.
- `wavelength::Number`: Wavelength of the incident light (units must match `stack[i].thickness`).
- `layer_j::Int`: (Optional) The index of a specific layer of intrest. If provided, the function 
  calculates partial matrices and coefficients split at this layer. 
- `S::Bool`: (Keyword) If `true`, computes coefficients for S-polarization (TE). 
  If `false` (default), computes for P-polarization (TM).

# Returns
A `Tuple` containing:
1. `C::SMatrix{2,2,ComplexF64}`: The total system transfer matrix.
2. `D0::SMatrix{2,2,ComplexF64}`: The partial transfer matrix from the ambient to `layer_j`.
3. `Dj::SMatrix{2,2,ComplexF64}`: The partial transfer matrix from `layer_j` to the substrate.
4. `t_total::ComplexF64`: The total Fresnel transmission coefficient for the entire stack.
5. `t0::ComplexF64`: Partial transmission coefficient from ambient to `layer_j`.
6. `tj::ComplexF64`: Partial transmission coefficient from `layer_j` to substrate.

# Mathematical Note
The transfer matrix \$ M\$ for an interface and subsequent layer propagation is defined such that:
\$\$\begin{pmatrix} E_{i-1}^+ \\ E_{i-1}^- \end{pmatrix} = M_i \begin{pmatrix} E_{i}^+ \\ E_{i}^- \end{pmatrix}\$\$
where \$ E^+\$ and \$ E^-\$ represent the forward and backward traveling electric field components.

This algorithm is based on the work by Koji Ohta and Hatsuo Ishida (DOI: 10.1364/ao.29.001952)

This method also assumes that the system is not magnetic
"""
function compute_transfer_coefficents(stack::Vector, theta::Number, wavelength::Number, layer_j::Int64 = 0; S::Bool = false)
    
    # Determine if P or S polarisation function should be used
    refl_func  = S ? Refl_Coeff_S : Refl_Coeff_P
    trans_func = S ? Trans_Coeff_S : Trans_Coeff_P
    
    # 2. Initialise Variables with Static Types for Efficiency
    # SMatrix{Rows, Cols, Type, TotalElements}
    identity_2x2 = one(SMatrix{2, 2, ComplexF64, 4})
    
    C  = identity_2x2
    D0 = identity_2x2
    Dj = identity_2x2

    # Define t_total as the holding varibales for the transmission coefficients
    t_total = 1.0 + 0.0im
    t0 = 1.0 + 0.0im # Partial Trans_Coeff below jth layer
    tj = 1.0 + 0.0im # Partial Trans_Coeff above jth layer
   
    # Holding Variables for Complex Refractive index of the zeroth level, layer below {i} and layer above {j}
    n0 = complex_n(stack[1].n, stack[1].k)
    
    
    
    #3. Iterative Computation over all layers in the stack
    for i in 1:length(stack)
        nj = complex_n(stack[i].n, stack[i].k)
        theta_j = snells_law(theta, n0, nj)

        # Calculate interface and propagation values
        if i == 1
            r_ij = refl_func(theta, n0, theta_j, nj)
            t_ij = trans_func(theta, n0, theta_j, nj)
            
            # Create a static matrix for the first interface
            M = SMatrix{2, 2, ComplexF64, 4}(
                1.0,  r_ij, 
                r_ij, 1.0
            )
        else
            ni = complex_n(stack[i-1].n, stack[i-1].k)
            theta_i = snells_law(theta, n0, ni)
            
            phi = phase_change(wavelength, ni, stack[i-1].thickness, theta_i)
            r_ij = refl_func(theta_i, ni, theta_j, nj)
            t_ij = trans_func(theta_i, ni, theta_j, nj)

            # Create a static matrix for propagation + interface
            # Note: SMatrix construction is column-major by default
            c_neg = cis(-phi)
            c_pos = cis(phi)
            M = SMatrix{2, 2, ComplexF64, 4}(
                c_neg,        r_ij * c_pos, 
                r_ij * c_neg, c_pos
            )
        end

        # Update total transfer matrix
        C = M * C
        t_total *= t_ij

        # 4. Handle Partial Coefficients
        if i == layer_j
            # Assign partial transefer matrix and transmission coeffiecents below target layer
            D0 = C
            t0 = t_total

        elseif i > layer_j && layer_j != 0
            # Assign partial transfer matrix and transmission coeffiecents above target layer
            Dj = M * Dj
            tj *= t_ij
        end
    end

    return C, D0, Dj, t_total, t0, tj
end

"""
    angular_ATR(stack, theta_range, wavelength, d=0.0, metal_layer=0; S=false)

Perform an angle dependant Attenuated Total Reflection (ATR) simulation, calculating 
reflectivity, transmittance, and the electric field enhancement at a specific location for each angle.

# Arguments
- `stack::Vector`: A vector of layer objects. Each must contain `.n`, `.k`, and `.thickness`. See `material_stack`
- `theta_range::AbstractVector`: A range or vector of incident angles in radians.
- `wavelength::Number`: The vacuum wavelength of the incident light.
- `d::Number`: The distance from the interface of the `target_layer` (in the same units 
  as `wavelength`) at which to compute the electric field.
- `metal_layer::Int`: An index used by `Target_Layer` to identify the specific layer 
  of interest (e.g., the gold layer in a Kretschmann configuration).
- `S::Bool`: If `true`, calculates for S-polarization (TE). If `false` (default), 
  calculates for P-polarization (TM).

# Returns
A `Tuple` containing eight vectors (each of length `length(theta_range)`):
1.  `reflectivity`: Intensity reflection coefficient (R).
2.  `transmittance`: Intensity transmission coefficient (T) into the final medium.
3.  `f_x`: Normalized field intensity component \$|E_x/E_0|^2\$.
4.  `f_y`: Normalized field intensity component \$|E_y/E_0|^2\$.
5.  `f_z`: Normalized field intensity component \$|E_z/E_0|^2\$.
6.  `f_p`: Total normalized field intensity enhancement \$|E_{total}/E_0|^2\$.


# Physics Note
Field enhancement is calculated as the ratio of the local field intensity to the 
incident field intensity.
- For **P-polarization**, the field exists in the plane of incidence: 
  \$E_{total} = \\sqrt{|E_x|^2 + |E_z|^2}\$.
- For **S-polarization**, the field is perpendicular to the plane of incidence: 
  \$E_{total} = |E_y|\$.

"""

function angular_ATR(stack::Vector, theta_range::AbstractVector, wavelength::Number, 
                             d::Number = 0.0, metal_layer::Int = 0; S::Bool = false)
    
    # 1. Pre-allocate output arrays
    num_pts      = length(theta_range)
    reflectivity = zeros(Float64, num_pts)
    transmittance = zeros(Float64, num_pts)
    f_x          = zeros(Float64, num_pts)
    f_y          = zeros(Float64, num_pts)
    f_z          = zeros(Float64, num_pts)
    f_p          = zeros(Float64, num_pts)

    n0    = Complex_n(stack[1].n, stack[1].k)
    n_end = Complex_n(stack[end].n, stack[end].k)
    layer_target, n_target = target_layer(stack, d, metal_layer)

    e0_f = 1.0 + 0.0im 

    for (i, angle) in enumerate(theta_range)
        angle_end    = snells_law(angle, n0, n_end)
        theta_target = snells_law(angle, n0, n_target)
        
        a = real((conj(n_end) * cos(angle_end)) / (conj(n0) * cos(angle)))

        C, D0, Dj, t_total, t0, tj = compute_transfer_coefficients(
            stack, angle, wavelength, layer_target; S=S
        )

        inv_c11 = 1.0 / C[1,1]
        reflectivity[i]  = abs2(C[2,1] * inv_c11)
        transmittance[i] = a * abs2(t_total * inv_c11)
          
        # Field Enhancement Logic
        kj_z = wave_vector_Z(wavelength, n_target, theta_target)
        ej_f = t0 * (Dj[1,1] * inv_c11) * e0_f
        ej_b = t0 * (Dj[2,1] * inv_c11) * e0_f

        e_f = ej_f * cis(kj_z * d)
        e_b = ej_b * cis(-kj_z * d)

        if S
            # S-polarization (TE): E-field is only along the y-axis
            # E_y = E_f + E_b
            f_y[i] = abs2(e_f + e_b)
            f_x[i] = 0.0
            f_z[i] = 0.0
            f_p[i] = f_y[i] # Total power is just Ey^2
        else
            # P-polarization (TM): E-field has x and z components
            f_x[i] = abs2((e_f - e_b) * cos(theta_target))
            f_z[i] = abs2((e_f + e_b) * sin(theta_target))
            f_y[i] = 0.0
            f_p[i] = f_x[i] + f_z[i] # Total power is Ex^2 + Ez^2
        end
        
        
    end

    return reflectivity, transmittance, f_x, f_y, f_z, f_p
end

"""
    wavelength_ATR(stack, theta, wavelength_range, d=0.0, metal_layer=0; S=false)

Perform a wavelength dependant Attenuated Total Reflection (ATR) simulation, calculating 
reflectivity, transmittance, and the electric field enhancement at a specific location for a fixed angle.

# Arguments
- `stack::Vector`: A vector of layer objects. Each must contain `.n`, `.k`, and `.thickness`. See `material_stack`
- `theta::Number`: A range or vector of incident angles in radians.
- `wavelength_range::AbstractVector`: The vacuum wavelength of the incident light.
- `d::Number`: The distance from the interface of the `target_layer` (in the same units 
  as `wavelength`) at which to compute the electric field.
- `metal_layer::Int`: An index used by `Target_Layer` to identify the specific layer 
  of interest (e.g., the gold layer in a Kretschmann configuration). Used for plasmonic systems
- `S::Bool`: If `true`, calculates for S-polarization (TE). If `false` (default), 
  calculates for P-polarization (TM).

# Returns
A `Tuple` containing eight vectors (each of length `length(theta_range)`):
1.  `reflectivity`: Intensity reflection coefficient (R).
2.  `transmittance`: Intensity transmission coefficient (T) into the final medium.
3.  `f_x`: Normalized field intensity component \$|E_x/E_0|^2\$.
4.  `f_y`: Normalized field intensity component \$|E_y/E_0|^2\$.
5.  `f_z`: Normalized field intensity component \$|E_z/E_0|^2\$.
6.  `f_p`: Total normalized field intensity enhancement \$|E_{total}/E_0|^2\$.


# Physics Note
Field enhancement is calculated as the ratio of the local field intensity to the 
incident field intensity.
- For **P-polarization**, the field exists in the plane of incidence: 
  \$E_{total} = \\sqrt{|E_x|^2 + |E_z|^2}\$.
- For **S-polarization**, the field is perpendicular to the plane of incidence: 
  \$E_{total} = |E_y|\$.

"""

function wavelength_ATR(stack::Vector, theta::Number, wavelength_range::AbstractVector, 
                             d::Number = 0.0, metal_layer::Int = 0; S::Bool = false)
    
    # 1. Pre-allocate output arrays
    num_pts      = length(wavelength_range)
    reflectivity = zeros(Float64, num_pts)
    transmittance = zeros(Float64, num_pts)
    f_x          = zeros(Float64, num_pts)
    f_y          = zeros(Float64, num_pts)
    f_z          = zeros(Float64, num_pts)
    f_p          = zeros(Float64, num_pts)

    n0    = Complex_n(stack[1].n, stack[1].k)
    n_end = Complex_n(stack[end].n, stack[end].k)
    layer_target, n_target = target_layer(stack, d, metal_layer)

    e0_f = 1.0 + 0.0im 

    for (i, wavelength) in enumerate(wavelength_range)
        angle_end    = snells_law(angle, n0, n_end)
        theta_target = snells_law(angle, n0, n_target)
        
        a = real((conj(n_end) * cos(angle_end)) / (conj(n0) * cos(angle)))

        C, D0, Dj, t_total, t0, tj = compute_transfer_coefficients(
            stack, theta, wavelength, layer_target; S=S
        )

        inv_c11 = 1.0 / C[1,1]
        reflectivity[i]  = abs2(C[2,1] * inv_c11)
        transmittance[i] = a * abs2(t_total * inv_c11)
          
        # Field Enhancement Logic
        kj_z = wave_vector_Z(wavelength, n_target, theta_target)
        ej_f = t0 * (Dj[1,1] * inv_c11) * e0_f
        ej_b = t0 * (Dj[2,1] * inv_c11) * e0_f

        e_f = ej_f * cis(kj_z * d)
        e_b = ej_b * cis(-kj_z * d)

        if S
            # S-polarization (TE): E-field is only along the y-axis
            # E_y = E_f + E_b
            f_y[i] = abs2(e_f + e_b)
            f_x[i] = 0.0
            f_z[i] = 0.0
            f_p[i] = f_y[i] # Total power is just Ey^2
        else
            # P-polarization (TM): E-field has x and z components
            f_x[i] = abs2((e_f - e_b) * cos(theta_target))
            f_z[i] = abs2((e_f + e_b) * sin(theta_target))
            f_y[i] = 0.0
            f_p[i] = f_x[i] + f_z[i] # Total power is Ex^2 + Ez^2
        end
        
        
    end

    return reflectivity, transmittance, f_x, f_y, f_z, f_p
end