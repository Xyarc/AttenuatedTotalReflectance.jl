"""
Struct object used to hold information on each layer in the overall system to be simulated

Args:

    - `material::String`: Name of the current layer

    - `n::Float64`: Refractive index of the current layer

    - `k::Float64`: Extinction coefficient of the current layer

    - `thickness::Float64`: Thickness of the current layer in meters [m]

Example:

    ```jldoctest
    
        julia> thin_film = layer("Silver", 0.051585, 3.9046, 20e-9)

        julia> thin_film.n

        julia> thin_film.thickness
    ```
"""
mutable struct layer
    material::String
    n::Float64
    k::Float64
    thickness::Float64
end

"""
Helper function to assemble an array of layer objects to hold information on the material layers of the system to be simulated with the addition of perfectly matched layer.

Args:

    - `number_layers::Int64`: Number of layers in the system 

Return:

    - `stack`: Array of layer objects of length number_layers+1

Mathematical Note:

Perfectly matched layer is required to absorb reflections from the last layer of the system. Both the PML must be larger than the thin films simulated to avoid oscillations in the output of the fresnel calculations

Example:

    ```jldoctest
    
        julia> thin_film_system = material_stack(3)

        julia> thin_film_system[2] = layer("Silver", 0.051585, 3.9046, 20e-9)

        julia> thin_film_system[2].material

        julia> thin_film_system[2].thickness
    ```
"""


function material_stack(number_layers::Int64)
    # Sets up frame work of arbitrary length
    # First and last layers are large compared to rest of stack
    # Last layer (number_layers) must be a Dielectric and PML must be added after this layer (number_layers+1)
    # Create an uninitialise array to hold Layer information
    stack = Array{Layer, 1}(undef, number_layers+1)
    
    # Assign layers to the stack
    for i in 1:number_layers
        if i == 1
            stack[i] = Layer("Substrate", 0.0, 0.0, 10e-6 +1e-7) # Make the substrate much larger than thin films
        elseif i == number_layers
            stack[i] = Layer("Last Layer", 0.0, 0.0, 5e-6)
        else
            stack[i] = Layer("Material Name", 0.0, 0.0, 0.0)
        end
        
    end
    stack[end] = Layer("PML", 0.0, 0.0, 1e-6)
    return stack
end


"""
Helper function to take n & k values to produce singular complex refractive index

Args:
    - n: Refractive index
    - k: Extinction coefficient

Return:
    - Complex index of refraction

Example:
        Taking the refractive index of Silver @ 586.6 nm: n = 0.051585 & k = 3.9046

        ```jldoctest
        julia> complex_n(0.051585, 3.9046)
        ```
"""
function complex_n(n::Float64, k::Float64=0.0)
    return n + k * 1.0im
end

"""
Implimentation of Snell's Law to compute the angle of refraction. Computes the resulting angle due to the refraction between two materials of differing refractive indices. Returned angle can be complex.

Args:

    - `theta_in::Number`: Angle of the incident light ray traveling through a material of refractive index, n_in, with repect to the normal prerpendicular to the interface.
    - `n_in::Number`: Refractive index of the first material the light passes through. 
    - `n_out::Number`: Refractive index of the second material the light passes through. 

    All Args may be complex

Return:

    Resulting angle due to refraction at the interface.

Example:

    Taking the refractive index of Silver @ 586.6 nm: e1 = -15.243 & e2 = 0.40284
    
    ```jldoctest
    julia> snells_law(pi/3, 1, 1.33)

    julia> snells_law(pi/4, 1.33, 1.0)
    ```    
"""
function snells_law(theta_in::Number, n_in::Number, n_out::Number)
    arg = (n_in / n_out) * sin(theta_in)
    # Adding "complex" handles ratios larger than 1 
    return asin(complex(arg))
end

"""
Helper function to convert dielectric constant into n and k.
    
Args:

    - `e1::Float64`: Real component of dielectric constant.

    - `e2::Float64`: Imaginary component of dielectric constant.

Return:

    - n: Refractive index

    - k: Extinction coefficient

Example:

    Taking the refractive index of Silver @ 586.6 nm: e1 = -15.243 & e2 = 0.40284
    
    ```jldoctest
    julia> epsilon_to_nk(-15.243, 0.40284)
    ```
"""

function epsilon_to_nk(e1::Float64, e2::Float64)
    n = sqrt(0.5 * (sqrt(e1^2 + e2^2) + e1))
    k = sqrt(0.5 * (sqrt(e1^2 + e2^2) - e1))

    return n, k
end

"""
Helper function to convert n and k into the components of the dielectric constant.

Args:

    - `n::Float64`: Refractive index
    - `k::Float64`: Extinction coefficient

Return:

    - e1: Real component of dielectric constant.
    - e2: Imaginary component of dielectric constant.

Example:

    Taking the refractive index of Silver @ 586.6 nm: n = 0.051585 & k = 3.9046
    
    ```jldoctest
    julia> nk_to_epsilon(0.051585, 3.9046)
    ```
"""

function nk_to_epsilon(n::Float64, k::Float64)
    e1 = n^2 - k^2
    e2 = 2 * n * k

    return e1, e2
end
