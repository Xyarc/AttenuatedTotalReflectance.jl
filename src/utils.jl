"""
Struct object used to hold information on each layer in the overall system to be simulated

Arguments:

    - `material::String`: Name of the current layer

    - `n::Float64`: Refractive index of the current layer

    - `k::Float64`: Extinction coefficient of the current layer

    - `thickness::Float64`: Thickness of the current layer in meters [m]. Must be positive.

Example:

    ```jldoctest
    
        julia> thin_film = layer("Silver", 0.051585, 3.9046, 20e-9)

        julia> thin_film.n
        0.051585
        julia> thin_film.thickness
        20e-9
    ```
"""
mutable struct layer
    material::String
    n::Float64
    k::Float64
    thickness::Float64

    function layer(material, n, k, thickness)
        thickness < 0 && throw(ArgumentError("Thickness cannot be negative"))
        new(material, n, k, thickness)
    end
end

"""
Helper function to assemble an array of layer objects to hold information on the material layers of the system to be simulated with the addition of perfectly matched layer.

Arguments:

    - `materials::Vector{<:NamedTuple}`: A vector containing information on the materials to use. (Name, n, k, thickness [m])

Return:

    - `stack`: Array of layer objects of length number_layers+1

Mathematical Note:

Perfectly matched layer is required to absorb reflections from the last layer of the system. Both the PML must be larger than the thin films simulated to avoid oscillations in the output of the fresnel calculations

Example:

    ```jldoctest
    
        julia> mats = [
            (material="Glass", n=1.5, k=0.0, thickness=5e-6),
            (material="Gold", n=0.18, k=3.0, thickness=50e-9),
            (material="Air", n=1.0, k=0.0, thickness=5e-6),
        ]

        julia> stack = material_stack(mats)

        julia> stack[1].material
        "Substrate: Glass"
        julia> stack[2].n
        0.18
        julia> stack[2].thickness
        50e-9
    ```
"""
function material_stack(materials::Vector{<:NamedTuple})
    # Sets up frame work of arbitrary length based on an vector of input materials
    # First and last layers normally should large compared to rest of stack
    # Last layer must be a Dielectric and PML must be added after this layer (number_layers+1)
    # Create an uninitialise array to hold Layer information
    
    stack = Array{layer, 1}(undef, length(materials)+1)
    
    # Assign layers to the stack
    for i in 1:length(materials)
        m = materials[i]
        if i == 1
            # Throw warning if substrate layer is thin
            if m.thickness < 1e-6
                 # Make the substrate much larger than thin films
                print("Warning: Substrate layer is thin!")
            end
        
            stack[i] = layer("Substrate: $(m.material)", m.n, m.k, m.thickness)
       
        else
            stack[i] = layer(m.material, m.n, m.k, m.thickness)
        end
        
    end

    if stack[end-1].thickness < 1e-6
        print("Warning: Final dielectric layer is thin!")
    end
    # Add PML to the end of the simulation space

    stack[end] = layer("PML", stack[end - 1].n + 1e-5, 0.0, 1e-6)
    return stack
end


"""
Helper function to take n & k values to produce singular complex refractive index

Arguments:
    - n: Refractive index
    - k: Extinction coefficient

Return:
    - Complex index of refraction

Example:
        Taking the refractive index of Silver @ 586.6 nm: n = 0.051585 & k = 3.9046

        ```jldoctest
        julia> complex_n(0.051585, 3.9046)
        0.051585 + 3.9046im
        ```
"""
function complex_n(n::Float64, k::Float64=0.0)
    return n + k * 1.0im
end

"""
Implimentation of Snell's Law to compute the angle of refraction. Computes the resulting angle due to the refraction between two materials of differing refractive indices. Returned angle can be complex.

Arguments:

    - `theta_in::Number`: Angle of the incident light ray traveling through a material of refractive index, n_in, with repect to the normal prerpendicular to the interface.
    - `n_in::Number`: Refractive index of the first material the light passes through. 
    - `n_out::Number`: Refractive index of the second material the light passes through. 

    All Arguments may be complex

Return:

    Resulting angle due to refraction at the interface.

Example:

    Taking the refractive index of Silver @ 586.6 nm: e1 = -15.243 & e2 = 0.40284
    
    ```jldoctest
    julia> snells_law(pi/3, 1, 1.33)
    0.7091 + 0.0im
    julia> snells_law(pi/4, 1.33, 1.0)
    1.2239 + 0.0im
    ```    
"""
function snells_law(theta_in::Number, n_in::Number, n_out::Number)
    arg = (n_in / n_out) * sin(theta_in)
    # Adding "complex" handles ratios larger than 1 
    return asin(complex(arg))
end

"""
Helper function to convert dielectric constant into n and k.
    
Arguments:

    - `e1::Float64`: Real component of dielectric constant.

    - `e2::Float64`: Imaginary component of dielectric constant.

Return:

    - n: Refractive index

    - k: Extinction coefficient

Example:

    Taking the refractive index of Silver @ 586.6 nm: e1 = -15.243 & e2 = 0.40284
    
    ```jldoctest
    julia> epsilon_to_nk(-15.243, 0.40284)
    (0.0516, 3.9045)
    ```
"""
function epsilon_to_nk(e1::Float64, e2::Float64)
    n = sqrt(0.5 * (sqrt(e1^2 + e2^2) + e1))
    k = sqrt(0.5 * (sqrt(e1^2 + e2^2) - e1))

    return n, k
end

"""
Helper function to convert n and k into the components of the dielectric constant.

Arguments:

    - `n::Float64`: Refractive index
    - `k::Float64`: Extinction coefficient

Return:

    - e1: Real component of dielectric constant.
    - e2: Imaginary component of dielectric constant.

Example:

    Taking the refractive index of Silver @ 586.6 nm: n = 0.051585 & k = 3.9046
    
    ```jldoctest
    julia> nk_to_epsilon(0.051585, 3.9046)
    (-15.243, 0.402)
    ```
"""
function nk_to_epsilon(n::Float64, k::Float64)
    e1 = n^2 - k^2
    e2 = 2 * n * k

    return e1, e2
end

"""
    target_layer(stack, d, metal_layer_index)

Helper function which evaluates the correct layer and its complex refractive index for field enhancement 
calculations.

# Arguments
- `stack::Vector{Layer}`: The multilayer stack.
- `depth::Number`: Distance above the layer of intrest (0.0 represents the interface).
- `layer_of_intrest::Int`: The index of the layer of interest (e.g., a gold thin film). 
  

# Returns
- `(layer_idx, n_complex)`: A tuple containing the resolved layer index and its 
  complex refractive index.
"""
function target_layer(stack::Vector, depth::Float64, layer_of_intrest::Int64)
    #Computes the layer in which depth belongs and return that layer + refractive index of that layer
    layer_of_intrest < 1 && throw(ArgumentError("Must provide a layer of interest to compute electric field enhancement!"))
    i = 1
    target_layer = layer_of_intrest
    stack_h = stack[target_layer+i].thickness

    while depth > stack_h
        i += 1
        stack_h += stack[target_layer+i].thickness
    end
    # Sets the target layer from the metal interface to the layer which depth z belongs and calculates refractive index of this layer
    target_layer += i
    n_target = complex_n(stack[target_layer].n, stack[target_layer].k) # Refractive index of the layer depth d belongs
    return target_layer, n_target
end