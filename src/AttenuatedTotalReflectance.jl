module AttenuatedTotalReflectance

export 
    layer,
    material_stack,
    complex_n,
    snells_law,
    epsilon_to_nk,
    nk_to_epsilon,
    compute_transfer_coefficents,
    angular_ATR,
    wavelength_ATR

include("Fresnel_Equations.jl")
include("utils.jl")

end
