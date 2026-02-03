module AttenuatedTotalReflectance

export 
    layer,
    material_stack,
    complex_n,
    snells_law,
    epsilon_to_nk,
    nk_to_epsilon,
    target_layer,
    refl_coeff_S,
    refl_coeff_P,
    trans_coeff_S,
    trans_coeff_P,
    compute_transfer_coefficients,
    angular_ATR,
    wavelength_ATR

include("Fresnel_Equations.jl")
include("utils.jl")

end
