"""
Features on/off of the SPAC model
"""
Base.@kwdef mutable struct SPACFeatures{FT<:AbstractFloat}
    "Whether to acclimate leaf Vcmax and Jmax TD"
    ACCLIMATE_T_VCMAX::Bool = true
    "Number of temperature memory (hours for hourly simulation)"
    ACCLIMATE_T_VCMAX_N::Int = 240
    "Allow leaf condensation"
    ALLOW_LEAF_CONDENSATION::Bool = false
    "Allow leaf regrowth when junction pressure is higher than -0.1 MPa"
    ALLOW_LEAF_REGROWTH::Bool = false
    "Allow leaf shedding in prescibe LAI mode to avoid numerical issues"
    ALLOW_LEAF_SHEDDING::Bool = false
    "Allow xylem to grow"
    ALLOW_XYLEM_GROWTH::Bool = false
    "Effective leaf spectra based on CI effect"
    EFFECTIVE_LEAF_SPECTRA::Bool = true
    "Enable the chemical energy related to photosynthesis and respiration"
    ENABLE_CHEMICAL_ENERGY::Bool = true
    "Enable drought legacy effect"
    ENABLE_DROUGHT_LEGACY::Bool = false
    "Whether to compute canopy reflectance"
    ENABLE_REF::Bool = true
    "Whether to prescribe a non-zero stem area index (SAI) per pixel; if false, SAI stays 0 everywhere"
    ENABLE_SAI::Bool = true
    "Whether to compute fluorescence"
    ENABLE_SIF::Bool = true
    "Apply sigmoid cutoff to leaf-level SIF emission PDF based on excitation wavelength"
    ENABLE_LEAF_SIF_SIGMOID::Bool = true
    "Rescale leaf-level SIF emission PDF to integrate to 1 after sigmoid cutoff"
    ENABLE_LEAF_SIF_RESCALE::Bool = true
    "Apply sigmoid cutoff to chloroplast-level SIF emission PDF based on excitation wavelength"
    ENABLE_CHL_SIF_SIGMOID::Bool = true
    "Rescale chloroplast-level SIF emission PDF to integrate to 1 after sigmoid cutoff"
    ENABLE_CHL_SIF_RESCALE::Bool = true
    "Whether to clamp p_sun_sensor (Pso) to min(Pso, Po, Ps); set false to match 4SAIL (allow non-physical overshoot)"
    ENFORCE_PSO_CLAMP::Bool = true
    "Fix the TD of η for Johnson-Berry model"
    FIX_ETA_TD::Bool = true
    "Threshold of the critical pressure or flow that trigger root disconnection"
    KR_ROOT_DISCONNECTION::FT = 0.5
    "Threshold of the critical pressure or flow that trigger a remainder of conductance"
    KR_THRESHOLD::FT = 0.001
    "Prescribe air layer information such as partial pressures"
    PRESCRIBE_AIR::Bool = true
    "Unlimited NSC pool in the plant"
    UNLIMITED_NSC_POOL::Bool = true

    "New C4 ∂A∂E and ∂Θ∂E methods"
    NEW_C4_STOMATAL_METHODS::Bool = false
end;
