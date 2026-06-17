"""

    leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}) where {FT}

Update the SIF conversion matrix of the leaf without reabsorption, given
- `config` SPAC configuration
- `bio` leaf biophysics

"""

function leaf_sif_matrices_chl! end;

leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}) where {FT} = leaf_sif_matrices_chl!(config, bio, cache, config.METHODS.FLUORESCENCE_SPECTRA_METHOD);

leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}, ::PlatespectFluorescenceSpectra) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    (; IΛ_SIF, IΛ_SIFE, ΔΛ_SIF, Λ_SIF, Λ_SIFE, Φ_PS) = SPECTRA;

    # update the SIF emission vector per excitation wavelength
    ϕ      = bio.auxil._ϕ_sif;
    factor = cache.cache_sif_1;

    for i in eachindex(IΛ_SIFE)
        ii = IΛ_SIFE[i];

        # read the SIF emission spectrum
        ϕ .= view(Φ_PS, IΛ_SIF);

        # tune SIF emission PDF based on the SIF excitation wavelength
        expsife = exp(Λ_SIFE[ii] / 10);
        @. factor = 1 / (1 + exp(-Λ_SIF / 10) * expsife);
        ϕ .*= factor;

        # rescale ϕ
        ϕ ./= ΔΛ_SIF' * ϕ;

        # matꜛ_chl = α * f_sife * ϕ / 2 (mat_b_chl + mat_f_chl; b/f partition cancels in total chl SIF)
        view(bio.auxil.matꜛ_chl, :, i) .= bio.auxil.α_leaf[ii] * bio.auxil.f_sife[ii] .* ϕ ./ 2;
    end;

    return nothing
);
