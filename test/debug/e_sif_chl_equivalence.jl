# Equivalence checks for simplified e_sif_chl Step 0 and direct matꜛ_chl leaf matrix.
# Run: julia --project=/path/to/SIFEsc/EmeraldTestV2 test/debug/e_sif_chl_equivalence.jl

using Test
using LinearAlgebra: mul!
using Statistics: mean
using PkgUtility.UniversalConstants: energy_to_photon!, photon_to_energy!

import Emerald.CanopyOptics as CO
import Emerald.SPAC as ESPAC
import Emerald.Namespace as ENS

const FT = Float64;
const RTOL = 1e-10;
const ATOL = 1e-10;

"""Expanded leaf chl matrices (mat_b_chl, mat_f_chl → M±) for reference."""
function mat_chl_expanded_reference!(config::ENS.SPACConfig{FT}, bio::ENS.LeafBio{FT}, cache::ENS.SPACCache{FT}) where {FT}
    (; SPECTRA) = config.CONSTANTS;
    (; IΛ_SIF, IΛ_SIFE, ΔΛ_SIF, Λ_SIF, Λ_SIFE, Φ_PS) = SPECTRA;

    ϕ = bio.auxil._ϕ_sif;
    factor = cache.cache_sif_1;
    mat_b_ref = similar(bio.auxil.matꜛ_chl);
    mat_f_ref = similar(bio.auxil.matꜛ_chl);
    matꜛ_ref = similar(bio.auxil.matꜛ_chl);
    matꜜ_ref = similar(bio.auxil.matꜛ_chl);

    for i in eachindex(IΛ_SIFE)
        ii = IΛ_SIFE[i];
        ϕ .= view(Φ_PS, IΛ_SIF);
        expsife = exp(Λ_SIFE[ii] / 10);
        @. factor = 1 / (1 + exp(-Λ_SIF / 10) * expsife);
        ϕ .*= factor;
        ϕ ./= ΔΛ_SIF' * ϕ;

        vec_b = view(bio.auxil.mat_b, :, i);
        vec_f = view(bio.auxil.mat_f, :, i);
        denom = vec_b .+ vec_f;
        mat_b_col = view(mat_b_ref, :, i);
        mat_f_col = view(mat_f_ref, :, i);
        direct = bio.auxil.α_leaf[ii] .* bio.auxil.f_sife[ii] .* ϕ ./ 2;
        for j in eachindex(denom)
            if denom[j] == 0
                mat_b_col[j] = direct[j];
                mat_f_col[j] = direct[j];
            else
                mat_b_col[j] = bio.auxil.α_leaf[ii] * bio.auxil.f_sife[ii] * ϕ[j] * vec_b[j] / denom[j];
                mat_f_col[j] = bio.auxil.α_leaf[ii] * bio.auxil.f_sife[ii] * ϕ[j] * vec_f[j] / denom[j];
            end;
        end;
    end;

    matꜛ_ref .= (mat_b_ref .+ mat_f_ref) ./ 2;
    matꜜ_ref .= (mat_b_ref .- mat_f_ref) ./ 2;

    return matꜛ_ref, matꜜ_ref
end;

function populate_phi_f!(config::ENS.SPACConfig{FT}, spac::ENS.BulkSPAC{FT}) where {FT}
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);
    (; DIM_AZI, DIM_INCL, DIM_PPAR_BINS) = config.DIMENSIONS;

    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        sen_geo.auxil.ϕ_f_shaded[irt] = leaf.photosystem.auxil.ϕ_f[end];
        if isnothing(DIM_PPAR_BINS)
            for i in 1:DIM_AZI
                sen_geo.auxil.ϕ_f_sunlit[irt][:, i] .= view(leaf.photosystem.auxil.ϕ_f, (i - 1) * DIM_INCL + 1:i * DIM_INCL);
            end;
        else
            for i in 1:DIM_INCL, j in 1:DIM_AZI
                sen_geo.auxil.ϕ_f_sunlit[irt][i, j] = leaf.photosystem.auxil.ϕ_f[sun_geo.auxil.ppar_index[i, j, irt]];
            end;
        end;
    end;

    return nothing
end;

"""Old Step-0 canopy assembly (up+down LIDF terms) → total e_sif_chl."""
function e_sif_chl_expanded_reference!(config::ENS.SPACConfig{FT}, spac::ENS.BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    sen_geo = spac.canopy.sensor_geometry;
    sun_geo = spac.canopy.sun_geometry;
    n_layer = length(leaves);
    (; SPECTRA) = config.CONSTANTS;

    a_leaf = spac.cache.cache_sife_1;
    a_stem = spac.cache.cache_sife_2;
    f_leaf = spac.cache.cache_sife_3;
    e_ref = zeros(FT, size(sun_geo.auxil.e_sif_chl)...);

    @inline local_lidf_weight(mat_0, mat_1) = (
        sun_geo.auxil._mat_incl_azi .= mat_0 .* mat_1;
        mul!(sun_geo.auxil._vec_azi, sun_geo.auxil._mat_incl_azi', can_str.auxil.p_incl_leaf);
        mean(sun_geo.auxil._vec_azi)
    );
    @inline local_lidf_weight(mat_1) = (
        mul!(sun_geo.auxil._vec_azi, mat_1', can_str.auxil.p_incl_leaf);
        mean(sun_geo.auxil._vec_azi)
    );
    _COS²_Θ_INCL_AZI = spac.cache.cache_incl_azi_1;
    _COS²_Θ_INCL_AZI .= (cosd.(config.DIMENSIONS.Θ_INCL) .^ 2);

    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaf = leaves[ilf];
        matꜛ_ref, matꜜ_ref = mat_chl_expanded_reference!(config, leaf.bio, spac.cache);

        a_leaf .= view(leaf.bio.auxil.α_leaf, SPECTRA.IΛ_SIFE) .* can_str.trait.δlai[irt];
        a_stem .= (1 .- view(SPECTRA.ρ_STEM, SPECTRA.IΛ_SIFE)) .* can_str.trait.δsai[irt];
        f_leaf .= a_leaf ./ (a_leaf .+ a_stem);

        sun_geo.auxil._e_dirꜜ_sife .= view(sun_geo.auxil.e_dirꜜ, SPECTRA.IΛ_SIFE, irt) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜜ_sife .= view(sun_geo.auxil.e_difꜜ, SPECTRA.IΛ_SIFE, irt) .* f_leaf .* SPECTRA.ΔΛ_SIFE;
        sun_geo.auxil._e_difꜛ_sife .= view(sun_geo.auxil.e_difꜛ, SPECTRA.IΛ_SIFE, irt + 1) .* f_leaf .* SPECTRA.ΔΛ_SIFE;

        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_dirꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜜ_sife);
        energy_to_photon!(SPECTRA.Λ_SIFE, sun_geo.auxil._e_difꜛ_sife);

        mul!(sun_geo.auxil._e_dirꜜ_sifꜛ, matꜛ_ref, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_dirꜜ_sifꜜ, matꜜ_ref, sun_geo.auxil._e_dirꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜛ, matꜛ_ref, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜜ_sifꜜ, matꜜ_ref, sun_geo.auxil._e_difꜜ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜛ, matꜛ_ref, sun_geo.auxil._e_difꜛ_sife);
        mul!(sun_geo.auxil._e_difꜛ_sifꜜ, matꜜ_ref, sun_geo.auxil._e_difꜛ_sife);

        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_dirꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜜ_sifꜜ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜛ);
        photon_to_energy!(SPECTRA.Λ_SIF, sun_geo.auxil._e_difꜛ_sifꜜ);

        ϕ_sunlit = sen_geo.auxil.ϕ_f_sunlit[irt];
        ϕ_shaded = sen_geo.auxil.ϕ_f_shaded[irt];
        sl_1_ = local_lidf_weight(ϕ_sunlit, 1);
        sh_1_ = local_lidf_weight(ϕ_shaded, 1);
        sl_θ² = local_lidf_weight(ϕ_sunlit, _COS²_Θ_INCL_AZI);
        sh_θ² = local_lidf_weight(ϕ_shaded, _COS²_Θ_INCL_AZI);
        sl_S_ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_abs);
        sl_sθ = local_lidf_weight(ϕ_sunlit, sun_geo.auxil.fs_cos_incl);

        sun_geo.auxil._sif_sunlitꜛ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;
        sun_geo.auxil._sif_sunlitꜜ_dif .= sun_geo.auxil._e_difꜜ_sifꜛ .* sl_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sl_θ² .+
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sl_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sl_θ²;
        sun_geo.auxil._sif_sunlitꜛ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .+ sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;
        sun_geo.auxil._sif_sunlitꜜ_dir .= sun_geo.auxil._e_dirꜜ_sifꜛ .* sl_S_ .- sun_geo.auxil._e_dirꜜ_sifꜜ .* sl_sθ;
        sun_geo.auxil._sif_shadedꜛ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;
        sun_geo.auxil._sif_shadedꜜ     .= sun_geo.auxil._e_difꜜ_sifꜛ .* sh_1_ .- sun_geo.auxil._e_difꜜ_sifꜜ .* sh_θ² .+
                                          sun_geo.auxil._e_difꜛ_sifꜛ .* sh_1_ .+ sun_geo.auxil._e_difꜛ_sifꜜ .* sh_θ²;

        ilai_direct = (1 - sun_geo.auxil.τ_ss_layer[irt]) / sun_geo.auxil.ks_leaf;
        ilai_diffuse = 1 - can_str.auxil.τ_dd_isotropic[irt];
        p_sun = sun_geo.auxil.p_sunlit[irt];
        e_down = sun_geo.auxil._sif_sunlitꜜ_dir .* ilai_direct .+
                 sun_geo.auxil._sif_sunlitꜜ_dif .* p_sun .* ilai_diffuse .+
                 sun_geo.auxil._sif_shadedꜜ .* (1 - p_sun) .* ilai_diffuse;
        e_up = sun_geo.auxil._sif_sunlitꜛ_dir .* ilai_direct .+
               sun_geo.auxil._sif_sunlitꜛ_dif .* p_sun .* ilai_diffuse .+
               sun_geo.auxil._sif_shadedꜛ .* (1 - p_sun) .* ilai_diffuse;
        e_ref[:, irt] .= e_down .+ e_up;
    end;

    return e_ref
end;

function setup_spac_for_sif(; air_bounds = collect(6.0:2.0:12.0))
    config = ENS.SPACConfig(FT; dataset = ENS.OLD_PHI_2021_1NM);
    @assert config.FEATURES.ENABLE_SIF;
    @assert config.METHODS.FLUORESCENCE_SPECTRA_METHOD isa ENS.PlatespectFluorescenceSpectra;

    spac = ENS.BulkSPAC(config; air_bounds = air_bounds);
    spac.meteo.rad_sw.e_dir .*= 800;
    spac.meteo.rad_sw.e_dif .*= 200;
    ESPAC.prescribe_traits!(config, spac; lai = 4.0, ci = 0.7, sai = 0.2);
    ESPAC.initialize_spac!(config, spac);
    for leaf in spac.plant.leaves
        leaf.photosystem.auxil.ϕ_f .= FT(0.02);
    end;

    CO.soil_albedo!(config, spac);
    CO.canopy_structure!(config, spac);
    CO.sun_geometry_aux!(config, spac);
    CO.sun_geometry!(config, spac);
    CO.longwave_radiation!(spac);
    CO.shortwave_radiation!(config, spac);
    CO.sensor_geometry_aux!(config, spac);
    CO.sensor_geometry!(config, spac);
    CO.reflection_spectrum!(config, spac);

    return config, spac
end;

@testset "e_sif_chl simplification equivalence" begin
    config, spac = setup_spac_for_sif();

    @testset "leaf matꜛ_chl equals (mat_b_chl + mat_f_chl) / 2" begin
        for leaf in spac.plant.leaves
            matꜛ_ref, matꜜ_ref = mat_chl_expanded_reference!(config, leaf.bio, spac.cache);
            @test isapprox(leaf.bio.auxil.matꜛ_chl, matꜛ_ref; rtol = RTOL, atol = ATOL);
            @test any(!iszero, matꜜ_ref);
        end;
    end;

    @testset "collapsed Step 0 matches expanded up+down sum" begin
        populate_phi_f!(config, spac);
        e_ref = e_sif_chl_expanded_reference!(config, spac);

        CO.fluorescence_spectrum!(config, spac);
        e_new = copy(spac.canopy.sun_geometry.auxil.e_sif_chl);

        @test isapprox(e_new, e_ref; rtol = RTOL, atol = ATOL);

        max_diff = maximum(abs.(e_new .- e_ref));
        rel_diff = max_diff / maximum(abs.(e_ref));
        println("Step 0 e_sif_chl max |Δ| = ", max_diff, "  max rel = ", rel_diff);
    end;
end;

println("All e_sif_chl equivalence tests passed.");
