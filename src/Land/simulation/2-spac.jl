"""

    site_spac(config::SPACConfig{FT}, gmd::Union{Dict,OrderedDict}; lai_layer_strategy::Bool = false) where {FT}

Create a un-initialized SPAC using the data from a grid (CHL, VCMAX25, LAI, and CI are not prescribed as these changes with time), given
- `config` Configurations for SPAC
- `gmd` Dictionary of GriddingMachine data in a grid
- `c3c4` String to specify whether the leaf is C3 or C4, default to "C3"
- `lai_layer_strategy` Whether to use LAI-based air layer strategy (at least 10 layers, more for dense canopies; default: false)
- `maxlai_nlayer` Whether to scale canopy resolution with MAX_LAI via `n_layer = max(10, ceil(MAX_LAI)*10)` — i.e. a minimum of 10 layers for MAX_LAI < 1, and `ceil(MAX_LAI)*10` layers otherwise (10, 20, 30, ...; default: false). Takes precedence over `lai_layer_strategy` when both are true.

"""
function site_spac(config::SPACConfig{FT}, gmd::Union{Dict,OrderedDict}; c3c4::String = "C3", lai_layer_strategy::Bool = false, maxlai_nlayer::Bool = false) where {FT}
    # compute air layer bounds based on maximum LAI
    zc = max(FT(0.05), gmd["CANOPY_HEIGHT"]);
    air_bounds = if maxlai_nlayer
        # minimum 10 layers (MAX_LAI < 1); otherwise ceil(MAX_LAI)*10 layers
        n_layer = max(10, Int(ceil(maximum(gmd["LAI"]))) * 10);
        collect(0:(2n_layer+2)) * zc / (2n_layer)
    elseif lai_layer_strategy
        n_layer = max(10, Int(ceil(maximum(gmd["LAI"]) * 10)));
        collect(0:(2n_layer+2)) * zc / (2n_layer)
    else
        collect(0:21) * zc / 20
    end;

    spac = BulkSPAC(
                config;
                air_bounds = air_bounds,
                c3c4 = c3c4,
                elevation = gmd["ELEVATION"],
                latitude = gmd["LATITUDE"],
                longitude = gmd["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, zc/2, zc]);

    # initialize soil color
    spac.soil_bulk.trait.color = gmd["SOIL_COLOR"];

    # set up LMA and infinite carbon pool (to avoid NSC pool)
    for i in eachindex(spac.plant.leaves)
        spac.plant.leaves[i].bio.trait.lma = gmd["LMA"];
    end;
    config.FEATURES.UNLIMITED_NSC_POOL ? spac.plant.pool.c_pool = Inf : nothing;

    # set up SAI (per-pixel from PFT+LAI when enabled, else 0)
    init_sai = config.FEATURES.ENABLE_SAI ? grid_sai(gmd)[1] : 0;
    prescribe_traits!(config, spac; lai = gmd["LAI"][1], sai = init_sai);

    # set up per-pixel stem reflectance (VIS/NIR) when SAI is enabled; otherwise
    # there are no stems, so leave the model default flat ρ_STEM untouched
    if config.FEATURES.ENABLE_SAI
        rho_vis, rho_nir = grid_stem_rho(gmd);
        spectra = config.CONSTANTS.SPECTRA;
        for i in eachindex(spectra.Λ)
            spectra.ρ_STEM[i] = spectra.Λ[i] < 700 ? FT(rho_vis) : FT(rho_nir);
        end;
    end;

    # update soil type information per layer
    for i in eachindex(spac.soils)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gmd["SOIL_α"][i]) && !isnan(gmd["SOIL_N"][i]) && !isnan(gmd["SOIL_ΘR"][i]) && !isnan(gmd["SOIL_ΘS"][i])
            spac.soils[i].trait.vc.α = gmd["SOIL_α"][i];
            spac.soils[i].trait.vc.N = gmd["SOIL_N"][i];
            spac.soils[i].trait.vc.M = 1 - 1 / spac.soils[i].trait.vc.N;
            spac.soils[i].trait.vc.Θ_RES = gmd["SOIL_ΘR"][i];
            spac.soils[i].trait.vc.Θ_SAT = gmd["SOIL_ΘS"][i];
            spac.soils[i].state.θ = spac.soils[i].trait.vc.Θ_SAT;
        end;
    end;

    # initialize the SPAC
    initialize_spac!(config, spac);

    return spac
end;
