# Per-pixel stem (bark) reflectance derived from CLM5 PFT fractions.
#
# Broadband VIS / NIR stem reflectance by CLM5 PFT index (1..17), where index 1
# is not_vegetated. Woody bark values are from Majasalmi & Bright (2019, GMD);
# grass/crop values are CLM5 Table 2.3.1 dead-herbaceous stem optics
# (see SIFEsc pft_stemrho_dictin.jl). Stem transmittance is 0 (bark ~ opaque).

const STEM_RHO_VIS_BY_IDX = [
    0.00,                          # 1  not_vegetated
    0.12, 0.12, 0.12,              # 2-4  needleleaf (conifer bark)
    0.21, 0.21, 0.21, 0.21, 0.21,  # 5-9  broadleaf trees (deciduous bark)
    0.21, 0.21, 0.21,              # 10-12 shrubs (deciduous bark, extrapolated)
    0.31, 0.31, 0.31,              # 13-15 c3 arctic grass, c3 grass, c4 grass
    0.31, 0.31,                    # 16-17 c3 crop, c3 irrigated
];

const STEM_RHO_NIR_BY_IDX = [
    0.00,                          # 1  not_vegetated
    0.36, 0.36, 0.36,              # 2-4  needleleaf (conifer bark)
    0.49, 0.49, 0.49, 0.49, 0.49,  # 5-9  broadleaf trees (deciduous bark)
    0.49, 0.49, 0.49,              # 10-12 shrubs (deciduous bark, extrapolated)
    0.53, 0.53, 0.53,              # 13-15 c3 arctic grass, c3 grass, c4 grass
    0.53, 0.53,                    # 16-17 c3 crop, c3 irrigated
];

# Default flat stem reflectance (model default ρ_STEM) returned when a pixel has
# no vegetated stem area to weight by.
const STEM_RHO_DEFAULT = 0.2;


"""

    stem_rho_from_pft(gmd::Union{Dict,OrderedDict})

Compute a per-pixel broadband stem reflectance pair `(ρ_vis, ρ_nir)` from the PFT
fractions in the grid dictionary, given
- `gmd` Dictionary of GriddingMachine data in a grid

The pixel value is a stem-area-weighted mean of per-PFT stem reflectance. Because
ρ is an intensive optical property, each PFT is weighted by its relative stem area
`w[p] = PFT_FRACTIONS[p] * SAI_LAI[p]` (the shared pixel LAI anchor cancels in the
normalized ratio). Returns `(0.2, 0.2)` for non-vegetated grids or pixels with no
woody/herbaceous stem area.

"""
function stem_rho_from_pft(gmd::Union{Dict,OrderedDict})
    pfts = gmd["PFT_FRACTIONS"];

    # guard non-vegetated / soil-stub grids (PFT_FRACTIONS not the full 17-vector)
    if length(pfts) < 17
        return (STEM_RHO_DEFAULT, STEM_RHO_DEFAULT)
    end;

    w_sum = 0.0;
    vis_sum = 0.0;
    nir_sum = 0.0;
    for p in 2:17
        w = pfts[p] * SAI_LAI_BY_IDX[p];
        w_sum += w;
        vis_sum += w * STEM_RHO_VIS_BY_IDX[p];
        nir_sum += w * STEM_RHO_NIR_BY_IDX[p];
    end;

    if !(w_sum > 0)
        return (STEM_RHO_DEFAULT, STEM_RHO_DEFAULT)
    end;

    return (vis_sum / w_sum, nir_sum / w_sum)
end;


"""

    grid_stem_rho(gmd::Union{Dict,OrderedDict})

Resolve the per-pixel stem reflectance pair `(ρ_vis, ρ_nir)` for a grid.
- Use stored `gmd["STEM_RHO_VIS"]` / `gmd["STEM_RHO_NIR"]` when both exist as finite scalars.
- Otherwise, fall back to `stem_rho_from_pft(gmd)`.

This keeps compatibility with old setup files that do not carry stem reflectance
traits, while allowing new setup JLD2 files to carry precomputed values.

"""
function grid_stem_rho(gmd::Union{Dict,OrderedDict})
    if haskey(gmd, "STEM_RHO_VIS") && haskey(gmd, "STEM_RHO_NIR")
        vis = gmd["STEM_RHO_VIS"];
        nir = gmd["STEM_RHO_NIR"];
        if vis isa Number && nir isa Number && isfinite(vis) && isfinite(nir)
            return (vis, nir)
        end;
    end;

    return stem_rho_from_pft(gmd)
end;
