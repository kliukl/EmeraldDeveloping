# Per-pixel stem area index (SAI) derived from CLM5 PFT fractions and LAI.
#
# Woody-to-total area ratio (alpha = WAI / PAI) by CLM5 PFT index (1..17), where
# index 1 is not_vegetated. Forest and shrub alpha are from Fang et al. (2019,
# Rev. Geophys.), grass and crop alpha are expert-assigned (see SIFEsc pft_sai_dictin.jl).
# The SAI:LAI multiplier is alpha / (1 - alpha).

const SAI_ALPHA_BY_IDX = [
    0.00,                          # 1  not_vegetated
    0.21, 0.21, 0.21,              # 2-4  needleleaf evergreen temperate/boreal, needleleaf deciduous boreal
    0.10, 0.19, 0.19, 0.21, 0.21,  # 5-9  broadleaf evergreen tropical/temperate, deciduous tropical/temperate/boreal
    0.27, 0.27, 0.27,              # 10-12 evergreen shrub, deciduous temperate/boreal shrub
    0.05, 0.05, 0.05,              # 13-15 c3 arctic grass, c3 grass, c4 grass
    0.10, 0.10,                    # 16-17 c3 crop, c3 irrigated
];

const SAI_LAI_BY_IDX = [a == 0 ? 0.0 : a / (1 - a) for a in SAI_ALPHA_BY_IDX];

# CLM5 PFT index groups: forest/shrub get a constant SAI (anchored to peak LAI),
# grass/crop get a SAI that tracks LAI through the year.
const FOREST_SHRUB_IDX = 2:12;
const GRASS_CROP_IDX = 13:17;


"""

    sai_from_pft(gmd::Union{Dict,OrderedDict})

Compute a per-pixel daily SAI vector (same length as `gmd["LAI"]`) from the PFT fractions and LAI
already present in the grid dictionary, given
- `gmd` Dictionary of GriddingMachine data in a grid

The pixel SAI is a vegetated-fraction-weighted mean of per-PFT SAI:
- forest/shrub PFTs contribute a constant `SAI_LAI[p] * P90(LAI)` (P90 = 90th percentile of daily LAI);
- grass/crop PFTs contribute `SAI_LAI[p] * LAI(t)` (tracks LAI).

so `SAI(t) = A + B * LAI(t)`. Returns zeros for non-vegetated grids.

"""
function sai_from_pft(gmd::Union{Dict,OrderedDict})
    pfts = gmd["PFT_FRACTIONS"];
    lai = gmd["LAI"];

    # guard non-vegetated / soil-stub grids (PFT_FRACTIONS not the full 17-vector)
    if length(pfts) < 17
        return zeros(eltype(lai), length(lai))
    end;

    veg = sum(@view pfts[2:17]);
    if !(veg > 0)
        return zeros(eltype(lai), length(lai))
    end;

    w = pfts ./ veg;

    finite_lai = filter(!isnan, lai);
    p90 = isempty(finite_lai) ? zero(eltype(lai)) : quantile(finite_lai, 0.9);

    # constant (forest/shrub) and LAI-tracking (grass/crop) coefficients
    a_const = 0.0;
    for p in FOREST_SHRUB_IDX
        a_const += w[p] * SAI_LAI_BY_IDX[p];
    end;
    a_const *= p90;

    b_coef = 0.0;
    for p in GRASS_CROP_IDX
        b_coef += w[p] * SAI_LAI_BY_IDX[p];
    end;

    return a_const .+ b_coef .* lai
end;


"""

    grid_sai(gmd::Union{Dict,OrderedDict})

Resolve the daily SAI vector for a grid.
- Use `gmd["SAI"]` when it exists as a vector and its length matches `gmd["LAI"]`.
- Otherwise, fall back to `sai_from_pft(gmd)`.

This keeps compatibility with old setup files where GriddingMachine stores scalar
`"SAI" => 0`, while allowing new setup JLD2 files to carry precomputed daily SAI.

"""
function grid_sai(gmd::Union{Dict,OrderedDict})
    lai = gmd["LAI"];

    if haskey(gmd, "SAI")
        sai = gmd["SAI"];
        if sai isa AbstractVector && length(sai) == length(lai)
            return sai
        end;
    end;

    return sai_from_pft(gmd)
end;
