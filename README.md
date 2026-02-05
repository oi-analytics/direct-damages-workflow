# Snakemake direct damages workflow

This is a simple snakemake workflow to carry out geospatial hazard–asset intersections and calculate various asset-level metrics, such as:
- maximum hazard intensity;
- units of asset damages; and
- rehabilitation cost of asset repair.

It has strict rules about the format of input files (hazard rasters, asset vectors, and management information like damage curves, rehabilitation costs, and design standards). Details of these rules are documented [below](#data-formatting-rules) (though the documentation may have some gaps, this is a WIP).

Once all the data is in the correct format, the [Quickstart](#quickstart) instructions describe how to run an analysis.

#### Other notes 

- This a WIP and has only been tested on Tanzania 2025. Some parameters for Tanzania might still be hardcoded.
- A sister repository for mobility modelling with the radiation model is in progress.
- This is currently a snakemake workflow rather than an importable module. This seemed most appropriate for the memory intensity, large number of I/O ops, and creation of many intermediate files.
- Minimal tests have been done, except to verify the intersection results.

## Quickstart

1. Clone the repository:
```bash
git clone git@github.com:oi-analytics/direct-damages-workflow.git
```

2. Create and activate the conda environment:
```bash
cd direct-damages-workflow
conda env create -f environment.yaml
conda activate oia-direct-damages
```

1. Download the input data from the shared drive (<add link>) and store it locally.

2. In `workflow/config.yaml`, set the `inputs` key to your input data location. The input data structure should be:

    - Vector asset files: `{inputs}/assets/raw/{asset}_{geom}.parquet`
    - Admin boundary files: `{inputs}/admin/level01.gpkg`
    - Hazard raster files: `{inputs}/hazards/raw/{hazard}_{epoch}_{scenario}_rp{rp}.tif`
    See the [file tree template](#filetree-in-datadir).

3. Run all intersections and damage calculations for the pluvial flood hazard on the Tanzania road edges asset (add `-n` flag for dry run):
```bash
cd workflow
export DATADIR=<path/to/project/datadir>
snakemake --cores 4 $DATADIR$/flags/tza_railway_edges/pluvial/.processed -n
```

---
# Data formatting rules

### Assets

Assets should be placed in `{input}/assets/raw` with format: `{asset}_{geom}.parquet` where:
- `asset`: asset type (e.g., `tza_airports`, `tza_railway`)
- `geom`: geometry type (e.g., `nodes`, `edges`, `polygon`)
- every file has as `asset_type` field that links entires to the csv files in [Management data](#management-data)

Rules in `rules/assets.smk` pre-process these to:
- Split by subregion
- Have a single geometry type per asset (e.g., LineString, Polygon, Point—no MultiLineStrings)
- Have WGS84 projection
- Have five columns: (unique) `id`, `asset_type`, `unit`, `unit_type`, `geometry`.

Processed assets are placed in `results/assets/` with format: `{asset}_{geom}/{subregion}.geoparquet`.

### Hazards

Hazard rasters should be placed in `results/hazards/raw/` with format:

```
{hazard}_{epoch}_{scenario}_rp{rp}.tif
```
- `hazard`: hazard name (e.g., `pluvial`, `cyclone`). Should match the hazard used in damage curves and rehabilitation costs
- `epoch`: time horizon (e.g., `2020`, `2050`, `2080`)
- `scenario`: climate scenario (e.g., `historical`, `ssp245`, `ssp585`)
- `rp`: return period (e.g., `00010`, `00050`, `00100`)

Hazard rasters should have WGS-84 projection and proper `NoData` values defined. Examples of preparing e.g., Fathom data will be put in the 2025 Tanzania project repo. 

Before the workflow begins, all hazards are aligned to a common grid and placed in `results/hazards/aligned/`. Rule: `rules/hazards.smk`.

### Management data

The workflow uses damage curves, design standards, and rehabilitation costs stored in `config/damage_curves`, `config/design_standards`, and `config/rehab_costs`:

| Data type            | Format | Location                                           |
|----------------------|--------|----------------------------------------------------|
| Damage curves        | CSV    | `config/damage_curves/{hazard}/{asset_type}.csv`   |
| Design standards     | CSV    | `config/design_standards/{hazard}.csv`             |
| Rehabilitation costs | CSV    | `config/rehab_costs/{hazard}.csv`                  |

Rehabilitation costs and design standards are indexed by `asset_type`.

Damage curves have an `intensity` column and three columns for damage fractions: `damage_fraction_max`, `damage_fraction_min`, `damage_fraction_mean`. Use `#` to note the units of intensity. Costs are specified in `costs_per_unit` with a separate `unit_type` column indicating `m` (LineStrings), `sqm` (Polygons), or `unit` (Points). Example processing scripts to get input costs and curves into the right format are in `analysis/scripts/` [check this link has probably changed].


## Filetree in datadir

### Before running workflow

```
.
├── admin
│   └── level01.gpkg
├── assets
│   └── raw
│       └── {iso}_{asset}_{geom}.geoparquet
├── config
│   ├── asset_types.csv
│   ├── damage_curves
│   │   └── {hazard}
│   │       └── {asset_type}.csv
│   ├── design_standards
│   │   └── {hazard}.csv
│   ├── rehab_costs
│   │   └── {hazard}.csv
│   └── subregions.txt
└── hazards
    └── raw
        ├── _reference.tif
        └── {hazard}_{epoch}_{scenario}_rp{rp}.tif
```

### After running workflow

```
.
├── admin
│   └── level01.gpkg
├── assets
│   ├── raw
│   │   └── {iso}_{asset}_{geom}.geoparquet
│   └── {iso}_{asset}_{geom}
│       └── {subregion}.geoparquet
├── config
│   ├── asset_types.csv
│   ├── damage_curves
│   │   └── {hazard}
│   │       └── {asset_type}.csv
│   ├── design_standards
│   │   └── {hazard}.csv
│   ├── rehab_costs
│   │   └── {hazard}.csv
│   └── subregions.txt
├── flags
│   └── {iso}_{asset}_{geom}
│       └── {hazard}
│           └── {subregion}.checked
├── hazards
│   ├── aligned
│   │   ├── _reference.tif
│   │   └── {hazard}_{epoch}_{scenario}_rp{rp}.tif
│   └── raw
│       ├── _reference.tif
│       └── {hazard}_{epoch}_{scenario}_rp{rp}.tif
└── intersections
    └── {iso}_{asset}_{geom}
        └── {hazard}
            └── {subregion}
                ├── profile.geoparquet
                └── splits.geoparquet
```

## To do

- [ ] Code to spatially aggregate and finalise outputs
- [ ] Code to make deliverable pivot tables of results
- [x] Separate hazard pre-processing into its own workflow
- [ ] (Optional) Add scripts for figures
- [ ] Option to interpolate to design standards
- [ ] More investigation simpson vs trapezoidal rule for expected value calculations
- [ ] Verify intersections.linestrings.unsplit() index matching logic
- [ ] Document the `exactextract` for polygons and add more damage metric functions.

