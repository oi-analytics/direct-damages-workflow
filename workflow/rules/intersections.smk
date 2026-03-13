from pathlib import Path
from warnings import warn


def get_subregions():
    subregions_file = Path(f"{TMPDIR}/config/subregions.txt")
    if not subregions_file.exists():
        return []
    with open(subregions_file) as f:
        subregions = [line.strip() for line in f if line.strip()]
    return subregions


rule intersect_subregion_hazard:
    """
    snakemake -c4 $RESDIR/intersections/tza_railway_edges/cyclone/dar_es_salaam/profile.geoparquet
    snakemake -c4 $RESDIR/intersections/tza_airports_polygons/cyclone/dar_es_salaam/profile.geoparquet
    snakemake -c4 $RESDIR/intersections/tza_roads_bridges_and_culverts_nodes/cyclone/dar_es_salaam/profile.geoparquet
    """
    input:
        asset_dir=f"{TMPDIR}/assets/{{asset}}_{{geom}}",
        hazard_dir=rules.align_hazard_rasters.output.outdir
    output:
        vector=f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    params:
        hazard="{hazard}",
        subregion="{subregion}",
        copy_raster_values=True,
        crs=LOCAL_CRS,
        damage_curve_dir=f"{INDIR}/config/damage_curves",
        rehab_cost_dir=f"{INDIR}/config/rehab_costs",
        protection_dir=f"{INDIR}/config/design_standards",
        splits_path=f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/splits.geoparquet"
    log:
        file="../logs/intersections/{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/intersections.py"


rule calculate_expected_values:
    """
    snakemake -c4 $RESDIR/intersections/tza_railway_edges/fluvial/kilimanjaro/expected.parquet
    """
    input:
        vector=f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    output:
        parquet=f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet"
    log:
        file="../logs/risk/expectations_{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/expectations.py"


def all_subregions(wildcards):
    checkpoints.determine_subregions.get()
    subregions = get_subregions()
    return expand(
        f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet",
        geom=wildcards.geom,
        asset=wildcards.asset,
        hazard=wildcards.hazard,
        subregion=subregions
    )


rule all_results_for_asset_and_hazard:
    """
    snakemake -c4 $TMPDIR/flags/tza_railway_edges/pluvial/.processed -n
    snakemake -c4 $TMPDIR/flags/tza_roads_bridges_and_culverts_nodes/pluvial/.processed -n
    snakemake -c4 $TMPDIR/flags/tza_airports_polygons/pluvial/.processed -n
    """
    input:
        all_subregions
    output:
        touch(f"{TMPDIR}/flags/{{asset}}_{{geom}}/{{hazard}}/.processed")


rule all_profiles:
    """
    snakemake -c2 all_profiles -n
    """
    input:
        lambda wildcards: expand(
            f"{RESDIR}/intersections/{{asset_geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
            asset_geom=ASSET_GEOMS,
            hazard=HAZARDS,
            subregion=(checkpoints.determine_subregions.get() or True) and get_subregions()
        )


rule all_expected:
    """
    snakemake -c2 all_expected -n
    """
    input:
        lambda wildcards: expand(
            f"{RESDIR}/intersections/{{asset_geom}}/{{hazard}}/{{subregion}}/expected.parquet",
            asset_geom=ASSET_GEOMS,
            hazard=HAZARDS,
            subregion=(checkpoints.determine_subregions.get() or True) and get_subregions()
        )


rule verify_intersections:
    """
    Check intersection results against input raster.

    snakemake -c1 $TMPDIR/flags/tza_railway_edges/cyclone/dar_es_salaam.checked
    snakemake -c1 $TMPDIR/flags/tza_roads_bridges_and_culverts_nodes/cyclone/dar_es_salaam.checked
    snakemake -c1 $TMPDIR/flags/tza_airports_polygons/cyclone/dar_es_salaam.checked
    """
    input:
        vector=f"{RESDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
        ref_dir=f"{TMPDIR}/assets/{{asset}}_{{geom}}",
        hazdir=f"{TMPDIR}/hazards"
    params:
        subregion="{subregion}",
        hazard="{hazard}"
    output:
        touch(f"{TMPDIR}/flags/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}.checked")
    script:
        "../scripts/verify_intersections.py"