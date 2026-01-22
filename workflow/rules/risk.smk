from pathlib import Path
from warnings import warn


def get_subregions():
    subregions_file = Path(f"{DATADIR}/config/subregions.txt")
    if not subregions_file.exists():
        return []
    with open(subregions_file) as f:
        subregions = [line.strip() for line in f if line.strip()]
    return subregions


rule intersect_subregion_hazard:
    """
    snakemake -c4 $DATADIR/intersections/tza_railway_edges/pluvial/kilimanjaro/profile.geoparquet
    snakemake -c4 $DATADIR/intersections/tza_airports_polygons/pluvial/kilimanjaro/profile.geoparquet
    snakemake -c4 $DATADIR/intersections/tza_roads_bridges_and_culverts_nodes/pluvial/kilimanjaro/profile.geoparquet

    snakemake -c4 $DATADIR/intersections/tza_railway_edges/hd35/kilimanjaro/profile.geoparquet
    snakemake -c4 $DATADIR/intersections/tza_airports_polygons/hd35/kilimanjaro/profile.geoparquet
    snakemake -c4 $DATADIR/intersections/tza_roads_bridges_and_culverts_nodes/hd35/kilimanjaro/profile.geoparquet
    """
    input:
        asset_dir=f"{DATADIR}/assets/{{asset}}_{{geom}}",
        hazard_dir=f"{DATADIR}/hazards/aligned"
    output:
        vector=f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    params:
        hazard="{hazard}",
        subregion="{subregion}",
        copy_raster_values=True,
        crs=LOCAL_CRS,
        damage_curve_dir=f"{DATADIR}/config/damage_curves",
        rehab_cost_dir=f"{DATADIR}/config/rehab_costs",
        protection_dir=f"{DATADIR}/config/design_standards",
        splits_path=f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/splits.geoparquet"
    log:
        file="../logs/intersections/{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/intersections.py"


rule calculate_expected_values:
    """
    snakemake -c4 $DATADIR/intersections/tza_railway_edges/pluvial/kilimanjaro/expected.parquet
    snakemake -c4 $DATADIR/intersections/tza_roads_bridges_and_culverts_nodes/pluvial/kilimanjaro/expected.parquet
    snakemake -c4 $DATADIR/intersections/tza_airports_polygons/pluvial/kilimanjaro/expected.parquet

    snakemake -c4 $DATADIR/intersections/tza_railway_edges/hd35/kilimanjaro/expected.parquet
    snakemake -c4 $DATADIR/intersections/tza_roads_bridges_and_culverts_nodes/hd35/kilimanjaro/expected.parquet
    snakemake -c4 $DATADIR/intersections/tza_airports_polygons/hd35/kilimanjaro/expected.parquet
    """
    input:
        vector=f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    output:
        parquet=f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet"
    log:
        file="../logs/risk/expectations_{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/expectations.py"


def all_subregions(wildcards):
    checkpoints.determine_subregions.get()
    subregions = get_subregions()
    return expand(
        f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet",
        geom=wildcards.geom,
        asset=wildcards.asset,
        hazard=wildcards.hazard,
        subregion=subregions
    )


rule all_results_for_asset_and_hazard:
    """
    snakemake -c4 $DATADIR/flags/tza_railway_edges/pluvial/.processed -n
    snakemake -c4 $DATADIR/flags/tza_roads_bridges_and_culverts_nodes/pluvial/.processed -n
    snakemake -c4 $DATADIR/flags/tza_airports_polygons/pluvial/.processed -n
    """
    input:
        all_subregions
    output:
        touch(f"{DATADIR}/flags/{{asset}}_{{geom}}/{{hazard}}/.processed")


rule all_intersections:
    """
    snakemake -c2 all_intersections -n
    """
    input:
        expand(
            f"{DATADIR}/intersections/{{asset_geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
            asset_geom=ASSET_GEOMS,
            hazard=HAZARDS,
            subregion=get_subregions()
        )


rule verify_intersections:
    """
    Check intersection results against input raster.

    snakemake -c1 $DATADIR/flags/tza_railway_edges/pluvial/kilimanjaro.checked
    snakemake -c1 $DATADIR/flags/tza_roads_bridges_and_culverts_nodes/pluvial/kilimanjaro.checked
    snakemake -c1 $DATADIR/flags/tza_airports_polygons/pluvial/kilimanjaro.checked
    snakemake -c1 $DATADIR/flags/tza_railway_edges/pluvial/shinyanga.checked
    """
    input:
        vector=f"{DATADIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
        ref_dir=f"{DATADIR}/assets/{{asset}}_{{geom}}",
        hazdir=f"{DATADIR}/hazards/aligned"
    params:
        subregion="{subregion}",
        hazard="{hazard}"
    output:
        touch(f"{DATADIR}/flags/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}.checked")
    script:
        "../scripts/verify_intersections.py"