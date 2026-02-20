from pathlib import Path
from warnings import warn


def get_subregions():
    subregions_file = Path(f"{INPUTDIR}/config/subregions.txt")
    if not subregions_file.exists():
        return []
    with open(subregions_file) as f:
        subregions = [line.strip() for line in f if line.strip()]
    return subregions


rule intersect_subregion_hazard:
    """
    snakemake -c4 $INPUTDIR/intersections/tza_railway_edges/cyclone/dar_es_salaam/profile.geoparquet
    snakemake -c4 $INPUTDIR/intersections/tza_airports_polygons/cyclone/dar_es_salaam/profile.geoparquet
    snakemake -c4 $INPUTDIR/intersections/tza_roads_bridges_and_culverts_nodes/cyclone/dar_es_salaam/profile.geoparquet
    """
    input:
        asset_dir=f"{INPUTDIR}/assets/{{asset}}_{{geom}}",
        hazard_dir=f"{INPUTDIR}/hazards/aligned"
    output:
        vector=f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    params:
        hazard="{hazard}",
        subregion="{subregion}",
        copy_raster_values=True,
        crs=LOCAL_CRS,
        damage_curve_dir=f"{INPUTDIR}/config/damage_curves",
        rehab_cost_dir=f"{INPUTDIR}/config/rehab_costs",
        protection_dir=f"{INPUTDIR}/config/design_standards",
        splits_path=f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/splits.geoparquet"
    log:
        file="../logs/intersections/{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/intersections.py"


rule calculate_expected_values:
    """
    snakemake -c4 $INPUTDIR/intersections/tza_railway_edges/cyclone/dar_es_salaam/expected.parquet
    snakemake -c4 $INPUTDIR/intersections/tza_airports_polygons/cyclone/dar_es_salaam/expected.parquet
    snakemake -c4 $INPUTDIR/intersections/tza_roads_bridges_and_culverts_nodes/cyclone/dar_es_salaam/expected.parquet
    """
    input:
        vector=f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet"
    output:
        parquet=f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet"
    log:
        file="../logs/risk/expectations_{geom}_{asset}_{subregion}_{hazard}.log"
    script:
        "../scripts/expectations.py"


def all_subregions(wildcards):
    checkpoints.determine_subregions.get()
    subregions = get_subregions()
    return expand(
        f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/expected.parquet",
        geom=wildcards.geom,
        asset=wildcards.asset,
        hazard=wildcards.hazard,
        subregion=subregions
    )


rule all_results_for_asset_and_hazard:
    """
    snakemake -c4 $INPUTDIR/flags/tza_railway_edges/pluvial/.processed -n
    snakemake -c4 $INPUTDIR/flags/tza_roads_bridges_and_culverts_nodes/pluvial/.processed -n
    snakemake -c4 $INPUTDIR/flags/tza_airports_polygons/pluvial/.processed -n
    """
    input:
        all_subregions
    output:
        touch(f"{INPUTDIR}/flags/{{asset}}_{{geom}}/{{hazard}}/.processed")


rule all_intersections:
    """
    snakemake -c2 all_intersections -n
    """
    input:
        expand(
            f"{INPUTDIR}/intersections/{{asset_geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
            asset_geom=ASSET_GEOMS,
            hazard=HAZARDS,
            subregion=get_subregions()
        )


rule verify_intersections:
    """
    Check intersection results against input raster.

    snakemake -c1 $INPUTDIR/flags/tza_railway_edges/cyclone/dar_es_salaam.checked
    snakemake -c1 $INPUTDIR/flags/tza_roads_bridges_and_culverts_nodes/cyclone/dar_es_salaam.checked
    snakemake -c1 $INPUTDIR/flags/tza_airports_polygons/cyclone/dar_es_salaam.checked
    """
    input:
        vector=f"{INPUTDIR}/intersections/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}/profile.geoparquet",
        ref_dir=f"{INPUTDIR}/assets/{{asset}}_{{geom}}",
        hazdir=f"{INPUTDIR}/hazards/aligned"
    params:
        subregion="{subregion}",
        hazard="{hazard}"
    output:
        touch(f"{INPUTDIR}/flags/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}.checked")
    script:
        "../scripts/verify_intersections.py"