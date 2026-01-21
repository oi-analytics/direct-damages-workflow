checkpoint determine_subregions:
    """
    snakemake --cores 1 $DATADIR/assets/subregions.txt
    """
    input:
        subregions=f"{DATADIR}/admin/level{ADMIN_LEVEL}.gpkg"
    output:
        subregions=f"{DATADIR}/config/subregions.txt",
    run:
        import geopandas as geopandas
        
        subregions = geopandas.read_file(input.subregions)
        subregions = subregions["subregion"].unique().tolist()
        with open(output.subregions, "w") as f:
            f.write('\n'.join(subregions))


rule process_input_edges:
    """
    snakemake --cores 4 $DATADIR/assets/tza_roads_edges
    snakemake --cores 4 $DATADIR/assets/tza_railway_edges
    """
    input:
        edges=f"{DATADIR}/assets/raw/{{asset}}_edges.geoparquet",
        admin=f"{DATADIR}/admin/level{ADMIN_LEVEL}.gpkg",
        subregions=f"{DATADIR}/config/subregions.txt"
    output:
        edgedir=directory(f"{DATADIR}/assets/{{asset}}_edges"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_edges.py"


rule process_input_nodes:
    """
    snakemake --cores 4 $DATADIR/assets/tza_roads_bridges_and_culverts_nodes
    """
    input:
        points=f"{DATADIR}/assets/raw/{{asset}}_nodes.geoparquet",
        admin=f"{DATADIR}/admin/level{ADMIN_LEVEL}.gpkg",
        subregions=f"{DATADIR}/config/subregions.txt"
    output:
        pointdir=directory(f"{DATADIR}/assets/{{asset}}_nodes"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_nodes.py"


rule process_input_polygons:
    """
    snakemake --cores 4 $DATADIR/assets/tza_airports_polygons
    snakemake --cores 4 $DATADIR/assets/tza_iww_ports_polygons
    snakemake --cores 4 $DATADIR/assets/tza_maritime_ports_polygons
    """
    input:
        polys=f"{DATADIR}/assets/raw/{{asset}}_polygons.geoparquet",
        admin=f"{DATADIR}/admin/level{ADMIN_LEVEL}.gpkg",
        subregions=f"{DATADIR}/config/subregions.txt"
    output:
        polydir=directory(f"{DATADIR}/assets/{{asset}}_polygons"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_polygons.py"


rule process_all_assets:
    input:
        expand(f"{DATADIR}/assets/{{asset_geom}}",
            asset_geom=ASSET_GEOMS
        )
