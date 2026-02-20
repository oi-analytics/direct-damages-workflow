checkpoint determine_subregions:
    """
    snakemake -c1 determine_subregions
    """
    input:
        subregions=f"{INPUTDIR}/admin/level{ADMIN_LEVEL}.parquet"
    output:
        subregions=f"{INPUTDIR}/config/subregions.txt",
    run:
        import geopandas as geopandas
        
        subregions = geopandas.read_parquet(input.subregions)
        subregions = subregions["subregion"].unique().tolist()
        with open(output.subregions, "w") as f:
            f.write('\n'.join(subregions))


rule process_input_edges:
    """
    snakemake --cores 4 $INPUTDIR/assets/tza_roads_edges
    snakemake --cores 4 $INPUTDIR/assets/tza_railway_edges
    """
    input:
        edges=f"{INPUTDIR}/assets/raw/{{asset}}_edges.geoparquet",
        admin=f"{INPUTDIR}/admin/level{ADMIN_LEVEL}.parquet",
        subregions=f"{INPUTDIR}/config/subregions.txt"
    output:
        edgedir=directory(f"{INPUTDIR}/assets/{{asset}}_edges"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_edges.py"


rule process_input_nodes:
    """
    snakemake --cores 4 $INPUTDIR/assets/tza_roads_bridges_and_culverts_nodes
    """
    input:
        points=f"{INPUTDIR}/assets/raw/{{asset}}_nodes.geoparquet",
        admin=f"{INPUTDIR}/admin/level{ADMIN_LEVEL}.parquet",
        subregions=f"{INPUTDIR}/config/subregions.txt"
    output:
        pointdir=directory(f"{INPUTDIR}/assets/{{asset}}_nodes"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_nodes.py"


rule process_input_polygons:
    """
    snakemake --cores 4 $INPUTDIR/assets/tza_airports_polygons
    snakemake --cores 4 $INPUTDIR/assets/tza_iww_ports_polygons
    snakemake --cores 4 $INPUTDIR/assets/tza_maritime_ports_polygons
    """
    input:
        polys=f"{INPUTDIR}/assets/raw/{{asset}}_polygons.geoparquet",
        admin=f"{INPUTDIR}/admin/level{ADMIN_LEVEL}.parquet",
        subregions=f"{INPUTDIR}/config/subregions.txt"
    output:
        polydir=directory(f"{INPUTDIR}/assets/{{asset}}_polygons"),
    params:
        local_crs=LOCAL_CRS
    script:
        "../scripts/process_input_polygons.py"


rule process_all_assets:
    input:
        expand(f"{INPUTDIR}/assets/{{asset_geom}}",
            asset_geom=ASSET_GEOMS
        )
