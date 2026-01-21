"""
Old rules to rework later.
"""
# rule check_asset_hazard_exposure:
#     """
#     Double-check intersection results against input raster.
# 
#     snakemake --cores 4 ../results/flags/tza_railway_edges/pluvial/kilimanjaro.checked
#     snakemake --cores 4 ../results/flags/tza_roads_bridges_and_culverts_nodes/pluvial/kilimanjaro.checked
#     snakemake --cores 4 ../results/flags/tza_airports_polygons/pluvial/kilimanjaro.checked
#     snakemake --cores 4 ../results/flags/tza_railway_edges/pluvial/shinyanga.checked
#     """
#     input:
#         vector="../results/risk/{asset}_{geom}/{hazard}/{subregion}/profile.geoparquet",
#         ref_dir="../results/assets/{asset}_{geom}",
#         hazdir="../results/hazards/aligned"
#     params:
#         subregion="{subregion}",
#         hazard="{hazard}"
#     output:
#         touch("../results/flags/{asset}_{geom}/{hazard}/{subregion}.checked")
#     script:
#         "../scripts/intersections_check.py"
# 
# def all_subregion_flags(wildcards):
#     checkpoints.determine_subregions.get()
#     subregions = get_subregions()
#     return expand(
#         f"{DATADIR}/flags/{{asset}}_{{geom}}/{{hazard}}/{{subregion}}.checked",
#         geom=wildcards.geom,
#         asset=wildcards.asset,
#         hazard=wildcards.hazard,
#         subregion=subregions
#     )


# # these are temporary workarounds for the Tanzania 2025 project
# rule calculate_expected_metrics_cleaned:
#     """
#     snakemake --cores 4 ../results/risk_final/tza_railway_edges/pluvial/kilimanjaro/expected.parquet
#     snakemake --cores 4 ../results/risk_final/tza_roads_bridges_and_culverts_nodes/pluvial/kilimanjaro/expected.parquet
#     snakemake --cores 4 ../results/risk_final/tza_airports_polygons/pluvial/kilimanjaro/expected.parquet
#     """
#     input:
#         vector="../results/risk_final/{asset}_{geom}/{hazard}/{subregion}/profile.geoparquet"
#     output:
#         parquet="../results/risk_final/{asset}_{geom}/{hazard}/{subregion}/expected.parquet"
#     log:
#         file="../logs/risk_final/expectations_{geom}_{asset}_{hazard}_{subregion}.log"
#     script:
#         "../scripts/expectations.py"
# 
# rule calculate_all_cleaned_metrics:
#     """
#     snakemake --cores 1 calculate_all_cleaned_metrics -n
#     """
#     input:
#         expand(
#             "../results/risk_final/{asset_geom}/{hazard}/{subregion}/expected.parquet",
#             asset_geom=ASSET_GEOMS,
#             hazard=[
#                 # "pluvial",
#                 # "fluvial",
#                 # "coastal",
#                 # "landslide",
#                 # "cyclone",
#                 # "hd35",
#                 "tasmax"
#                 ],
#             subregion=get_subregions()
#         )