import os
import geopandas as gpd
import logging
from tqdm import tqdm


__all__ = ["format_subregion_name", "undo_subregion_formatting"]


def load_subregions(path):
    print(f"Loading subregions from: {path}.")
    with open(path, "r") as f:
        subregions = [line.strip() for line in f.readlines()]
    return subregions


def check_geoms(edges):
    geom_types = edges.geometry.geom_type.unique().tolist()
    assert all([geom_type in ["LineString"] for geom_type in geom_types]), \
        f"Found unexpected geometry types in edges: {geom_types}."


def calculate_overlap(row, admin:dict):
    clipped = row.geometry.intersection(
        admin[row["subregion"]]
    )
    return clipped.length, clipped


def prepare_admin_data(admin, subregions):
    admin = admin[["subregion", "geometry"]].copy()

    print("Checking admin subregions against expected list.")
    admin_subregions = set(admin["subregion"].unique().tolist())
    missing_subregions = admin_subregions - set(subregions)
    assert len(missing_subregions) == 0, \
        f"Admin data contains unexpected subregions: {missing_subregions}."

    admin = admin.dissolve(by="subregion", as_index=False)
    return admin


def check_for_duplicates(gdf, column="id"):
    n = len(gdf)
    n_unique = gdf[column].nunique()
    assert n == n_unique, f"Found {n - n_unique} duplicate indices in GeoDataFrame."


def intersect_by_overlap(edges, admin):
    check_for_duplicates(edges, column="id")
    check_for_duplicates(admin, column="subregion")

    edges_with_admin = gpd.sjoin(
        edges, admin, how="inner", predicate="intersects"
    )
    admin_dict = dict(zip(admin['subregion'], admin.geometry))
    tqdm.pandas(desc="Calculating overlaps")
    edges_with_admin[["overlap", "geometry"]] = edges_with_admin.progress_apply(
        calculate_overlap, axis=1, admin=admin_dict, result_type="expand"
    )
    edges_with_admin = edges_with_admin[edges_with_admin["overlap"] > 0].copy()
    edges_with_admin = edges_with_admin.explode(index_parts=True)
    edges_with_admin["id"] = edges_with_admin["id"] + "_" + edges_with_admin.index.get_level_values(1).astype(str)
    edges_with_admin = edges_with_admin.reset_index(drop=True)
    return edges_with_admin


def main(input, output, params):
    edges = gpd.read_parquet(input.edges).to_crs(params.local_crs)
    admin = gpd.read_parquet(input.admin).to_crs(params.local_crs)
    subregions = load_subregions(input.subregions)
    logging.info(f"Using local projection EPSG:{params.local_crs}.")

    check_geoms(edges)

    admin = prepare_admin_data(admin, subregions)
    edges = intersect_by_overlap(edges, admin)

    logging.info(f"Found {len(subregions)} subregions.")

    os.makedirs(output.edgedir, exist_ok=True)
    logging.info(f"Saving subregions to: {output.edgedir}.")

    empty = []
    for subregion in (pbar := tqdm(subregions)):
        pbar.set_postfix({'subregion': subregion})

        edges_subregion = edges[edges["subregion"] == subregion].copy()
        edges_subregion = edges_subregion.drop(columns=["index_right", "subregion"])
        assert "asset_type" in edges_subregion.columns, \
            "All assets must have an asset_type column."

        if len(edges_subregion) == 0:
            empty.append(subregion)

        edges_subregion = edges_subregion.to_crs("EPSG:4326")
        edges_subregion.to_parquet(
            os.path.join(output.edgedir, f"{subregion}.geoparquet"),
            index=False
        )
    
    if len(empty) > 0:
        logging.warning(f"No edges found for subregions: {*empty,}")
    logging.info("Done. Saved output with EPSG:4326 projection.")


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO
    )

    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    main(input, output, params)
