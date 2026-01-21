import logging
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
import numpy as np
import tempfile
from pathlib import Path

import snail.intersection as snint

from . import utils
from .. import naming


def vectorised_damage_calculation(
    defended_values: np.ndarray,
    damage_function: callable,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns:
        damage_units: damaged units (binary: damaged area)
        damage_fraction: damage fraction for cost calculation
    """
    damage_frac = np.vectorize(damage_function)(defended_values)
    return damage_frac


def intersect(
        vector:gpd.GeoDataFrame,
        rasters:list[str],
        damage_curves:dict,
        rehab_costs:dict,
        design_standards:dict,
        splits_path: str = None,
    ) -> gpd.GeoDataFrame:
    """Main intersection function for point geometries."""

    assert len(rasters) > 0, "No rasters provided for intersection."

    logging.info("Constructing grid from rasters...")
    grid, window = utils.process_raster_grid(rasters, vector)
    raster_basenames = utils.make_raster_basenames(rasters)

    logging.info("Processing point geometries...")
    vector = vector.reset_index(drop=True)
    vector_splits = vector.copy()  # No splitting needed for points

    logging.info("Finding indices...")
    vector_splits = snint.apply_indices(
        vector_splits, grid, index_i="raster_i", index_j="raster_j"
    )

    logging.info("Setting all measurements to 1 unit...")
    vector_splits["unit"] = 1
    vector_splits["unit_type"] = "unit"

    # make a temporary multi-band VRT to do all raster reads in one go
    with tempfile.TemporaryDirectory() as temp_dir:
        vrt_path = utils.create_multiband_vrt(rasters, output_dir=temp_dir)
        vector_splits = utils.copy_raster_values_multiband(
            vector_splits, vrt_path, raster_basenames, window
        )

    asset_types = list(vector_splits["asset_type"].unique())
    hazard_cols = [f"hazard-{Path(r).stem}" for r in rasters]

    asset_type_damages = []
    damage_cols = set()
    cost_cols = set()

    for asset_type in (pbar_asset := tqdm(asset_types)):
        pbar_asset.set_postfix({'Processing asset type': asset_type})
        vector_asset = vector_splits[vector_splits["asset_type"] == asset_type].copy()
        
        new_columns = {}
        
        for hazard_col in (pbar_haz := tqdm(hazard_cols)):
            pbar_haz.set_postfix({'Processing hazard column': hazard_col})
            hazard = naming.get_hazard_from_colname(hazard_col)

            design_standard_df = design_standards[hazard]
            design_hazard: str = design_standard_df.loc[asset_type, "design_hazard"]

            defended_col = hazard_col.replace("hazard-", "defended-")
            if design_hazard is None or pd.isna(design_hazard):
                logging.warning(
                    f"\nNo design standard provided for asset type '{asset_type}' from hazard '{hazard}'. "
                    "Skipping subtraction.\n")
                new_columns[defended_col] = vector_asset[hazard_col]
            else:
                design_standard_col = "hazard-" + design_hazard
                if design_standard_col not in vector_asset.columns:
                    raise ValueError(
                        f"\nDesign standard hazard column '{design_standard_col}' not found "
                        f"for asset type '{asset_type}'.\n"
                    )
                thresholds = vector_asset[design_standard_col].values
                defended_values = (
                    vector_asset[hazard_col].values - thresholds
                ).clip(min=0.0)
                new_columns[defended_col] = defended_values
                logging.info(
                    f"\nDesign standards: subtracted '{design_standard_col}' from '{hazard_col}' "
                    f"for '{asset_type}'.\n")
                
            # start vectorized damage and cost calculations
            defended_array = new_columns[defended_col]
            
            for suffix in ["mean", "min", "max"]:
                damage_function = damage_curves[(hazard, asset_type)][suffix]
                damage_col = defended_col.replace("defended-", "damage-") + "_" + suffix
                
                damage_frac = vectorised_damage_calculation(
                    defended_array, damage_function
                )
                
                new_columns[damage_col] = damage_frac
                damage_cols.add(damage_col)
                
                cost = rehab_costs[hazard].loc[asset_type, f"{suffix}_cost_usd"]
                cost_col = damage_col.replace("damage-", "cost-") + "_" + suffix
                new_columns[cost_col] = cost * damage_frac
                cost_cols.add(cost_col)
        
        new_columns_df = pd.DataFrame(new_columns, index=vector_asset.index)
        vector_asset = pd.concat([vector_asset, new_columns_df], axis=1)
        asset_type_damages.append(vector_asset)
    
    vector = pd.concat(asset_type_damages, axis=0).set_index("id")

    if splits_path is not None:
        logging.warning(f"Splits saving not relevant for point geometries, ignoring {splits_path}")

    return vector
