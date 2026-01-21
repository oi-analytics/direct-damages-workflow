import logging
from pathlib import Path
from osgeo import gdal
from tqdm import tqdm
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
import snail.intersection as snint


def make_raster_basenames(raster_files):
    raster_basenames = []
    for raster_path in raster_files:
        basename = Path(raster_path).stem
        raster_basenames.append(basename)
    return raster_basenames


def grid_from_window(raster_file, bounds, verbose=False) -> snint.GridDefinition:
    """Create a snint.GridDefinition.from_raster for window defined by bounds."""
    with rasterio.open(raster_file) as src:
        window = rasterio.windows.from_bounds(
            bounds[0], bounds[1], bounds[2], bounds[3],
            transform=src.transform
        ).round()
        logging.info(f"Computed window from bounds: {window}")
        window_transform = rasterio.windows.transform(window, src.transform)

    grid = snint.GridDefinition(
        width=int(window.width),
        height=int(window.height),
        transform=window_transform,
        crs=src.crs.to_string()
    )
    return grid, window


def process_raster_grid(
        raster_files:list[str], vector:gpd.GeoDataFrame, verify_consistency=False
        ) -> snint.GridDefinition:
    """Make a grid for list of rasters, based on vector bounds."""
    bounds = vector.total_bounds
    grid, window = grid_from_window(raster_files[0], bounds)
    logging.info(f"{grid=}")

    if len(raster_files) > 1 and verify_consistency:
        logging.info("Checking raster grid consistency")
        for raster_path in raster_files[1:]:
            other_grid, _ = grid_from_window(raster_path, bounds)
            if other_grid != grid:
                raise AttributeError(
                    (
                        f"Raster attribute mismatch in file {raster_path}:\n"
                        f"Height: expected={grid.height}; actual={other_grid.height}\n"
                        f"Width: expected={grid.width}; actual={other_grid.width}\n"
                        f"Transform equal? {other_grid.transform == grid.transform}\n"
                        f"Transform expected= {grid.transform}\n"
                        f"Transform actual= {other_grid.transform}\n"
                        f"CRS equal? {other_grid.crs == grid.crs}"
                    )
                )
    
    return grid, window


def create_multiband_vrt(raster_files: list[str], output_dir: str = None):
    """
    Create a multi-band VRT from a list of raster files.
    """
    if output_dir is None:
        logging.warning("No output directory specified, using current directory.")
        output_dir = "."

    output_vrt = Path(output_dir) / "hazard_stack.vrt"
    vrt_options = gdal.BuildVRTOptions(separate=True)

    _vrt = gdal.BuildVRT(str(output_vrt), raster_files, options=vrt_options)
    _vrt = None
    
    logging.info(f"Created multi-band VRT with {len(raster_files)} bands at {output_vrt}")
    return str(output_vrt)


def copy_raster_values_multiband(
    vector_splits: gpd.GeoDataFrame, 
    vrt_path: str,
    raster_basenames: list[str],
    window: rasterio.windows.Window = None
) -> gpd.GeoDataFrame:
    """Copy windowed raster values from multi-band VRT to split geometries."""

    logging.info("Reading all raster values from multi-band VRT")
    
    raster_data: dict[str, np.ndarray] = {}
    
    with rasterio.open(vrt_path) as src:
        if window is not None:
            logging.info(f"Reading windowed area: {window}")
            data = src.read(window=window, masked=True)
        else:
            logging.info(f"Reading full raster: {src.height}x{src.width}")
            data = src.read(masked=True)
        
        logging.info(f"Read {data.shape[0]} bands from VRT")
        
        for band_idx, basename in enumerate(tqdm(raster_basenames, desc="Extracting values")):
            colname = f"hazard-{basename}"
            band_data = data[band_idx]
            raster_data[colname] = snint.get_raster_values_for_splits(
                vector_splits, band_data, index_i="raster_i", index_j="raster_j"
            )
    
    raster_data = pd.DataFrame(raster_data)
    vector_splits = pd.concat([vector_splits, raster_data], axis="columns")
    assert len(raster_data) == len(vector_splits)
    return vector_splits