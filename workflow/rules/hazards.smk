from pathlib import Path
import os


def get_all_input_hazards(wildcards):
    """Input function that runs at execution time"""
    hazards_dir = Path(f"{INPUTDIR}/hazards/raw")
    hazards = []
    for root, dirs, files in os.walk(hazards_dir):
        for file in files:
            if file.endswith(".tif"):
                if file.startswith("."):
                    continue
                hazards.append(os.path.join(root, file))
    if len(hazards) == 0:
        raise ValueError(f"No input hazard rasters found in {hazards_dir}")
    return hazards

    
rule align_hazard_rasters:
    """
    Align all hazard rasters to a common grid, corresponding to the
    highest resolution hazard data.

    snakemake -c4 align_hazard_rasters
    """
    input:
        reference_raster=f"{INPUTDIR}/hazards/raw/_reference.tif",
        rasters=get_all_input_hazards,
    output:
        outdir=directory(f"{INPUTDIR}/hazards/aligned")
    log:
        file="../logs/hazards/align_rasters.log"
    script:
        "../scripts/align_rasters.py"
