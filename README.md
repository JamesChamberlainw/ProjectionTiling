# Projecting Tiling Tool  
`projTiler` Class  

This tool aids in the creation of projection-aligned tiles across an assigned region defined by an input polygon. Generated tiles can be exported for direct upload to Google Earth Engine (GEE) as a set of shapefiles. Multiple shapefiles are produced due to GEE upload and processing limitations with overly large files, and the tool serves as a workaround by performing locally the operations that would otherwise consistantly fail over some regions with GEE direclty (where using `ee.ReduceToVector`).

The aim of this tool, developed with [PlotToSat](https://github.com/Art-n-MathS/PlotToSat) in mind, is to create a set of polygons (tiles) whose edges align with the satellite’s native resolution, ensuring they do not overlap adjacent data points. This allows data summarised by PlotToSat to remain as close as possible to the raw values, avoiding blurring or edge bleeding caused by overlaps between a polygon and the satellite’s native sampling grid.

(GIF BEING CREATED TONIGHT (14/04/2026) - IF YOU SEE THIS MESSAGE COME BACK IN A SHORT TIME AND THIS WILL UPDATED)

## Example Files

Two example notebooks are provided:
- `tiler_example_basics.ipynb` -- contains the bare essentials for using the tool
- `tiler_interactive_example.ipynb` -- a more in-depth guide covering all inputs, including an interactive example using geemap

An additional example demonstrates one approach to handling and automating PlotToSat with multiple shapefiles using `papermill`. The corresponding `requirements.txt` and notebook can be found in `/PlotToSat_example/.`

# Installation Guide 

Install the required dependencies using:

```
pip install -r .\requirements.txt
```

or alternativly just:

```
pip install earthengine-api geemap geopandas ipywidgets notebook ipykernel
``` 

## Requirements

Earth Engine API — `ee`  
Default Python libraries — `re`, `math`, `time`
