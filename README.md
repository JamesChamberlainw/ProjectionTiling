# Projecting Tiling Tool  
`projTiler` Class  

This tool aids in the creation of projection-aligned tiles across an assigned region defined by an input polygon. Generated tiles can be exported for direct upload to Google Earth Engine (GEE) as a set of shapefiles. Multiple shapefiles are produced due to GEE upload and processing limitations with overly large files, and the tool serves as a workaround by performing locally the operations that would otherwise consistantly fail over some regions with GEE direclty (where using `ee.ReduceToVector`).

The aim of this tool, developed with PlotToSat in mind, is to create a set of polygons (tiles) whose edges align with the satellite’s native resolution, ensuring they do not overlap adjacent data points. This allows data summarised by PlotToSat to remain as close as possible to the raw values, avoiding blurring or edge bleeding caused by overlaps between a polygon and the satellite’s native sampling grid.

(GIF BEING CREATED TONIGHT (14/04/2026) - IF YOU SEE THIS MESSAGE COME BACK IN A SHORT TIME AND THIS WILL UPDATED)

### Additional Notes

This was originally created with [PlotToSat](https://github.com/Art-n-MathS/PlotToSat) in mind. However, feel free to adapt this to your needs, as it is intentionally kept simple with adjustable parts for ease of reuse - see `tiler_interactive_example.ipynb`.

For larger regions, multiple shapefiles may be generated (depending on shapefile size limits and Earth Engine upload/processing limits). I recommend using `papermill` with a simple Earth Engine queue status checker to automate PlotToSat notebook execution. I can provide an example if you contact me.

## Requirements

Earth Engine API — `ee`  
Default Python libraries — `re`, `math`, `time`
