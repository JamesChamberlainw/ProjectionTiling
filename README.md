# Projecting Tiling Tool  
`projTiler` Class  

This tool takes an input polygon and projection (as a resolution) and creates projection-aligned tiles across the region. The tiles are exported as shapefiles and uploaded to a Google Earth Engine Asset.

### Additional Notes

This was originally designed with [PlotToSat](https://github.com/Art-n-MathS/PlotToSat) in mind. However, feel free to adapt this to your needs, as it is intentionally kept simple with adjustable parts for ease of reuse - see `tiler_interactive_example.ipynb`.

For larger regions, multiple shapefiles may be generated (depending on shapefile size limits and Earth Engine upload/processing limits). I recommend using `papermill` with a simple Earth Engine queue status checker to automate PlotToSat notebook execution. I can provide an example if you contact me.

## Requirements

Earth Engine API — `ee`  
Default Python libraries — `re`, `math`, `time`
