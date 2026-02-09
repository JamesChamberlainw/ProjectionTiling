import ee
import time
import math
import re




def export_tiles_to_gee(tiles, assetId_prefix, i_offset=0, sleep_time=5):
    """
        Exports tiles to GEE asset (as Asset)

        tiles: list of ee.FeatureCollection tiles
        assetId_prefix: str, prefix for assetId (asset folder + filename prefix) e.g., "projects/jameswilliamchamberlain/assets/activechunks_reduced/chunk_"
        i_offset: int, offset for tile numbering
        sleep_time: int, time to wait between exports to avoid failures
    """

    if assetId_prefix is None:
        print("No assetId_prefix provided, skipping export.")
        return

    # upload each tile (and wait to avoid failures)
    for i, tile in enumerate(tiles):
        # chunk_name = f"chunk_{i+1}"
        task = ee.batch.Export.table.toAsset(
            collection=tile,
            description=f"export_{i_offset+i+1}",
            assetId=f"{assetId_prefix}{i_offset+i+1}"
        )

        task.start()
        
        time.sleep(sleep_time)




class projTiler():
    def __init__(self, polygon, projection, width_px, assetId_prefix=None, shpfile_max_polygons=71**2, export_fn=export_tiles_to_gee, identifier_col="identifier", identifier_name_fn=lambda col, row, lon_min, lon_max, lat_min, lat_max: f"col{col}_row{row}_lonmin{lon_min}__lonmax{lon_max}_latmin{lat_min}__latmax{lat_max}"):
        """
            Tiles a given polygon into tiles that algin with the provided projection

            args (required): 
                polygon: ee.Geometry.Polygon, the polygon to tile 
                projection: ee.Projection, the projection to align tiles with (e.g., extract from an ee.Image.projection()
                width_px: int, the width of each tile in pixels (in this case a "pixel" is defined within the provided projection)
            args (optional):
                assetId_prefix: str, prefix for assetId for export_fn (asset folder + filename prefix) e.g., "projects/jameswilliamchamberlain/assets/activechunks_reduced/chunk_"
                shpfile_max_polygons: int, maximum number of polygons per shapefile when splitting large polygons (default 71**2 ~ 5041)
                export_fn: function, function to export tiles (default is export_tiles_to_gee MUST accept positional args: tiles (lst of ee.FeatureCollection), assetId_prefix (str))
                identifier_name_fn: function, function to name each tile (unqiue identifier) (default is a lambda that takes in col, row, lon_min, lon_max, lat_min, lat_max)
                identifier_col: str, the name of the unique identifier column for each generated polyogon tile (default "identifier")
        """  
        
        self.polygon = polygon
        self.proj = projection
        self.width_px = width_px
        self.shpfile_max_polygons = shpfile_max_polygons
        self.name_fn = identifier_name_fn
        self.identifier_col = identifier_col
        
        self.crs = self.proj.crs().getInfo()
        info = self.proj.getInfo()
        T_raw = info.get("transform") or self.proj.transform().getInfo() # or proj.crsTransform().getInfo() # not stored in pixel but in meters on a and e
        self.T = self.normalize_transform(T_raw)
        scale = float(self.proj.nominalScale().getInfo())
        a,b,c,d,e,f = self.T
        self.T = (a/scale, b/scale, c, d/scale, e/scale, f)      # rescale down to pixels  
        self.scale = scale

        self.tile(width_px, shpfile_max_polygons=shpfile_max_polygons)
        # export_fn(self.tile_lst, assetId_prefix="projects/jameswilliamchamberlain/assets/activechunks_reduced/chunk_", i_offset=0, sleep_time=5)
        export_fn(self.tile_lst, assetId_prefix)

    def normalize_transform(self, T):
        """
            Normalise GEE transfomrs to GDAL / rasterio affine transform format (a,b,c,d,e,f) 
        """
        if isinstance(T, (list, tuple)) and len(T) == 6 and not isinstance(T[0], (list, tuple)):
            return tuple(map(float, T))
        if isinstance(T, (list, tuple)) and len(T) == 2 and all(len(r) == 3 for r in T):
            (a,b,c),(d,e,f)=T; return (float(a),float(b),float(c),float(d),float(e),float(f))
        if isinstance(T, (list, tuple)) and len(T) == 3 and all(len(r) == 3 for r in T):
            (a,b,c),(d,e,f),_ = T; return (float(a),float(b),float(c),float(d),float(e),float(f))
        if isinstance(T, (list, tuple)) and len(T) == 9:
            a,b,c,d,e,f, *_ = T; return (float(a),float(b),float(c),float(d),float(e),float(f))
        if isinstance(T, str) and "PARAM_MT" in T:
            vals = {}
            for name in ("A","B","C","D","E","F"):
                m = re.search(rf'PARAM_MT\["{name}"\s*,\s*([-+0-9.eE]+)\]', T); vals[name]=float(m.group(1))
            return (vals["A"],vals["B"],vals["C"],vals["D"],vals["E"],vals["F"])
        raise ValueError(f"Unrecognized transform format: {T}")

    def invert2x2(self, a, b, d, e):
        det=a*e-b*d
        if det==0: raise ValueError("Non-invertible affine")
        return (e/det, -b/det, -d/det, a/det)

    def world_to_pixel(self, T, x, y):
        a,b,c,d,e,f=T
        ia,ib,id_,ie=self.invert2x2(a,b,d,e)
        dx,dy=x-c,y-f
        return ia*dx+ib*dy, id_*dx+ie*dy

    def pixel_to_world(self, T, col, row):
        a,b,c,d,e,f=T
        return a*col+b*row+c, d*col+e*row+f

    def snap_point(self, T, x, y):
        col,row = self.world_to_pixel(T,x,y)
        return self.pixel_to_world(T, math.floor(col), math.floor(row)) # round or floor? both seem to work 

    def tile_from_pixel(self, T, col_ref, row_ref, width_px, height_px):
        c0, r0 = int(round(col_ref)), int(round(row_ref)) # upper-left  
        c1, r1 = c0+width_px, r0+height_px  # with default +width_px +height_px positions else it moves around the following is true: 
        x0,y0 = self.pixel_to_world(T,c0,r0) # c0 r0 - bottom left 
        x1,y1 = self.pixel_to_world(T,c1,r0) # c1 r0 - bottom right 
        x2,y2 = self.pixel_to_world(T,c1,r1) # c1 r1 - top right
        x3,y3 = self.pixel_to_world(T,c0,r1) # c0 r1 - top left
        return [(x0,y0),(x1,y1),(x2,y2),(x3,y3),(x0,y0)]

    def tile_from_world_point(self, T, x_ref, y_ref, width_px, height_px):
        xs,ys = self.snap_point(T,x_ref,y_ref)
        col,row = self.world_to_pixel(T,xs,ys)
        return self.tile_from_pixel(T,col,row,width_px,height_px)
    
    def create_bbox_coords(self, points):
        """
            Create bounding box pixel coordinates from world coordinates.
        """

        cols, rows = [], []
        for x, y in points:
            col, row = self.world_to_pixel(self.T, x, y)
            cols.append(col)
            rows.append(row)

        col_min, col_max = min(cols), max(cols)
        row_min, row_max = min(rows), max(rows)
        return col_min, col_max, row_min, row_max

    def tile_polygon(self, width_px, height_px, polygon=None, clip_polygon=None):
        """ Tile the polygon into tiles of size (width_px, height_px) aligned with transform T without Google Earth Engine. """

        # replacement of polygon (optional otherwise assumes self.polygon - useful for larger shapefiles)
        if polygon is None:
            polygon = self.polygon

        # replacement of clip polygon (optional otherwise assumes its the maximum outer bounds of the provided polygon)
        if clip_polygon is None:
            clip_polygon = polygon 

        bounds = polygon.transform(self.proj, 1).bounds(proj=self.proj, maxError=1)
        bounds_coords = bounds.coordinates().getInfo()[0]

        col_min, col_max, row_min, row_max = self.create_bbox_coords(bounds_coords)

        # print("Bounding Box in Pixel Coordinates:")
        # print("col_min:", col_min, "col_max:", col_max)
        # print("row_min:", row_min, "row_max:", row_max)
        # n_cols = math.ceil((col_max - col_min) / width_px)
        # n_rows = math.ceil((row_max - row_min) / height_px)
        # print("tiles to create:", n_cols * n_rows)

        features = []
        for col in range(int(col_min), int(col_max), int(width_px)):      # step by tile width
            for row in range(int(row_min), int(row_max), int(height_px)): # step by tile height
                tile_ring = self.tile_from_pixel(self.T, col, row, int(width_px), int(height_px))
                # IMPORTANT: pass a list of rings
                poly = ee.Geometry.Polygon([tile_ring], proj=self.proj, geodesic=False)
                # name polygon 
                # poly.set("identifier", self.name_fn(col, row, tile_ring[0][0], tile_ring[1][0], tile_ring[0][1], tile_ring[3][1]))
                features.append(ee.Feature(poly, {'col': col, 'row': row}).set(self.identifier_col, self.name_fn(abs(col), abs(row), tile_ring[0][0], tile_ring[1][0], tile_ring[0][1], tile_ring[3][1])))

        # ensure tiles is a FeatureCollection
        tiles = ee.FeatureCollection(features)

        # filter tiles to fit polygon
        tiles = tiles.filterBounds(clip_polygon)
        print("Tiles created:", tiles.size().getInfo())

        return tiles
    
    def tile(self, width_px, shpfile_max_polygons=5000):
        # split into reasonable polygons that would result in manageable tiles
        # files = []

        # find max number of width_px and height_px that can fit into shpfile_max_polygons  
        polygon = self.polygon
        shapefile_bounds = []

        max_width_px = math.sqrt(shpfile_max_polygons) * width_px # assume square-ish
        print("Max width px per shapefile:", max_width_px)
        print(shpfile_max_polygons, "max polygons per shapefile")
        print(width_px, "tile width px")
        bounds = self.tile_polygon(max_width_px, max_width_px)

        # create tiles for each given polygon in shapefile
        polygon_shapefile_lst = []
        bounds_list = bounds.toList(bounds.size())

        # return bounds
        for i in range(bounds.size().getInfo()):
            bound = bounds_list.get(i)
            polygon_shapefile_lst.append(ee.Feature(bound).geometry())

        tile_lst = []

        for i, poly in enumerate(polygon_shapefile_lst):
            tiles = self.tile_polygon(width_px, width_px, polygon=poly, clip_polygon=self.polygon)
            tile_lst.append(tiles)

        self.tile_lst = tile_lst
        self.bounds = bounds

        return tile_lst, bounds