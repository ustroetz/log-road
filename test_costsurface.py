from logroad import CostSurface as cs

if __name__ == "__main__":
    riverfn = 'testdata/rivers.shp' # EPSG32610
    slopefn = 'testdata/Slope.tif'  # EPSG32610
     
    cs.main(riverfn, slopefn)


