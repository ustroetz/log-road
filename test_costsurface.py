from logroad import CostSurface as cs

if __name__ == "__main__":
    riverfn = 'testdata/rivers.shp'
    slopefn = 'testdata/Slope.tif'
     
    cs.main(riverfn, slopefn)


