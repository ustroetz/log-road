from logroad import RoadModel as rm

if __name__ == "__main__":
    standsfn = 'testdata/test_stand.shp'
    costSurfacefn = 'testdata/CostSurface.tif'
    newRoadsfn = 'testdata/newRoad.shp'
    
    rm.main(standsfn,costSurfacefn,newRoadsfn)


