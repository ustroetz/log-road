from logroad import RoadModel as rm

if __name__ == "__main__":
    standsfn = 'testdata/test_stand.shp'       # EPSG3857
    costSurfacefn = 'testdata/CostSurface.tif' # EPSG32610
    newRoadsfn = 'testdata/newRoad.shp'
    
    rm.main(standsfn,costSurfacefn,newRoadsfn)


