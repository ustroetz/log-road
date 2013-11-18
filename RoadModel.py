import ogr, gdal, osr
import os
import numpy as np
from skimage.graph import route_through_array
import requests

def buffer(standsfn, bufferfn, Dist=250):
    stands = ogr.Open(standsfn)
    lyrStands = stands.GetLayer()
    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(bufferfn):
        shpdriver.DeleteDataSource(bufferfn)
    ds = shpdriver.CreateDataSource(bufferfn)
    lyrStandsBuffer = ds.CreateLayer(bufferfn, geom_type=ogr.wkbPolygon)
    featureDefn = lyrStandsBuffer.GetLayerDefn()

    featureStand = lyrStands.GetNextFeature()
    while featureStand:
        geomStand = featureStand.GetGeometryRef()
        geomBufferStand = geomStand.Buffer(Dist)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBufferStand)
        lyrStandsBuffer.CreateFeature(outFeature)
        outFeature.Destroy()
        featureStand.Destroy()
        featureStand = lyrStands.GetNextFeature()

def selectCell(gridfn,bufferfn):
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    gridList = range(lyrGrid.GetFeatureCount())
    
    buffer = ogr.Open(bufferfn)
    lyrBuffer = buffer.GetLayer()
    bufferList = range(lyrBuffer.GetFeatureCount())
    
    standDict = {}
    
    for i in gridList:
        # loop grid
        featGrid = lyrGrid.GetFeature(i)
        geomGrid = featGrid.GetGeometryRef()
        
        cellStands = []
        for j in bufferList:
            featBuffer = lyrBuffer.GetFeature(j)
            geomBuffer = featBuffer.GetGeometryRef()
            
            if geomGrid.Intersects(geomBuffer):
                cellStands.append(j)
    
        standDict[i] = cellStands

    selectedCell = max(standDict, key=lambda x:len(standDict[x]))
    return selectedCell
def osm2tif(standsfn,bbox):      
    def reprojectToWGS84(standsfn,standsWGSfn):  # creates 'standsWGS.geojson'
        if os.path.exists(standsWGSfn):
              os.remove(standsWGSfn)
        command = "ogr2ogr -f GeoJSON -t_srs EPSG:4326 %s %s" %(standsWGSfn, standsfn)
        os.system(command)
    
    def osmRoadsAPI(standsWGSfn):
        extent =  getBbox(standsWGSfn)
        bboxCoords = str(extent[0]) + ',' + str(extent[2]) + ',' + str(extent[1]) + ',' + str(extent[3])

        url = 'http://www.overpass-api.de/api/xapi?way[highway=*][bbox=%s]' % bboxCoords
        osm = requests.get(url)

        # write osm_data to osm file
        file = open(r'OSMroads.osm', 'w')
        file.write(osm.text)
        file.close()

    def osm2geojson(osmRoadsGeoJSONfn):
        roadsDs = ogr.Open('OSMroads.osm')
        inLayer = roadsDs.GetLayer(1) # layer 1 for ways

        outDriver = ogr.GetDriverByName('GeoJSON')
    
        if os.path.exists(osmRoadsGeoJSONfn):
            outDriver.DeleteDataSource(osmRoadsGeoJSONfn)

        outDataSource = outDriver.CreateDataSource(osmRoadsGeoJSONfn)
        outLayer = outDataSource.CreateLayer(osmRoadsGeoJSONfn, geom_type=ogr.wkbLineString )

        # create the input SpatialReference
        sourceSR = inLayer.GetSpatialRef()

        # create the output SpatialReference
        targetSR = osr.SpatialReference()
        targetSR.ImportFromWkt('PROJCS["Albers Equal Area",GEOGCS["grs80",DATUM["unknown",SPHEROID["Geodetic_Reference_System_1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",43],PARAMETER["standard_parallel_2",48],PARAMETER["latitude_of_center",34],PARAMETER["longitude_of_center",-120],PARAMETER["false_easting",600000],PARAMETER["false_northing",0],UNIT["Meter",1]]')

        # create transform
        coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)

        # Get the output Layer's Feature Definition
        featureDefn = outLayer.GetLayerDefn()

        # loop through the input features
        inFeature = inLayer.GetNextFeature()
        while inFeature:

            # get the input geometry
            geom = inFeature.GetGeometryRef()
            # reproject the geometry
            geom.Transform(coordTrans)
        
            # create a new feature
            outFeature = ogr.Feature(featureDefn)

            # set new geometry
            outFeature.SetGeometry(geom)
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)

            # destroy the features and get the next input feature
            outFeature.Destroy
            inFeature.Destroy
            inFeature = inLayer.GetNextFeature()

        # Close DataSources
        roadsDs.Destroy()
        outDataSource.Destroy()

    def geojson2geotiff(rasterfn,bbox,osmRoadsGeoJSONfn,osmRoadsTiffn):    
        xmin,xmax,ymin,ymax = bbox
        xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    
        source_ds = ogr.Open(osmRoadsGeoJSONfn)
        source_layer = source_ds.GetLayer()
    
        target_ds = gdal.GetDriverByName('GTiff').Create(osmRoadsTiffn, xsize, ysize, gdal.GDT_Byte)
        target_ds.SetGeoTransform((xmin, pixelWidth, 0, ymax, 0, pixelHeight))
        band = target_ds.GetRasterBand(1)
        NoData_value = 255
        band.SetNoDataValue(NoData_value)
        band.FlushCache()
        gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])   

        target_dsSRS = osr.SpatialReference()
        target_dsSRS.ImportFromWkt('PROJCS["Albers Equal Area",GEOGCS["grs80",DATUM["unknown",SPHEROID["Geodetic_Reference_System_1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",43],PARAMETER["standard_parallel_2",48],PARAMETER["latitude_of_center",34],PARAMETER["longitude_of_center",-120],PARAMETER["false_easting",600000],PARAMETER["false_northing",0],UNIT["Meter",1]]')
        target_ds.SetProjection(target_dsSRS.ExportToWkt())
        
    def main(standsfn,bbox):
        standsWGSfn = 'standsWGS.geojson'
        osmRoadsGeoJSONfn = 'OSMroads.geojson'
        osmRoadsTiffn = 'OSMroads.tif'
        reprojectToWGS84(standsfn,standsWGSfn)  # creates 'standsWGS.geojson'
        osmRoadsAPI(standsWGSfn) # creates 'roads.osm'
        osm2geojson(osmRoadsGeoJSONfn) # creates 'osmroads.geojson'
        geojson2geotiff(costSurfacefn,bbox,osmRoadsGeoJSONfn,osmRoadsTiffn) # creates 'OSMroads.tif'
        os.remove('standsWGS.geojson')
        os.remove('OSMroads.osm')
    
    if __name__ == "__main__":
        main(standsfn,bbox)
     
def getBbox(shpfn):
    ds = ogr.Open(shpfn)
    lyr = ds.GetLayer()
    bbox = lyr.GetExtent()
    return bbox

def coord2pixelOffset(rasterfn,x,y):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    xOffset = int((x - originX)/pixelWidth)
    yOffset = int((y - originY)/pixelHeight)
    return xOffset,yOffset,pixelWidth,pixelHeight

def bbox2pixelOffset(rasterfn,bbox):
    xmin,xmax,ymin,ymax = bbox
    pxmin,pymin,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmin,ymin)
    pxmax,pymax,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmax,ymax)
    xsize = abs(pxmax - pxmin)
    ysize = abs(pymax - pymin)
    xoff = pxmin
    yoff = pymin-ysize
    return xoff,yoff,xsize,ysize,pixelWidth,pixelHeight #xoff,yoff are the counts from the origion of the raster. xsize (rasterwidth),ysize(rasterheight) are the cols,rows of the raster
    
def raster2array(rasterfn,standsfn):
    
    bbox = getBbox(standsfn)
    xmin,xmax,ymin,ymax = bbox
    offset = 200.0
    xmin -= offset
    xmax += offset
    ymin -= offset
    ymax += offset
    bbox = xmin,xmax,ymin,ymax
    raster = gdal.Open(rasterfn)
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray(xoff, yoff, xsize, ysize)
    return array, bbox  
        
def array2raster(newCostSurfacefn,rasterfn,bbox,array):
    xmin,xmax,ymin,ymax = bbox
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    raster = gdal.Open(rasterfn)
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newCostSurfacefn, xsize, ysize, gdal.GDT_Byte)
    outRaster.SetGeoTransform((xmin, pixelWidth, 0, ymax, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
        
def createPath(CostSurfacefn,costSurfaceArray,selectedCell,gridfn):   
    # get index of selected cell
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    featSelectedCell = lyrGrid.GetFeature(selectedCell)
    geomSelectedCell = featSelectedCell.GetGeometryRef()
    centroidSelectedCell = geomSelectedCell.Centroid()
    selectedCellX = centroidSelectedCell.GetX()
    selectedCellY = centroidSelectedCell.GetY()
    selectedCellIndexX,selectedCellIndexY,pixelWidth,pixelHeight = coord2pixelOffset(CostSurfacefn,selectedCellX,selectedCellY)

    # get index of existing road (first point that occurs)
    roadIndexX, roadIndexY = np.where(costSurfaceArray == 0)
    # create path
    indices, weight = route_through_array(costSurfaceArray, (selectedCellIndexX,selectedCellIndexY), (roadIndexX[0],roadIndexY[0]))
    indices = np.array(indices).T
    path = np.zeros_like(costSurfaceArray)
    path[indices[0], indices[1]] = 1
    
    # merge path with existing roads
    costSurfaceArray[path == 1] = 0
    return costSurfaceArray
 
    
def main(standsfn,bufferfn,gridfn,costSurfacefn,osmRoadsfn,newCostSurfacefn):
    buffer(standsfn, bufferfn) # creates 'buffer_stands.shp'
    selectedCell = selectCell(gridfn,bufferfn) # creates string 'selectedCell'
    costSurfaceArray, bbox = raster2array(costSurfacefn,standsfn) # creates array 'costSurfaceArray' and float 'bbox'
    
    osm2tif(standsfn,bbox) # creates 'OSMroads.tif' (existing OSM roads)
    
    osmRoadsArray, bbox = raster2array(osmRoadsfn,standsfn) # creates array 'osmRoadsArray' and float 'bbox'
    costSurfaceArray[osmRoadsArray == 1.0] = 0 # updates array 'costSurfaceArray'    
    costSurfaceArray = createPath(costSurfacefn,costSurfaceArray,selectedCell,gridfn) # updates array 'costSurfaceArray'
    
    array2raster(newCostSurfacefn,costSurfacefn,bbox,costSurfaceArray) # writes 'costSurfaceArray' to 'newCostSurface.tif'
    
    
    
        
if __name__ == "__main__":
    standsfn = 'stands.shp'
    gridfn = 'grid.shp'
    bufferfn = 'buffer_stands.shp'
    costSurfacefn = 'slope.tif'
    osmRoadsfn = 'osmRoads.tif'
    newCostSurfacefn = 'newCostSurface.tif'
    main(standsfn,bufferfn,gridfn,costSurfacefn,osmRoadsfn,newCostSurfacefn)