import ogr, gdal, osr
import os
import numpy as np
from skimage.graph import route_through_array
import requests
from math import ceil

def createProjectBbox(standsfn,offsetBbox):
    bbox = getBbox(standsfn)
    xmin,xmax,ymin,ymax = bbox
    xmin -= offsetBbox
    xmax += offsetBbox
    ymin -= offsetBbox
    ymax += offsetBbox
    bbox = xmin,xmax,ymin,ymax
    return bbox
    
def createGrid(gridfn,bbox,gridHeight,gridWidth):

    xmin,xmax,ymin,ymax = bbox
    
    
    # get rows
    rows = ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(gridfn):
        os.remove(gridfn)
    outDataSource = outDriver.CreateDataSource(gridfn)
    outLayer = outDataSource.CreateLayer(gridfn,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Close DataSources
    outDataSource.Destroy()
    
    # update bbox based on new grid
    bbox = getBbox(gridfn)
    return bbox


def createBuffer(standsfn, bufferfn, Dist=1):
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
        
    bufferList = range(lyrStandsBuffer.GetFeatureCount())
    return bufferList

def selectCell(gridfn,bufferfn,bufferList):
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    gridList = range(lyrGrid.GetFeatureCount())
    
    buffer = ogr.Open(bufferfn)
    lyrBuffer = buffer.GetLayer()
    
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
    gridList.remove(selectedCell)
    
    return selectedCell
def osm2tif(standsfn,bbox,osmRoadsfn):
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
        
    def main(standsfn,bbox,osmRoadsfn):
        standsWGSfn = 'standsWGS.geojson'
        osmRoadsGeoJSONfn = 'OSMroads.geojson'
        reprojectToWGS84(standsfn,standsWGSfn)  # creates 'standsWGS.geojson'
        osmRoadsAPI(standsWGSfn) # creates 'roads.osm'
        osm2geojson(osmRoadsGeoJSONfn) # creates 'osmroads.geojson'
        geojson2geotiff(costSurfacefn,bbox,osmRoadsGeoJSONfn,osmRoadsfn) # creates 'OSMroads.tif'
        os.remove('standsWGS.geojson')
        os.remove('OSMroads.osm')
    
    if __name__ == "__main__":
        main(standsfn,bbox,osmRoadsfn)
     
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
    
def raster2array(rasterfn,bbox):
    
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray(xoff, yoff, xsize, ysize)
    return array  

def pixelOffset2coord(rasterfn,xOffset,yOffset):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    coordX = originX+pixelWidth*xOffset
    coordY = originY+pixelHeight*yOffset    
    return coordX, coordY
    
def array2raster(newRasterfn,rasterfn,bbox,array):
    xmin,xmax,ymin,ymax = bbox
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    raster = gdal.Open(rasterfn)
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, xsize, ysize, gdal.GDT_Byte)
    outRaster.SetGeoTransform((xmin, pixelWidth, 0, ymax, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(9999)
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

def tif2shp(costSurfacefn,bbox,array,outSHPfn):
   
    # update array to set No_data value
    array[array != 0] = 9999
    array[array == 0] = 1
    
    # array2tif
    FinalCostSurfacefn = 'FinalCostSurface.tif'
    array2raster(FinalCostSurfacefn,costSurfacefn,bbox,array) # writes 'costSurfaceArray' to 'FinalCostSurface.tif'
    raster = gdal.Open(FinalCostSurfacefn)
    band = raster.GetRasterBand(1)
        
    #create new shp file
    newSHPfn = 'intermediate.shp'
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(newSHPfn):
        shpDriver.DeleteDataSource(newSHPfn)
    outDataSource = shpDriver.CreateDataSource(newSHPfn)
    outLayer = outDataSource.CreateLayer(newSHPfn, geom_type=ogr.wkbPolygon )
    
    # polygonize
    gdal.Polygonize(band, None, outLayer, -1) 
    
def RemoveBuffer(bufferfn,rasterfn,costSurfaceArray):

    standsBuffer = ogr.Open(bufferfn)
    lyrStandsBuffer = standsBuffer.GetLayer()
    bufferList = range(lyrStandsBuffer.GetFeatureCount())
    
    # create list of costSurfaceArray value = 0
    roadList = np.where(costSurfaceArray == 0)
    count = 0
    for item in roadList[0]:
        indexX = item 
        indexY = roadList[1][count]
    
    removeBufferList = []
    for j in bufferList:
        featureBuffer = lyrStandsBuffer.GetFeature(j)
        geomBuffer = featureBuffer.GetGeometryRef()      

        count = 0
        for indexX in roadList[0]:
            indexY = roadList[1][count]
            Xcoord, Ycoord = pixelOffset2coord(rasterfn,indexX,indexY)
            roadPoint = ogr.Geometry(ogr.wkbPoint)
            roadPoint.AddPoint(Xcoord, Ycoord)

            if roadPoint.Within(geomBuffer):
                if j not in removeBufferList:
                    removeBufferList.append(j)   
                    
            count += 1  
                    
      
    bufferList = (set(bufferList) - set(removeBufferList))
    return bufferList
    
def main(standsfn,costSurfacefn,newCostSurfacefn):
    offsetBbox = 200
    bbox = createProjectBbox(standsfn,offsetBbox) # creates bbox that extents standsfn bbox by specified offset
    
    bufferfn = 'buffer.shp'
    bufferList = createBuffer(standsfn, bufferfn) # creates 'buffer_stands.shp' and list of buffer features

    gridfn = 'grid.shp'
    gridHeight = gridWidth = 100
    bbox = createGrid(gridfn,bbox,gridHeight,gridWidth) # creates 'grid.shp' and updates bbox based on grid's extent
    
    costSurfaceArray = raster2array(costSurfacefn,bbox) # creates array 'costSurfaceArray' and float 'bbox'
    
    osmRoadsfn = 'osmRoads.tif'
    osm2tif(standsfn,bbox,osmRoadsfn) # creates 'OSMroads.tif' (existing OSM roads)    
    osmRoadsArray = raster2array(osmRoadsfn,bbox) # creates array 'osmRoadsArray' and 'bbox'
    
    costSurfaceArray[osmRoadsArray == 1.0] = 0 # updates array 'costSurfaceArray'  
    
    OSMCostSurfacefn = 'OSMCostSurface.tif'
    array2raster(OSMCostSurfacefn,costSurfacefn,bbox,costSurfaceArray) # writes costSurfaceArray to 'OSMCostSurface.tif'
    
    bufferList = RemoveBuffer(bufferfn,OSMCostSurfacefn,costSurfaceArray) # removes buffers touching OSM roads
    
    while len(bufferList) > 0:
        
        selectedCell = selectCell(gridfn,bufferfn,bufferList) # creates string 'selectedCell'
        print 'selected cell:', selectedCell
        
        costSurfaceArray = createPath(costSurfacefn,costSurfaceArray,selectedCell,gridfn) # updates array 'costSurfaceArray'
                
        bufferList = RemoveBuffer(bufferfn,OSMCostSurfacefn,costSurfaceArray) # removes buffers touching new road
        print 'buffers left: ',len(bufferList)        
        print bufferList

    outSHPfn = 'Roads.shp'
    tif2shp(costSurfacefn,bbox,costSurfaceArray,outSHPfn) # writes final roads in shapefile
    
#   os.remove('buffer_stands.shp')
#   os.remove('grid.shp')
#   os.remove('osmRoads.tif')
    
    
    
        
if __name__ == "__main__":
    standsfn = 'stands1.shp'
    costSurfacefn = 'slope1.tif'
    newCostSurfacefn = 'newCostSurface.tif'
    main(standsfn,costSurfacefn,newCostSurfacefn)