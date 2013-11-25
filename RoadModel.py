import ogr, gdal, osr
import os
import numpy as np
from skimage.graph import route_through_array
import requests
from math import ceil, sqrt
import itertools


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
        


def array2shp(array,outSHPfn,newCostSurfacefn):
    
    # max distance between points
    raster = gdal.Open(newCostSurfacefn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))

    # array2dict
    count = 0
    roadList = np.where(array == 0)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(newCostSurfacefn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

    # dict2wkbMultiLineString
    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    for i in itertools.combinations(pointDict.values(), 2):
        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(i[0][0],i[0][1])
        point2 = ogr.Geometry(ogr.wkbPoint)
        point2.AddPoint(i[1][0],i[1][1])

        distance = point1.Distance(point2)

        if distance < maxDistance:
            line = ogr.Geometry(ogr.wkbLineString)
            line.AddPoint(i[0][0],i[0][1])
            line.AddPoint(i[1][0],i[1][1])
            multiline.AddGeometry(line)

    # wkbMultiLineString2shp
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(outSHPfn, geom_type=ogr.wkbMultiLineString )
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multiline)
    outLayer.CreateFeature(outFeature)
    
def bbox2pixelOffset(rasterfn,bbox):
    xmin,xmax,ymin,ymax = bbox
    pxmin,pymin,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmin,ymin)
    pxmax,pymax,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmax,ymax)
    xsize = abs(pxmax - pxmin)
    ysize = abs(pymax - pymin)
    xoff = pxmin
    yoff = pymin-ysize
    return xoff,yoff,xsize,ysize,pixelWidth,pixelHeight #xoff,yoff are the counts from the origion of the raster. xsize (rasterwidth),ysize(rasterheight) are the cols,rows of the raster

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

        
def createGrid(gridfn,bbox,gridHeight,gridWidth,offsetBbox):

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
    bbox = createProjectBbox(standsfn,offsetBbox)
    return bbox


def createBuffer(standsfn, bufferfn, skidDist):
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
        standCentroid = geomStand.Centroid()
        
        bufferDistList = []
        for i in xrange(geomStand.GetGeometryCount()):
            ring = geomStand.GetGeometryRef(i)
            for j in xrange(ring.GetPointCount()):
                standEdge = ogr.Geometry(ogr.wkbPoint)
                standEdge.AddPoint(ring.GetPoint(j)[0], ring.GetPoint(j)[1])
                dist = standCentroid.Distance(standEdge)
                bufferDistList.append(dist)
        distStandEdge = max(bufferDistList)
        bufferDist = max(skidDist - distStandEdge, 0)
        geomBufferStand = geomStand.Buffer(bufferDist)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBufferStand)
        lyrStandsBuffer.CreateFeature(outFeature)
        outFeature.Destroy()
        featureStand.Destroy()
        featureStand = lyrStands.GetNextFeature()
        
    bufferList = range(lyrStandsBuffer.GetFeatureCount())
    return bufferList
    
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
    StraightLineDict = {}
    count = 0
    roadIndexY, roadIndexX = np.where(costSurfaceArray == 0)
    while count < len(roadIndexY):
        dist =  sqrt((selectedCellIndexY-roadIndexY[count])**2+(selectedCellIndexX-roadIndexX[count])**2)
        index = (roadIndexY[count],roadIndexX[count])        
        StraightLineDict[index] = dist
        count +=1
        
    roadIndex = min(StraightLineDict, key=StraightLineDict.get)

    # create path
    indices, weight = route_through_array(costSurfaceArray, (selectedCellIndexY,selectedCellIndexX), roadIndex)
    indices = np.array(indices).T
    path = np.zeros_like(costSurfaceArray)
    path[indices[0], indices[1]] = 1
    
    # merge path with existing roads
    costSurfaceArray[path == 1] = 0
    return costSurfaceArray

def createProjectBbox(standsfn,offsetBbox):
    bbox = getBbox(standsfn)
    xmin,xmax,ymin,ymax = bbox
    xmin -= offsetBbox
    xmax += offsetBbox
    ymin -= offsetBbox
    ymax += offsetBbox
    bbox = xmin,xmax,ymin,ymax
    return bbox    
    
def getBbox(shpfn):
    ds = ogr.Open(shpfn)
    lyr = ds.GetLayer()
    bbox = lyr.GetExtent()
    return bbox

def getGridWidth(costSurfacefn,gridWidth):
    if gridWidth == None:
        raster = gdal.Open(costSurfacefn)
        geotransform = raster.GetGeoTransform()
        gridWidth = geotransform[1]
        
    return gridWidth
            
def osm2tif(bbox,costSurfacefn,osmRoadsTiffn):
    def reprojectToWGS84(bbox):  
        leftbottom = ogr.Geometry(ogr.wkbPoint)
        leftbottom.AddPoint(bbox[0], bbox[2])
        righttop = ogr.Geometry(ogr.wkbPoint)
        righttop.AddPoint(bbox[1], bbox[3])
        inSpatialRef = osr.SpatialReference()
        inSpatialRef.ImportFromWkt('PROJCS["Albers Equal Area",GEOGCS["grs80",DATUM["unknown",SPHEROID["Geodetic_Reference_System_1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",43],PARAMETER["standard_parallel_2",48],PARAMETER["latitude_of_center",34],PARAMETER["longitude_of_center",-120],PARAMETER["false_easting",600000],PARAMETER["false_northing",0],UNIT["Meter",1]]')
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromEPSG(4326)

        coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        leftbottom.Transform(coordTransform)
        righttop.Transform(coordTransform)

        bboxWGS84 = (leftbottom.GetX(), righttop.GetX(),leftbottom.GetY(), righttop.GetY())
        return bboxWGS84
    
    def osmRoadsAPI(bboxWGS84):
        bboxCoords = str(bboxWGS84[0]) + ',' + str(bboxWGS84[2]) + ',' + str(bboxWGS84[1]) + ',' + str(bboxWGS84[3])
        url = 'http://www.overpass-api.de/api/xapi?way[highway=*][bbox=%s]' % bboxCoords
        osm = requests.get(url)
        file = open(r'OSMroads.osm', 'w')
        file.write(osm.text)
        file.close()

    def osm2shp(osmRoadsSHPfn):
        roadsDs = ogr.Open('OSMroads.osm')
        inLayer = roadsDs.GetLayer(1) # layer 1 for ways

        outDriver = ogr.GetDriverByName('ESRI Shapefile')
    
        if os.path.exists(osmRoadsSHPfn):
            outDriver.DeleteDataSource(osmRoadsSHPfn)

        outDataSource = outDriver.CreateDataSource(osmRoadsSHPfn)
        outLayer = outDataSource.CreateLayer(osmRoadsSHPfn, geom_type=ogr.wkbLineString )

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
            inFeature = inLayer.GetNextFeature()

        # Close DataSources
        roadsDs = None
        outDataSource = None

    def shp2geotiff(rasterfn,bbox,osmRoadsSHPfn,osmRoadsTiffn):    
        xmin,xmax,ymin,ymax = bbox
        xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    
        source_ds = ogr.Open(osmRoadsSHPfn)
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
        
    def main(bbox,costSurfacefn,osmRoadsTiffn):
        osmRoadsSHPfn = 'OSMroads.shp'
        bboxWGS84 = reprojectToWGS84(bbox)  # reprojects bbox to WGS84
        #osmRoadsAPI(bboxWGS84) # creates 'roads.osm'
        osm2shp(osmRoadsSHPfn) # creates 'osmroads.shp'
        shp2geotiff(costSurfacefn,bbox,osmRoadsSHPfn,osmRoadsTiffn) # creates 'OSMroads.tif'
        #os.remove('OSMroads.osm')
    
    if __name__ == "__main__":
        main(bbox,costSurfacefn,osmRoadsTiffn)

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

def raster2array(rasterfn,bbox):
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray(xoff, yoff, xsize, ysize)
    return array  

  
def removeBuffer(bufferfn,bufferList,rasterfn,costSurfaceArray):

    standsBuffer = ogr.Open(bufferfn)
    lyrStandsBuffer = standsBuffer.GetLayer()
    
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
        for indexY in roadList[0]:
            indexX = roadList[1][count]
            Xcoord, Ycoord = pixelOffset2coord(rasterfn,indexX,indexY)
            roadPoint = ogr.Geometry(ogr.wkbPoint)
            roadPoint.AddPoint(Xcoord, Ycoord)

            if roadPoint.Within(geomBuffer):
                if j not in removeBufferList:
                    removeBufferList.append(j)   
                    
            count += 1  
                    
      
    bufferList = (set(bufferList) - set(removeBufferList))
    return bufferList
    
def selectCell(gridfn,bufferfn,bufferList):
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    gridList = range(lyrGrid.GetFeatureCount())
    
    buffer = ogr.Open(bufferfn)
    lyrBuffer = buffer.GetLayer()
    
    gridDict = {}
    
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
    
        gridDict[i] = cellStands

    selectedCell = max(gridDict, key=lambda x:len(gridDict[x]))
    removeBufferList = gridDict[selectedCell]
    bufferList = (set(bufferList) - set(removeBufferList))
    
    return selectedCell, bufferList



     


def main(standsfn,costSurfacefn,newRoadsfn,gridWidth=None,skidDist=100):
    offsetBbox = 20
    bufferfn = 'buffer.shp'
    gridfn = 'grid.shp'
    osmRoadsTiffn = 'OSMRoads.tif'
    newCostSurfacefn = 'newCostSurface.tif' 
    gridHeight = gridWidth = getGridWidth(costSurfacefn,gridWidth) # creates gridWidth&gridHeight from raster width if not specified      

    bbox = createProjectBbox(standsfn,offsetBbox) # creates bbox that extents standsfn bbox by specified offset
    
    bufferList = createBuffer(standsfn, bufferfn, skidDist) # creates 'buffer_stands.shp' and list of buffer features

    bbox = createGrid(gridfn,bbox,gridHeight,gridWidth,offsetBbox) # creates 'grid.shp' and updates bbox based on grid's extent
    
    costSurfaceArray = raster2array(costSurfacefn,bbox) # creates array 'costSurfaceArray' and float 'bbox'
    
    osm2tif(bbox,costSurfacefn,osmRoadsTiffn) # creates 'OSMroads.tif' (existing OSM roads)    
    osmRoadsArray = raster2array(osmRoadsTiffn,bbox) # creates array 'osmRoadsArray' and 'bbox'
    costSurfaceArray[osmRoadsArray == 1.0] = 0 # updates array 'costSurfaceArray'  
    array2raster(newCostSurfacefn,costSurfacefn,bbox,costSurfaceArray) # writes costSurfaceArray to 'OSMCostSurface.tif'
    
    print "Original buffers: ", bufferList
    bufferList = removeBuffer(bufferfn,bufferList,newCostSurfacefn,costSurfaceArray) # removes buffers touching OSM roads
    print "Buffers without existing roads: ", bufferList

    while bufferList:
        
        selectedCell, bufferList = selectCell(gridfn,bufferfn,bufferList) # creates string 'selectedCell'
        print 'Selected cell: ', selectedCell
        
        costSurfaceArray = createPath(newCostSurfacefn,costSurfaceArray,selectedCell,gridfn) # updates array 'costSurfaceArray'

        bufferList = removeBuffer(bufferfn,bufferList,newCostSurfacefn,costSurfaceArray) # removes buffers touching new road
        print "Buffers without new roads: ", bufferList
        
        costSurfaceArray[costSurfaceArray != 0] = 9999 # updates array 'costSurfaceArray'

    array2shp(costSurfaceArray,newRoadsfn,newCostSurfacefn) # writes final roads in shapefile
        

    
        
if __name__ == "__main__":
    standsfn = 'stands2.shp'
    costSurfacefn = '/Volumes/GIS/Basedata/PNW/terrain/slope'
    newRoadsfn = 'newRoads.shp'
    
    main(standsfn,costSurfacefn,newRoadsfn)