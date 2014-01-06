import ogr, gdal, osr
import os, sys, fnmatch
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
    outRaster = driver.Create(newRasterfn, xsize, ysize, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((xmin, pixelWidth, 0, ymax, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(9999)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
        


def array2shp(costSurfaceArray,osmRoadsArray,outSHPfn,newCostSurfacefn):
    
    
    ## New Road to geometry
    # max distance between points
    raster = gdal.Open(newCostSurfacefn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))

    # array2dict
    costSurfaceArray[osmRoadsArray == 1.0] = 1 # remove existing roads from new roads
    count = 0
    roadList = np.where(costSurfaceArray == 0)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(newCostSurfacefn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

    # dict2wkbMultiLineString
    multilineNewRoad = ogr.Geometry(ogr.wkbMultiLineString)
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
            multilineNewRoad.AddGeometry(line)
            
    # calculate length of line
    length = multilineNewRoad.Length()
    
    ## OSM road to geometry
    # max distance between points
    raster = gdal.Open(newCostSurfacefn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))

    # array2dict
    count = 0
    roadList = np.where(osmRoadsArray == 1.0)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(newCostSurfacefn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

    # dict2wkbMultiLineString
    multilineOSMRoad = ogr.Geometry(ogr.wkbMultiLineString)
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
            multilineOSMRoad.AddGeometry(line)
    
    # Transform to EPSG 3857
    inSpatialRef = osr.SpatialReference()      
    inSpatialRef.ImportFromProj4('+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(3857)
    
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    multilineNewRoad.Transform(coordTrans)
    multilineOSMRoad.Transform(coordTrans)
    
    # wkbMultiLineString2shp
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(outSHPfn, geom_type=ogr.wkbMultiLineString )
    # create a field
    idField = ogr.FieldDefn('type', ogr.OFTString)
    outLayer.CreateField(idField)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multilineNewRoad)
    outFeature.SetField('type', 'new')
    outLayer.CreateFeature(outFeature)
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multilineOSMRoad)
    outFeature.SetField('type', 'existing')
    outLayer.CreateFeature(outFeature)

    
    return length
    
def bbox2pixelOffset(rasterfn,bbox):
    xmin,xmax,ymin,ymax = bbox
    pxmin,pymin,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmin,ymin)
    pxmax,pymax,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,xmax,ymax)
    xsize = abs(pxmax - pxmin)
    ysize = abs(pymax - pymin)
    xoff = pxmin
    yoff = pymax
    return xoff,yoff,xsize,ysize,pixelWidth,pixelHeight #xoff,yoff are the counts from the origion of the raster. xsize (rasterwidth),ysize(rasterheight) are the cols,rows of the raster

def cellID2cellIndex(selectedCell,gridfn,CostSurfacefn):
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    featSelectedCell = lyrGrid.GetFeature(selectedCell)
    geomSelectedCell = featSelectedCell.GetGeometryRef()
    centroidSelectedCell = geomSelectedCell.Centroid()
    selectedCellX = centroidSelectedCell.GetX()
    selectedCellY = centroidSelectedCell.GetY()
    selectedCellIndexX,selectedCellIndexY,pixelWidth,pixelHeight = coord2pixelOffset(CostSurfacefn,selectedCellX,selectedCellY)
    return selectedCellIndexX,selectedCellIndexY
  
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

        
def checkfile(fn):
    try:
        open(fn)
    except:
        sys.exit('ERROR: File can not be found in the file system.')
def createGrid(gridfn,bbox,offsetBbox,gridHeight,gridWidth):
    xmin,xmax,ymin,ymax = bbox
    xmin -= offsetBbox
    xmax += offsetBbox
    ymin -= offsetBbox
    ymax += offsetBbox
    
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
    bbox = createProjectBbox(gridfn)
    return bbox

def createGridDict(gridfn,bufferfn):
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    gridList = range(lyrGrid.GetFeatureCount())
    
    buffer = ogr.Open(bufferfn)
    lyrBuffer = buffer.GetLayer()
    bufferList = range(lyrBuffer.GetFeatureCount())
    
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
        
    gridDict = {i: cellStands for i, cellStands in gridDict.items() if cellStands}
    return gridDict

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
        bufferDist = max(0,skidDist - distStandEdge)
        geomBufferStand = geomStand.Buffer(bufferDist)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBufferStand)
        lyrStandsBuffer.CreateFeature(outFeature)
        outFeature.Destroy()
        featureStand.Destroy()
        featureStand = lyrStands.GetNextFeature()
    

        
    
def createPath(CostSurfacefn,costSurfaceArray,selectedCell,gridfn):   
    # get index of selected cell
    selectedCellIndexX,selectedCellIndexY = cellID2cellIndex(selectedCell,gridfn,CostSurfacefn)
    
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

def createProjectBbox(standsfn):
    bbox = getBbox(standsfn)
    return bbox    
    
def getBbox(shpfn):
    ds = ogr.Open(shpfn)
    lyr = ds.GetLayer()
    bbox = lyr.GetExtent()
    return bbox

def getGridWidth(costSurfacefn):
    raster = gdal.Open(costSurfacefn)
    geotransform = raster.GetGeoTransform()
    gridWidth = geotransform[1]
        
    return gridWidth
            
def osm2tif(bbox,costSurfacefn,osmRoadsSHPfn):
    def reprojectToWGS84(bbox,offsetBbox):  
        xmin,xmax,ymin,ymax = bbox
        xmin -= offsetBbox
        xmax += offsetBbox
        ymin -= offsetBbox
        ymax += offsetBbox
        bbox = xmin,xmax,ymin,ymax
        leftbottom = ogr.Geometry(ogr.wkbPoint)
        leftbottom.AddPoint(bbox[0], bbox[2])
        righttop = ogr.Geometry(ogr.wkbPoint)
        righttop.AddPoint(bbox[1], bbox[3])
        inSpatialRef = osr.SpatialReference()
        inSpatialRef.ImportFromProj4('+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
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
        targetSR.ImportFromProj4('+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

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
        
        featureCount = outLayer.GetFeatureCount()
            
        # Close DataSources
        roadsDs = None
        outDataSource = None
        
        return featureCount
        
    def main(bbox,costSurfacefn,osmRoadsSHPfn):
        offsetBbox = 2
        bboxWGS84 = reprojectToWGS84(bbox,offsetBbox)  # reprojects bbox to WGS84
        osmRoadsAPI(bboxWGS84) # creates 'roads.osm'
        featureCount = osm2shp(osmRoadsSHPfn) # creates 'osmroads.shp'
        while featureCount == 0:
            offsetBbox **= 2
            bboxWGS84 = reprojectToWGS84(bbox,offsetBbox)
            osmRoadsAPI(bboxWGS84) # creates 'roads.osm'
            featureCount = osm2shp(osmRoadsSHPfn) # creates 'osmroads.shp'
        
        offsetBbox += 100
        return offsetBbox

    offsetBbox = main(bbox,costSurfacefn,osmRoadsSHPfn)
    return offsetBbox

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



def printList(myDict):
    myList = []
    for item in myDict:
        for jtem in myDict[item]:
            if jtem not in myList:
                myList.append(jtem)
    print myList
def poly2line(input_poly,output_line):

    source_ds = ogr.Open(input_poly)
    source_layer = source_ds.GetLayer()

    # polygon2geometryCollection
    geomcol =  ogr.Geometry(ogr.wkbGeometryCollection)
    for feat in source_layer:
        geom = feat.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        geomcol.AddGeometry(ring)
        

    # geometryCollection2shp
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(output_line):
    	shpDriver.DeleteDataSource(output_line)
    outDataSource = shpDriver.CreateDataSource(output_line)
    outLayer = outDataSource.CreateLayer(output_line, geom_type=ogr.wkbMultiLineString)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(geomcol)
    outLayer.CreateFeature(outFeature)
def raster2array(rasterfn,bbox):
    xoff, yoff, xsize, ysize, pixelWidth, pixelHeight = bbox2pixelOffset(rasterfn,bbox)
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray(xoff, yoff, xsize, ysize)
    return array  

  
def reprojectFrom3857(inputfn,outputfn):
    ds = ogr.Open(inputfn)
    inLayer = ds.GetLayer()
    
    inSpatialRef = osr.SpatialReference()    
    inSpatialRef.ImportFromEPSG(3857)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4('+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
    
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outputfn):
        shpDriver.DeleteDataSource(outputfn)
    outDataSet = shpDriver.CreateDataSource(outputfn)
    outLayer = outDataSet.CreateLayer(outputfn, geom_type=ogr.wkbPolygon)
    outLayerDefn = outLayer.GetLayerDefn()

    inFeature = inLayer.GetNextFeature()
    while inFeature:
        geom = inFeature.GetGeometryRef()
        geom.Transform(coordTrans)
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(geom)
        outLayer.CreateFeature(outFeature)
        inFeature = inLayer.GetNextFeature()

def removeBuffer(bufferfn,gridDict,rasterfn,costSurfaceArray):

    standsBuffer = ogr.Open(bufferfn)
    lyrStandsBuffer = standsBuffer.GetLayer()
    
    # create bufferList
    bufferList = []
    for cellID in gridDict:
        for bufferID in gridDict[cellID]:
            if bufferID not in bufferList:
                bufferList.append(bufferID)
    
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
                    
    gridDict = {x:[z for z in y if z not in removeBufferList] for x,y in gridDict.items()}
    gridDict = {i: cellStands for i, cellStands in gridDict.items() if cellStands}        
    return gridDict




def selectCell(gridfn,bufferfn,gridDict,costSurfaceArray,rasterfn):

    selectedCell = max(gridDict, key=lambda x:len(gridDict[x]))
    # in case stands don't intersect
    if len(gridDict[selectedCell]) == 1:
        # find closest cell of selected buffer to road

        # find closest road to centroid of buffer
        standsBuffer = ogr.Open(bufferfn)
        lyrStandsBuffer = standsBuffer.GetLayer()
        featureBuffer = lyrStandsBuffer.GetFeature(gridDict[selectedCell][0])
        geomBuffer = featureBuffer.GetGeometryRef()
        centroidBuffer = geomBuffer.Centroid()
        bufferX = centroidBuffer.GetX()
        bufferY = centroidBuffer.GetY()
        bufferIndexX,bufferIndexY,pixelWidth,pixelHeight = coord2pixelOffset(rasterfn,bufferX,bufferY)
        
        roadIndexY, roadIndexX = np.where(costSurfaceArray == 0)
        count = 0
        StraightLineDict = {}
        while count < len(roadIndexY):
            dist =  sqrt((bufferIndexY-roadIndexY[count])**2+(bufferIndexX-roadIndexX[count])**2)
            index = (roadIndexX[count],roadIndexY[count])        
            StraightLineDict[index] = dist
            count +=1
        selectedRoadIndex = min(StraightLineDict, key=StraightLineDict.get)
        selectedRoadCoordX, selectedRoadCoordY = pixelOffset2coord(rasterfn,selectedRoadIndex[0],selectedRoadIndex[1])
        selectedRoadCoord = ogr.Geometry(ogr.wkbPoint)
        selectedRoadCoord.AddPoint(selectedRoadCoordX, selectedRoadCoordY)
        
        # find closest cell to selected road 
        gridDs = ogr.Open(gridfn)
        lyrGrid = gridDs.GetLayer()
        subGridDict = {}
        count = 0
        for cell in gridDict:
            if (gridDict[cell][0] == gridDict[selectedCell][0]):
                featureCell = lyrGrid.GetFeature(cell)
                geomCell = featureCell.GetGeometryRef()
                centroidCell = geomCell.Centroid()
                cellX = centroidCell.GetX()
                cellY = centroidCell.GetY()
                cellCoords = ogr.Geometry(ogr.wkbPoint)
                cellCoords.AddPoint(cellX, cellY)
                dist = cellCoords.Distance(selectedRoadCoord)
                subGridDict[cell] = dist
        selectedCell = min(subGridDict, key=subGridDict.get)

    removeBufferList = gridDict[selectedCell]        
    gridDict = {x:[z for z in y if z not in removeBufferList] for x,y in gridDict.items()}
    gridDict = {i: cellStands for i, cellStands in gridDict.items() if cellStands}
    return selectedCell, gridDict
def shp2raster(inputSHPfn,rasterfn,outputRasterfn):    
    
    source_ds = ogr.Open(inputSHPfn)
    source_layer = source_ds.GetLayer()
    
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    target_ds = gdal.GetDriverByName('GTiff').Create(outputRasterfn, cols, rows, 1, gdal.GDT_Float32) 
    target_ds.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    band = target_ds.GetRasterBand(1)
    NoData_value = 0
    band.SetNoDataValue(NoData_value)
    band.FlushCache()        
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])   

    target_dsSRS = osr.SpatialReference()
    target_dsSRS.ImportFromWkt('PROJCS["Albers Equal Area",GEOGCS["grs80",DATUM["unknown",SPHEROID["Geodetic_Reference_System_1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",43],PARAMETER["standard_parallel_2",48],PARAMETER["latitude_of_center",34],PARAMETER["longitude_of_center",-120],PARAMETER["false_easting",600000],PARAMETER["false_northing",0],UNIT["Meter",1]]')
    target_ds.SetProjection(target_dsSRS.ExportToWkt())

def shp2array(inputSHPfn,rasterfn):  
    outputRasterfn = 'intermediate.tif'  
    shp2raster(inputSHPfn,rasterfn,outputRasterfn)
    
    raster = gdal.Open(outputRasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    os.remove(outputRasterfn)
    return array

def main(standsfn,costSurfacefn,newRoadsfn,skidDist=0):
    bufferfn = 'buffer.shp'
    gridfn = 'grid.shp'
    osmRoadsSHPfn = 'OSMroads.shp'
    newCostSurfacefn = 'newCostSurface.tif' 
    standsLinefn = 'standsLine.shp'
    reprostandsfn = 'standsReprojected.shp'
    checkfile(standsfn)
    checkfile(costSurfacefn)
    
    reprojectFrom3857(standsfn,reprostandsfn)
    
    gridHeight = gridWidth = getGridWidth(costSurfacefn) # creates gridWidth&gridHeight from raster width if not specified      

    createBuffer(reprostandsfn, bufferfn, skidDist) # creates 'buffer_stands.shp' and list of buffer features
    print 'Buffer created'
    
    bbox = createProjectBbox(bufferfn) # creates bbox based on bufferfn's bbox
    
    offsetBbox = osm2tif(bbox,costSurfacefn,osmRoadsSHPfn) # creates 'OSMroads.tif' (existing OSM roads)    
    
    bbox = createGrid(gridfn,bbox,offsetBbox,gridHeight,gridWidth) # creates 'grid.shp' and updates bbox based on grid's extent
    gridDict = createGridDict(gridfn,bufferfn) # creates dict (key=cellID; value: List of buffers intersecting cell)
    #gridDict = eval((open('gridDict.txt', 'r')).read())    
    #gridDictFile = open('gridDict.txt', 'w')
    #gridDictFile.write(str(gridDict))    
    print 'Grid created'
    
    # update array (based on OSM streets, areas outside property, stand boundaries)
    costSurfaceArray = raster2array(costSurfacefn,bbox) # creates array 'costSurfaceArray' 
    array2raster(newCostSurfacefn,costSurfacefn,bbox,costSurfaceArray) # writes costSurfaceArray to 'newCostSurface.tif'
  
    osmRoadsArray = shp2array(osmRoadsSHPfn,newCostSurfacefn) # creates array 'osmRoadsArray'
    costSurfaceArray[osmRoadsArray == 1.0] = 0 # updates array 'costSurfaceArray'
    
    standsarray = shp2array(reprostandsfn,newCostSurfacefn) # creates array from stands
    costSurfaceArray[standsarray == 0] **= 2     # updates array, area outside of stands
  
    poly2line(reprostandsfn,standsLinefn) # creates stands line shapefile
    standLinesarray = shp2array(standsLinefn,newCostSurfacefn) # creates array from stands line
    costSurfaceArray[standLinesarray == 1] /= 2    # updates array, stands boundaries    
    
    
    gridDict = removeBuffer(bufferfn,gridDict,newCostSurfacefn,costSurfaceArray) # removes buffers touching OSM roads
    
    while gridDict:
        print "remaining buffer:"
        printList(gridDict) 
            
        selectedCell, gridDict = selectCell(gridfn,bufferfn,gridDict,costSurfaceArray,newCostSurfacefn) # creates 'selectedCell'
        print selectedCell

        costSurfaceArray = createPath(newCostSurfacefn,costSurfaceArray,selectedCell,gridfn) # updates array 'costSurfaceArray'

        gridDict = removeBuffer(bufferfn,gridDict,newCostSurfacefn,costSurfaceArray) # removes buffers touching new road from gridDict
        
        
    length = array2shp(costSurfaceArray,osmRoadsArray,newRoadsfn,newCostSurfacefn) # writes final roads in shapefile and returns length of new roads in meters
        
    print 'Length new roads (miles): ', length*0.000621371 
            
    for filename in os.listdir('.'):
        for pattern in ['buffer*','grid*','newCostSurface*','standsLine*','standsReprojected*']:
            if fnmatch.fnmatch(filename, pattern):
                os.remove(filename)
    
        
if __name__ == "__main__":
    standsfn = 'testdata/test_stand.shp'
    costSurfacefn = 'testdata/CostSurface.tif'
    newRoadsfn = 'newRoad.shp'
    
    main(standsfn,costSurfacefn,newRoadsfn)