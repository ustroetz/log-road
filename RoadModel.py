import ogr, gdal
import os
import numpy as np
from skimage.graph import route_through_array

def Buffer(standsfn, bufferfn, Dist=250):
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

def getBbox(shpfn):
    ds = ogr.Open(shpfn)
    lyr = ds.GetLayer()
    bbox = lyr.GetExtent()
    return bbox

def array2raster(array,rasterfn,newrasterfn,bbox):
 
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    
    originX, originY, rasterWidth,rasterHeight = bbox2pixelOffset(rasterfn,bbox)
    xmin,xmax,ymin,ymax = bbox

    
    xOffset = int((xmin - originX)/pixelWidth)
    yOffset = int((ymin - originY)/pixelHeight)
    print xOffset, yOffset
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelHeight)+1

    print xcount, ycount
    dataset = gdal.GetDriverByName('GTiff').Create(newrasterfn, xcount, ycount, gdal.GDT_Byte)
    dataset.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
        ))
            
    rasterSRS.ImportFromWkt(raster.GetProjectionRef())
    dataset.SetProjection(rasterSRS.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache() 


def coord2pixelOffset(rasterfn,x,y):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    xOffset = int((x - originX)/pixelWidth)
    yOffset = int((y - originY)/pixelHeight)
    return xOffset,yOffset

def bbox2pixelOffset(rasterfn,bbox):
    xmin,xmax,ymin,ymax = bbox
    pxmin,pymin = coord2pixelOffset(rasterfn,xmin,ymin)
    pxmax,pymax = coord2pixelOffset(rasterfn,xmax,ymax)
    rasterWidth = abs(pxmax - pxmin)
    rasterHeight = abs(pymax - pymin)
    originX = pxmin
    originY = pymin-rasterHeight
    return originX, originY, rasterWidth,rasterHeight 
    
def raster2array(rasterfn,shpfn):
    bbox = getBbox(shpfn)
    xmin,xmax,ymin,ymax = bbox
    bboxOffset = 10.0
    xmin = xmin-bboxOffset
    xmax = xmax+bboxOffset
    ymin = ymin-bboxOffset
    ymax = ymax+bboxOffset
    bbox = xmin,xmax,ymin,ymax
    raster = gdal.Open(rasterfn)
    originX, originY, rasterWidth, rasterHeight = bbox2pixelOffset(rasterfn,bbox)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray(originX, originY, rasterWidth, rasterHeight)
    
    
    newrasterfn = 'blbla.tif'
    array2raster(array,rasterfn,newrasterfn,bbox)    
    quit()
    return array
    
def createPath(costSurfacefn,costSurfaceArray,selectedCell,gridfn):   
    # get index of selected cell
    grid = ogr.Open(gridfn)
    lyrGrid = grid.GetLayer()
    featSelectedCell = lyrGrid.GetFeature(selectedCell)
    geomSelectedCell = featSelectedCell.GetGeometryRef()
    centroidSelectedCell = geomSelectedCell.Centroid()
    selectedCellX = centroidSelectedCell.GetX()
    selectedCellY = centroidSelectedCell.GetY()
    
    selectedCellIndexX,selectedCellIndexY = coord2pixelOffset(costSurfacefn,selectedCellX,selectedCellY)
    
    # get index of existing road (first point that occurs)
    roadIndex = np.nonzero(costSurfaceArray == 1)
    
    # create path
    indices, weight = route_through_array(costSurfaceArray, (selectedCellIndexX,selectedCellIndexY), tuple(roadIndex[0]))
    indices = np.array(indices).T
    path = np.zeros_like(image)
    path[indices[0], indices[1]] = 1
    
    # merge path with existing roads
    costSurfaceArray[path == 1] = 0
    
def main(standsfn,bufferfn,gridfn,costSurfacefn,osmRoadsfn):
    Buffer(standsfn, bufferfn)
    selectedCell = selectCell(gridfn,bufferfn) 
    costSurfaceArray = raster2array(costSurfacefn,standsfn)
    osmRoadsArray = raster2array(osmRoadsfn,standsfn)
    costSurfaceArray[osmRoadsArray == 1.0] = 0
    createPath(costSurfacefn,costSurfaceArray,selectedCell,gridfn)
    
        
if __name__ == "__main__":
    standsfn = 'stands.shp'
    gridfn = 'grid.shp'
    bufferfn = 'buffer_stands.shp'
    costSurfacefn = 'slopeBig.tif'
    osmRoadsfn = 'osmRoads.tif'
    main(standsfn,bufferfn,gridfn,costSurfacefn,osmRoadsfn)