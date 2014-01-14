import gdal, ogr, osr, os
import numpy as np

def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    return array  
    
def shp2raster(inputfn,baseRasterfn):    
    outputfn = 'rasterized.tif'
    
    source_ds = ogr.Open(inputfn)
    source_layer = source_ds.GetLayer()
    
    raster = gdal.Open(baseRasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    target_ds = gdal.GetDriverByName('GTiff').Create(outputfn, cols, rows, 1, gdal.GDT_Byte) 
    target_ds.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    band = target_ds.GetRasterBand(1)
    NoData_value = 0
    band.SetNoDataValue(NoData_value)
    band.FlushCache()        
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])   

    target_dsSRS = osr.SpatialReference()
    target_dsSRS.ImportFromEPSG(32610)
    target_ds.SetProjection(target_dsSRS.ExportToWkt())


def array2raster(rasterfn,costSurfaceArray):
    newRasterfn = 'testdata/CostSurface.tif'     
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3] 
    pixelWidth = geotransform[1] 
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(costSurfaceArray)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
        

def main(riverfn,slopefn):
    shp2raster(riverfn,slopefn)  
    riverArray = raster2array('rasterized.tif')
    os.remove('rasterized.tif')

    costSurfaceArray = raster2array(slopefn)

    costSurfaceArray[costSurfaceArray > 50] **=2
    costSurfaceArray[riverArray == 1] **=5
    costSurfaceArray[costSurfaceArray < 0.001] = 1 # change value of 0 to 1 otherwise they mix up with 0 values of roads
    
    array2raster(slopefn,costSurfaceArray)
    
if __name__ == "__main__":
    riverfn = 'rivers.shp'
    slopefn = 'Slope.tif'
 
    main(riverfn, slopefn)