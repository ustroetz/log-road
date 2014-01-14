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
    target_dsSRS.ImportFromWkt('PROJCS["Albers Equal Area",GEOGCS["grs80",DATUM["unknown",SPHEROID["Geodetic_Reference_System_1980",6378137,298.257222101],TOWGS84[0,0,0,0,0,0,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",43],PARAMETER["standard_parallel_2",48],PARAMETER["latitude_of_center",34],PARAMETER["longitude_of_center",-120],PARAMETER["false_easting",600000],PARAMETER["false_northing",0],UNIT["Meter",1]]')
    target_ds.SetProjection(target_dsSRS.ExportToWkt())


def array2raster(rasterfn,costSurfaceArray):
    newRasterfn = 'testdata/CostSurface2.tif'     
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