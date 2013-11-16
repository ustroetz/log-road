import ogr, gdal, osr
import requests
import os

def reprojectToWGS84(stands):  # creates 'standsWGS.geojson'
    if os.path.exists('standsWGS.geojson'):
          os.remove('standsWGS.geojson')
    command = "ogr2ogr -f GeoJSON -t_srs EPSG:4326 %s %s" %('standsWGS.geojson', stands)
    os.system(command)

def osmRoadsAPI():
    stands = ogr.Open('standsWGS.geojson')
    standsLayer = stands.GetLayer()

    extent =  standsLayer.GetExtent()
    bboxCoords = str(extent[0]) + ',' + str(extent[2]) + ',' + str(extent[1]) + ',' + str(extent[3])

    url = 'http://www.overpass-api.de/api/xapi?way[highway=*][bbox=%s]' % bboxCoords
    osm = requests.get(url)

    # write osm_data to osm file
    file = open(r'OSMroads.osm', 'w')
    file.write(osm.text)
    file.close()

def osm2geojson():
    roadsDs = ogr.Open('OSMroads.osm')
    inLayer = roadsDs.GetLayer(1) # layer 1 for ways

    # Create the output Layer
    outDriver = ogr.GetDriverByName('GeoJSON')

    # Remove output file if it already exists
    if os.path.exists('OSMroads.geojson'):
        outDriver.DeleteDataSource('OSMroads.geojson')

    # Create the output GeoJSON
    outDataSource = outDriver.CreateDataSource('OSMroads.geojson')
    outLayer = outDataSource.CreateLayer('OSMroads.geojson', geom_type=ogr.wkbLineString )

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

def geojson2geotiff(stands):
  pixel_size = 25
  NoData_value = 255

  # Open the data source
  source_ds = ogr.Open('OSMroads.geojson')
  source_layer = source_ds.GetLayer()

  # get extent of stands shapefile
  stands = ogr.Open(stands)
  standsLayer = stands.GetLayer()
  x_min, x_max, y_min, y_max = standsLayer.GetExtent()

  # Create the destination data source
  x_res = int((x_max - x_min) / pixel_size)
  y_res = int((y_max - y_min) / pixel_size)
  target_ds = gdal.GetDriverByName('GTiff').Create('OSMroads.tif', x_res, y_res, gdal.GDT_Byte)
  target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
  band = target_ds.GetRasterBand(1)
  band.SetNoDataValue(NoData_value)

  # Rasterize
  gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])   
    
def main(stands):
    reprojectToWGS84(stands)  # creates 'standsWGS.geojson'
    osmRoadsAPI() # creates 'roads.osm'
    osm2geojson() # creates 'osmroads.geojson'
    geojson2geotiff(stands) # creates 'OSMroads.tif'

    os.remove('standsWGS.geojson')
    os.remove('OSMroads.osm')

if __name__ == "__main__":
    stands = 'stands.shp'
    main(stands)