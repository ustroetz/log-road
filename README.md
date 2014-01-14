Log Road
===========
### Overview
Log Road models logging roads based on existing OSM street data, forest stands, and a cost surface.

####Input

OSM street data | Forest Stands | Cost Surface
--- | --- | ---
OSM | Shapefile <sub><sup>(EPSG3857)</sup></sub> | GeoTIFF <sub><sup>(EPSG32610)</sup></sub>
![Alt text](/Images/InputOSM.png) | ![Alt text](/Images/InputStands.png) | ![Alt text](/Images/InputCostSurface.png)
 
####Process     
<p align="center">
  <img src="/Images/Process.gif" />
</p>

* Buffer from the stands are created based on a specified skidding distance.
* Fishnet grid covering the buffer extent is created.
* Dictionary is created. For each grid cell the intersecting buffers are recorderd: `{GridCellID:[BufferID1,BufferID2,...]}`
* OSM data are downloaded and get burned into the Cost Surface
* Cost Surface gets updated: Areas outside of the stands are more expensive, areas at the boundary of the stands are less expensive
* Buffers intersecting with OSM streets are removed from the dicitionary. 
* Loop over Dictionary:
    * Cell from dictionary with most intersecting buffers is selected.
    * New road from cell to closest road is created.
    * Buffers intersecting with new road are removed from the dicitionary. 
* New roads are written to output shapefile.

####Output
New Roads |
--- | 
Shapefile <sub><sup>(EPSG3857)</sup></sub>  |
![Alt text](/Images/Output.png) |

### Installation

Requires `python-gdal`, `skimage.graph`, `requests`, and `numpy`

To install, simply `python setup.py install` or work directly from the root directory.

### Cost Surface
Ideally the CostSurface consists for actuall $ US values in order to estimate the actual costs of construction to roads. 
But for the porpuse to identify the length and location of the new roads, a relativ coast surface is sufficient.

`test_costsurface.py` with `CostSurface.py` can be used to create a cost surface based on variables slope and river (buffered rivers). 
In `CostSurface.py` the relation between the two can be specified. We received the best test results with the following parameters:
```
    costSurfaceArray[costSurfaceArray > 50] **=2
    costSurfaceArray[riverArray == 1] **=5
```

