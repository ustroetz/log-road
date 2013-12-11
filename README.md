Log Road
===========
### Overview
Log Road models logging roads based on existing OSM street data, forest stands, and a cost surface.

####Input

OSM street data | Forest Stands | Cost Surface
--- | --- | ---
OSM | Shapefile | GeoTIFF
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
Shapefile |
![Alt text](/Images/Output.png) |

### Installation

Requires `python-gdal`, `skimage.graph`, `requests`, and `numpy`

To install, simply `python setup.py install` or work directly from the root directory.
