# unit test suite for RoadModel.py
# run 'py.test unit_test.py'

# missing test functions: functions creating files and arrays

from RoadModel import *

def test_createProjectBbox():
    
    # Input
    standsfn = 'stands.shp'
    offsetBbox = 200

    # Expected Output
    bbox = (274057.7229842426, 278445.30065590295, 1064147.15267288, 1067090.7083014385)

    assert(createProjectBbox(standsfn,offsetBbox)) == bbox
    

def test_createGrid():
    
    # Input
    gridfn = 'grid.shp'
    bbox = (274057.7229842426, 278445.30065590295, 1064147.15267288, 1067090.7083014385)

    gridHeight = gridWidth = 100

    # Expected Output
    outbbox = (274057.7229842426, 278457.7229842426, 1064090.7083014385, 1067090.7083014385)

    
    assert(createGrid(gridfn,bbox,gridHeight,gridWidth)) == outbbox

    
def test_selectCell():
    
    # Input
    gridfn = 'grid.shp'
    bufferfn = 'buffer.shp'
    
    # Expected Output
    selectedCell = 623
    
    assert(selectCell(gridfn,bufferfn)) == selectedCell

def test_getBbox():
    
    # Input
    shpfn = 'stands.shp'
    
    # Expected Output
    bbox = (274257.7229842426, 278245.30065590295, 1064347.15267288, 1066890.7083014385)
    
    assert(getBbox(shpfn)) == bbox
    
def test_coord2pixelOffset():
    
    # Input 
    rasterfn = 'slope.tif'
    x = 274057.722984 
    y = 1064090.7083
    
    # Expected Output
    xOffset = 82  
    yOffset = 498
    pixelWidth = 9.082676986761792
    pixelHeight = -9.082676986761792
    
    assert(coord2pixelOffset(rasterfn,x,y)) == (xOffset, yOffset, pixelWidth, pixelHeight)
    
def test_bbox2pixelOffset():
    
    # Input
    rasterfn = 'slope.tif'
    bbox = (274057.7229842426, 278457.7229842426, 1064090.7083014385, 1067090.7083014385)
    
    # Expected Output
    xoff = 82
    yoff = 168
    xsize = 484
    ysize = 330
    pixelWidth = 9.082676986761792
    pixelHeight = -9.082676986761792
    
    assert(bbox2pixelOffset(rasterfn,bbox)) == (xoff,yoff,xsize,ysize,pixelWidth,pixelHeight)
    

    
