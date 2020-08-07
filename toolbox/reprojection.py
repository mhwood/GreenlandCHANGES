
from pyproj import Proj, transform
import numpy as np

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    inProj = Proj(init='epsg:'+str(inputCRS))
    outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transform(inProj, outProj, polygon_array[:,x_column], polygon_array[:,y_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon