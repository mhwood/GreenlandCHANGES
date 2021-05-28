
from pyproj import Transformer
from pyproj import Proj, transform
import numpy as np



def run_reprojection_test(inputCRS,outputCRS):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    if inputCRS == 4326 and outputCRS == 3413:
        x = -45
        y = 70
        x_test, y_test = transformer.transform(y, x)
        if round(x_test) == 0.0 and round(y_test) == round(-2187927.649279021):
            a = 1  # everything is fine
        else:
            raise ValueError('The reproject_polygon script is not working as expected')
    elif inputCRS == 3413 and outputCRS == 4326:
        x = 474989
        y = -2282155
        y_test, x_test = transformer.transform(x, y)
        if round(x_test) == round(-33.24277887223666) and round(y_test) == round(68.71949001916302):
            a = 1  # everything is fine
        else:
            raise ValueError('The reproject_polygon script is not working as expected')
    else:
        raise ValueError('This epsg is not available in the reprojection script test')





def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1,run_test = True):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       For 3413->4326, you put in (x,y) and get out (y,x)
    #       For 326XX->3413, you put in (x,y) and get out (x,y)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,x_column], polygon_array[:,y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    if run_test:
        run_reprojection_test(inputCRS,outputCRS)

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon