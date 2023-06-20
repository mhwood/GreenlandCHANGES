#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Collection IDs
#
# MEaSUREs Annual Antarctic Ice Velocity Maps V001
# Short Name: MEaSUREs Antarctic Annual Velocity
# https://nsidc.org/data/nsidc-0720/versions/1
# https://search.earthdata.nasa.gov/search?q=NSIDC-0720%20V001
# 
# MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005
# Short Name: MEaSUREs Greenland Quarterly Velocity
# https://nsidc.org/data/nsidc-0727/versions/5
# https://search.earthdata.nasa.gov/search?q=NSIDC-0727%20V005
# 
# MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5
# Short Name: MEaSUREs Greenland Monthly Velocity
# https://nsidc.org/data/nsidc-0731/versions/5
# https://search.earthdata.nasa.gov/search?q=NSIDC-0731+V005
#  
# ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2 (ATL14) -- (Note: Cloud data also avail-- See Earthdata for more info)
# Short Name: ATL14 Antarctic LIH
# https://nsidc.org/data/atl14/versions/2
# https://search.earthdata.nasa.gov/search?q=ATL14%20V002
#
# ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2 (ATL15)
# Short Name: ATL15 Antarctic LIH Change
# https://nsidc.org/data/atl15/versions/2
# https://search.earthdata.nasa.gov/search?q=ATL15%20V002
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aaivm  = 'C2245171699-NSIDC_ECS'    # MEaSUREs Annual Antarctic Ice Velocity Maps V001
gqisvm = 'C2627036252-NSIDC_ECS'    # MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005
gmisvm = 'C2627046644-NSIDC_ECS'    # MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5
atl14   = 'C2500138845-NSIDC_ECS'   # ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2 (ATL14)
atl15   = 'C2500140833-NSIDC_ECS'   # ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2 (ATL15) 

# Each collection has a long name first and a short name second -- REQUIRED
collections = {
            'MEaSUREs Annual Antarctic Ice Velocity Maps V001': aaivm, 
            'MEaSUREs Antarctic Annual Velocity': aaivm,

            'MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005': gqisvm,
            'MEaSUREs Greenland Quarterly Velocity': gqisvm,

            'MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5': gmisvm,
            'MEaSUREs Greenland Monthly Velocity': gmisvm,

            'ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2': atl14,
            'ATL14 Antarctic Elevation': atl14,
            
            'ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2': atl15,
            'ATL15 Antarctic Elevation': atl15,
            }

# Get collection ID from long or short name
def collection(collection_key):    
    if collection_key not in collections:
        raise ValueError('Collection key not recognized')
    return collections[collection_key]

# Print available collections in dictionary
def print_collections():
    print('Available collections: \n')
    for i in range(len(collections)):
        if i % 2 == 0:
            print(list(collections.keys())[i], end='\n')
            print("Short name: " + list(collections.keys())[i+1], end='\n')
            print()
            i += 1
    return

# Get long and short name from a collection ID
def get_names(collection_key):
    if list(collections.keys()).index(collection_key) % 2 == 0:
        long_name = collection_key
        short_name = list(collections.keys())[list(collections.keys()).index(collection_key)+1]
    else:
        short_name = collection_key
        long_name = list(collections.keys())[list(collections.keys()).index(collection_key)-1]
    return long_name, short_name


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO - decide if we want to keep this functionality
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def add_to_dict(collection_key_long, collection_key_short, collection_id):
    collections[collection_key_long] = collection_id
    collections[collection_key_short] = collection_id
    return

def remove_from_dict(collection_id):
    for key in list(collections.keys()):
        if collections[key] == collection_id:
            del collections[key]
    return