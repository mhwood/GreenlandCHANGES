#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Collection IDs
#
# MEaSUREs Annual Antarctic Ice Velocity Maps V001
# Short Name: Antarctic Annual Velocity
# https://nsidc.org/data/nsidc-0720/versions/1
# https://search.earthdata.nasa.gov/search?q=NSIDC-0720%20V001
# 
# MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005
# Short Name: Greenland Quarterly Velocity
# https://nsidc.org/data/nsidc-0727/versions/5
# https://search.earthdata.nasa.gov/search?q=NSIDC-0727%20V005
# 
# MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5
# Short Name: Greenland Monthly Velocity
# https://nsidc.org/data/nsidc-0731/versions/1
# https://search.earthdata.nasa.gov/search?q=NSIDC-0731%20V005
#  
# ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2 (ATL14) -- (Cloud data also avail-- See Earthdata for more info)
# Short Name: ATL14
# https://nsidc.org/data/atl14/versions/2
# https://search.earthdata.nasa.gov/search?q=ATL14%20V002
#
# ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2 (ATL14) -- (Cloud data also avail-- See Earthdata for more info)
# Short Name: ATL14ch
# https://nsidc.org/data/atl14/versions/2
# https://search.earthdata.nasa.gov/search?q=ATL14%20V002
#
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aaivm  = 'C2245171699-NSIDC_ECS' # MEaSUREs Annual Antarctic Ice Velocity Maps V001
gqisvm = 'C2627036252-NSIDC_ECS' # MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005
gmisvm = 'C2627046644-NSIDC_ECS' # MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5
atl14   = 'C2500138845-NSIDC_ECS' # ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2 (ATL14)
atl14ch = 'C2500140833-NSIDC_ECS' # ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2 (ATL14) 

# each collection has a long name and a short name -- REQUIRED
collections = {
            'MEaSUREs Annual Antarctic Ice Velocity Maps V001': aaivm, 
            'MEaSUREs Antarctic Annual Velocity': aaivm,

            'MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat V005': gqisvm,
            'MEaSUREs Greenland Quarterly Velocity': gqisvm,

            'MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat, Version 5': gmisvm,
            'MEaSUREs Greenland Monthly Velocity': gmisvm,

            'ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height, Version 2 (ATL14)': atl14,
            'ATL14 Antarctic LIH': atl14,
            
            'ATLAS/ICESat-2 L3B Gridded Antarctic and Arctic Land Ice Height Change, Version 2 (ATL14)': atl14ch,
            'ATL14 Antarctic LIH Change': atl14ch
            }


def collection(collection_key):    
    if collection_key not in collections:
        raise ValueError('Collection key not recognized')
    # if the collection key is indexed in an even position
    return collections[collection_key]

def print_collections():
    print('Available collections: \n')
    for i in range(len(collections)):
        if i % 2 == 0:
            print(list(collections.keys())[i], end='\n')
            print("Short name: " + list(collections.keys())[i+1], end='\n')
            print()
            i += 1
    return

def get_names(collection_key):
    if list(collections.keys()).index(collection_key) % 2 == 0:
        long_name = collection_key
        short_name = list(collections.keys())[list(collections.keys()).index(collection_key)+1]
    else:
        short_name = collection_key
        long_name = list(collections.keys())[list(collections.keys()).index(collection_key)-1]
    return long_name, short_name