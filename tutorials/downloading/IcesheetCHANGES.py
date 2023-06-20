# this script defines a class used to hold all pertinent information to
# compile the elevation and velocity

import time
import datetime
import os
import itertools
import requests
import collection

class IcesheetCHANGES:
    def __init__(self, project_folder, data_folder, collection_key, collection_id, region_name, projection):

        #this initiates the global domain and all of its associated parameters
        self.project_folder = project_folder
        self.data_folder = data_folder
        self.collection_key = collection_key
        self.collection_id = collection_id
        self.region_name = region_name
        self.projection = projection
        self.long_name, self.short_name = collection.get_names(self.collection_key)
        self.download_path = os.path.join(self.data_folder,self.region_name,'Velocity',self.short_name,'Data')
        self.velocity_grid_x = []
        self.velocity_grid_y = []
        self.elevation_grid_y = []
        self.elevation_grid_x = []


    # Projections are defined by EPSG codes (https://epsg.io/)
    def set_projection(self, projection):
        self.projection = projection
    def get_projection(self):
        return self.projection

class AntarcticCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, collection_key, collection_id, region_name, projection = 3031):
        super().__init__(project_folder, data_folder, collection_key, collection_id, region_name, projection)
        #self.downloaded_data_path = os.path.join(self.data_folder,self.region_name,'Velocity',self.short_name,'Data')
        pass

    def test(self):
        print('testAC')

class GreenlandCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, collection_key, collection_id, region_name, projection = 3413):
        super().__init__(project_folder, data_folder, collection_key, collection_id, region_name, projection)
        #self.downloaded_data_path = os.path.join(self.data_folder,self.region_name,'Velocity',self.short_name,'Data')
        pass 

    def test(self):
        print('testGC')

# this function is from the ICESat2 data page
def cmr_filter_urls(search_results):
    """Select only the desired data files from CMR response."""
    if 'feed' not in search_results or 'entry' not in search_results['feed']:
        return []
    
    entries = [e['links']
                for e in search_results['feed']['entry']
                if 'links' in e]
    # Flatten "entries" to a simple list of links
    links = list(itertools.chain(*entries))

    urls = []
    unique_filenames = set()
    for link in links:
        if 'href' not in link:
            # Exclude links with nothing to download
            continue    # continue jumps to next iteration in the loop
        if 'inherited' in link and link['inherited'] is True:
            # Why are we excluding these links?
            continue
        if 'rel' in link and 'data#' not in link['rel']:
            # Exclude links which are not classified by CMR as "data" or "metadata"
            continue
        if 'title' in link and 'opendap' in link['title'].lower():
            # Exclude OPeNDAP links--they are responsible for many duplicates
            # This is a hack; when the metadata is updated to properly identify
            # non-datapool links, we should be able to do this in a non-hack way
            continue

        filename = link['href'].split('/')[-1]
        if filename in unique_filenames:
            # Exclude links with duplicate filenames (they would overwrite)
            continue
        unique_filenames.add(filename)

        urls.append(link['href'])
        
    return urls

## write available file names and links to csv file
def write_to_csv(AC, availible_file_names, availible_file_links):
    if not os.path.exists(os.path.join(AC.project_folder,AC.region_name,'Velocity','Metadata')):
        os.makedirs(os.path.join(AC.project_folder,AC.region_name,'Velocity','Metadata'))
    f = open(os.path.join(AC.project_folder,AC.region_name,'Velocity','Metadata', AC.collection_key + ' ' + AC.collection_id + '.csv'),'w')

    f.write('File_Name,URL')
    for ea in range(len(availible_file_names)):
        f.write('\n'+availible_file_names[ea]+','+availible_file_links[ea])
    f.close()


def cmr_api_url(AC):
    # Build and call the CMR API URL
    cmr_query_url = 'https://cmr.earthdata.nasa.gov/search/granules.json?echo_collection_id=' + AC.collection_id + '&page_size=2000'
    response = requests.get(cmr_query_url)

    # print error code based on response
    if response.status_code != 200:
        print('ERROR: {}'.format(response.status_code))
    search_page = response.json()
    # If JSON contains an error message, print the message at the key, 'error'
    if 'errors' in search_page:
        print(search_page['errors'])
    else: 
        urls = cmr_filter_urls(search_page)
        print("Successfully obtained {} URLs.".format(len(urls)))
        # Get only .nc files and also store the file names in a seperate list
        availible_file_links = []
        availible_file_names = []
        for url in urls:
            if '.nc'in url and not 'xml' in url:
                availible_file_links.append(url)
                url_parts = url.split('/')
                file_name = url_parts[-1]
                availible_file_names.append(file_name)
        write_to_csv(AC, availible_file_names, availible_file_links)

    return availible_file_names, availible_file_links



# check for existing files, get a list of files not on disk
def obtain_download_list(AC, file_names, file_links):
    
    output_folder=os.path.join(AC.data_folder,'Velocity','MEaSUREs','Annual','Data')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    download_list_names=[]
    download_list_links=[]
    for fil in file_names:
        download_file = True
        if fil not in os.listdir(os.path.join(output_folder)):
            for existing_fil in os.listdir(os.path.join(output_folder)):
                if fil == existing_fil:  # this allows for older versions to be kept 
                    download_file = False
        else:
            download_file=False
        if download_file:
            download_list_names.append(fil)
            download_list_links.append(file_links[file_names.index(fil)])
    #n_existing_files = len(file_names)-len(download_list)
    return(download_list_names, download_list_links) #,n_existing_files)


def download_nc_direct(AC, file_names, file_links):
    not_downloaded_names = []
    not_downloaded_links = []
    print('Downloading ' + str(len(file_names)) + ' file(s)...')
    for i in range(len(file_names)):
        # progress report
        print('Downloading ' + str(i + 1) + '/' + str(len(file_names)) + ': ' + file_names[i], end='\r')
        # get request (download file)
        r = requests.get(file_links[i], allow_redirects=True)
        if r.status_code != 200:    # 200 is the standard response for successful HTTP requests
            print('ERROR: ' + str(r.status_code) + '\n')
            not_downloaded_names.append(file_names[i])
            not_downloaded_links.append(file_links[i])
        # write content to file
        open(os.path.join(AC.data_folder,'Velocity','MEaSUREs','Annual','Data',file_names[i]), 'wb').write(r.content)
    return(not_downloaded_names, not_downloaded_links)

def run_download_nc_direct(AC, dl_list_names, dl_list_links):
    # Download the files
    not_downloaded_names, not_downloaded_links = download_nc_direct(AC, dl_list_names, dl_list_links)

    # Attempt to download any files that were not downloaded the first time
    if len(not_downloaded_names) > 0:
        print('The following files were not downloaded:')
        for i in range(len(not_downloaded_names)):
            print(not_downloaded_names[i] + ': ' + not_downloaded_links[i])
            
        print('Attempting to download these files once more...')
        not_downloaded_names, not_downloaded_links = download_nc_direct(AC, not_downloaded_names, not_downloaded_links)
        if len(not_downloaded_names) > 0:
            print('Please download these files manually and place them in the following folder:')
            print(os.path.join(AC.data_folder,'Velocity','MEaSUREs','Annual','Data'))
