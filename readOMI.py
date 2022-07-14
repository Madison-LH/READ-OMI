#imports
from grp import struct_group
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from regex import P
import xarray as xr
import h5py
import netCDF4
#import pdb


 #single file
#dataDIR = './OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5'
#'/Users/madisonhonore/Documents/VS code/OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958-HDFEOS-SWATHS-ColumnAmountNO2-Data_Fields-AmfTropClear.nc'
#replace the above with the location of your OMI data
#DS = np.fromfile('./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5', dtype=float)

#f1= h5py.File('./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5', 'r')
#df= pd.read_hdf('./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5')


#datasetNames= [n for n in f1.keys()]
#for n in datasetNames:
    #print(n)

#pdb.set_trace()

#hf= h5py.File('./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5', 'r')
#print(hf.keys())
#n1= hf.get('HDFEOS INFORMATION')
#n1
#n1= np.array(n1)
#n1.shape
#print(np.array(n1))

#above was test code to get the names of groups in the .he5 file

#initalize file directory
filename = "./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5"

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    #print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    #print(type(f[a_group_key])) 

    # If a_group_key is a group name, 
    # this gets the object names in the group and returns as a list
    data = list(f[a_group_key])

    # If a_group_key is a dataset name, 
    # this gets the dataset values and returns as a list
    data = list(f[a_group_key])
    # preferred methods to get dataset values:
    ds_obj = f[a_group_key]      # returns as a h5py dataset object
    #ds_arr = f[a_group_key][()]  # returns as a numpy array
    dataset= f.get('HDFEOS INFORMATION')
    print(dataset)
    G1= f.get('HDFEOS INFORMATION')
    G1_items= list(G1.items())
    print('items in group 1: ', G1_items)
    archmetadatagrp= np.array(G1.get('ArchivedMetadata.0'))
    print(archmetadatagrp)
    coremetadatagrp= np.array(G1.get('CoreMetadata.0'))
    print(coremetadatagrp)
    structmetadatagrp= np.array(G1.get('StructMetadata.0'))
    print(structmetadatagrp)
    
    print(np.array(dataset))
    dataset2= f.get('HDFEOS')
    print(dataset2)
    print(np.array(dataset2))
    additionalgrp= np.array(G1.get('ADDITIONAL'))
    print(additionalgrp)
    SWATHSgrp= np.array(G1.get('SWATHS'))
    print(SWATHSgrp)
    

    #plotting NO2 column data
    #imports
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import h5py

FILE_NAME = "./OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958.he5"
path = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'
DATAFIELD_NAME = path + 'CloudFraction'
with h5py.File(FILE_NAME, mode='r') as f:
    dset = f[DATAFIELD_NAME]
    data =dset[:].astype(np.float64)

    # Retrieve any attributes that may be needed later.
    # String attributes actually come in as the bytes type and should
    # be decoded to UTF-8 (python3).
    scale = f[DATAFIELD_NAME].attrs['ScaleFactor']
    offset = f[DATAFIELD_NAME].attrs['Offset']
    missing_value = f[DATAFIELD_NAME].attrs['MissingValue']
    fill_value = f[DATAFIELD_NAME].attrs['_FillValue']
    title = f[DATAFIELD_NAME].attrs['Title'].decode()
    units = f[DATAFIELD_NAME].attrs['Units'].decode()

    # Retrieve the geolocation data.
    path = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'
    latitude = f[path + 'Latitude'][:]
    longitude = f[path + 'Longitude'][:]

    data[data == missing_value] = np.nan
    data[data == fill_value] = np.nan
    data = scale * (data - offset)
    datam = np.ma.masked_where(np.isnan(data), data)

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat = 90,
                llcrnrlon=-180, urcrnrlon = 180)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.scatter(longitude, latitude, c=datam, s=1, cmap=plt.cm.jet,
             edgecolors=None, linewidth=0)    
    cb = m.colorbar()
    cb.set_label(units)


    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title), fontsize=8)
    fig = plt.gcf()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)





