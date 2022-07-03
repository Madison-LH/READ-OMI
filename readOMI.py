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
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    print(type(f[a_group_key])) 

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
    
    #print(np.array(dataset))
    dataset2= f.get('HDFEOS')
    print(dataset2)
    print(np.array(dataset2))
    additionalgrp= np.array(G1.get('ADDITIONAL'))
    print(additionalgrp)
    SWATHSgrp= np.array(G1.get('SWATHS'))
    print(SWATHSgrp)
    
    #plotting NO2 column data
    #imports
    from netCDF4 import Dataset as NetCDFFile 
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    import rioxarray

    #initalize directory for no2 data
    fn= './OMI-Aura_L2-OMNO2_2004m1001t1810-o01143_v003-2019m0814t172958-HDFEOS-SWATHS-ColumnAmountNO2-Data_Fields-ColumnAmountNO2Strat.nc'
    xr = rioxarray.open_rasterio(fn).plot()
    print(xr)




