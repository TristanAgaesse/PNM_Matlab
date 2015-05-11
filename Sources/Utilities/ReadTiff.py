#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:11:32 2015

@author: ta240184
"""
import sys
import hdf5storage
import os
import numpy as np

scriptDirectory = sys.argv[1]
outputDirectory = sys.argv[2]
outputName = sys.argv[3]

sys.path.append(scriptDirectory)
import tifffile as tff


matFile = outputDirectory+os.sep+'temp_'+outputName+'.mat'

data= hdf5storage.loadmat(matFile)
image=data['image']

outputFileName = outputDirectory+os.sep+outputName+'.tiff'

tff.imsave(outputFileName,image.astype(np.uint16))
