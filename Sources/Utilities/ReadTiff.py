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
inputDirectory = sys.argv[2]
inputName = sys.argv[3]

sys.path.append(scriptDirectory)
import tifffile as tff

#Read tiff file
inputFileName = inputDirectory+os.sep+inputName

image = tff.imread(inputFileName)
image = image.astype(np.uint16)

#Write temporary mat file
matFile = inputDirectory+os.sep+'temp_'+inputName+'.mat'
print(matFile)
hdf5storage.savemat(matFile,{'temp_image':image})


