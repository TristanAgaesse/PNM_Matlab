# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:00:04 2015

@author: ta240184
"""


import numpy as np
from VirtualMaterials.Utilities  import tifffile as tff



foldername='/home/270.12-Modeling_PEMFC_Li/theseTristan/PSI_24BA/3DSamples/'

img=tff.imread(foldername+'PSI_FusionImages_2540.tif')

newImg=np.zeros(img.shape,dtype=np.uint8)

GDL= (img == 50)
newImg[ GDL ] = 255

holder= (img == 200)
newImg[ holder ] = 255

membraneSup= (img == 180)
newImg[ membraneSup ] = 255

membraneInf= (img == 75)
newImg[ membraneInf ] = 255

tff.imsave(foldername+'PSI_24BA_2540_GDL.tif',newImg)


