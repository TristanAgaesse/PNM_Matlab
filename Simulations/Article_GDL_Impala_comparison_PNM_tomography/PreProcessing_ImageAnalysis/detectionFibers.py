# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 15:30:43 2015

@author: ta240184
"""
import numpy as np
#from scipy import ndimage
#from skimage import morphology
import sys
sys.path.append('/home/ta240184/Documents/GitHub_Repo/Image_Computations')
import tifffile as tff
import time
import SimpleITK as sitk


#Script to decompose solid into fibers and binder
beginTime=time.time()
foldername='/home/ta240184/theseTristan/PSI_24BA/PreTraitement/2540_pixels/'

img=tff.imread(foldername+'PSI_boundaries_AllDetected_2540.tif').astype(np.uint8)

solid = img==255

#ball=morphology.ball(12)
#binder = ndimage.morphology.binary_opening(solid, structure=ball)
ballsize=7
itkImage = sitk.GetImageFromArray(solid.astype(np.uint8))
itkImage = sitk.BinaryErode(itkImage, int(ballsize), sitk.sitkBall, 0.0, 1.0,  False)
itkImage = sitk.BinaryDilate(itkImage, int(ballsize), sitk.sitkBall, 0.0, 1.0,  False)   
binder=sitk.GetArrayFromImage(itkImage).astype(np.bool)   


fibers=np.logical_and(solid,np.logical_not(binder))

foldername='/home/ta240184/theseTristan/PSI_24BA/3DSamples/'
sample=tff.imread(foldername+'PSI_sampleDrainage_2540.tif').astype(np.uint8)
sample[ binder ] = 225
sample[ fibers ] = 255

print('saving results')
tff.imsave(foldername+'24BA_sampleDrainage_2540_fibersIdentified_7.tif',sample.astype(np.uint8))

endTime=time.time()
print("Time spent : {} s".format(endTime-beginTime))