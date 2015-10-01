# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:15:30 2015

@author: ta240184
"""

virtMat_dir="/home/ta240184/Documents/GitHub_Repo/Virtual_Materials"
import numpy as np
#import sys
#sys.path.append(virtMat_dir)
#sys.path.append(virtMat_dir+"/Utilities")
#import tifffile as tff
#sys.path.append(virtMat_dir+"/PhysicalComputations")
#import SimpleITK as sitk

import VirtualMaterials.Utilities.tifffile as tff
import VirtualMaterials.Simulation.FullMorphology as fm
import math
#----------------------------------------------------------
folder = "/home/270.12-Modeling_PEMFC_Li/theseTristan/PSI_24BA/"

inputFileName = folder+"3DSamples/PSI_sampleDrainage_2540.tif"
outputFileName = folder+"ResultsFullMorphology/24BA_2540_FullMorpho_theta117.tif"

inputImage=tff.imread(inputFileName).astype(np.uint8)


#Parameters for Full Morphology
voxelLength = 2.2e-6
gamma = 72e-3
   
   
#memoryType=np.float16
#itkimage = sitk.GetImageFromArray(np.logical_not(inputImage==0).astype(np.uint8))
#itkdistanceMap = sitk.DanielssonDistanceMap( itkimage )
#distanceMap=sitk.GetArrayFromImage(itkdistanceMap).astype(memoryType)   
   
#maxBallRadius = int(distanceMap.max())
#   
#pressureList =   [2*gamma/(ballRadius*voxelLength)  for ballRadius in range(1,maxBallRadius,4) ]
#pressureCode =   [100+ballRadius for ballRadius in range(1,maxBallRadius,4)]


pressureList=-np.asarray([1400,2200,2800,3900,5300])/math.cos(117*math.pi/180.0)    
pressureCode = [114,122,128,139,153] 


#Perform full morphology    
outputImage = fm.FullMorphology(inputImage,inletFace=1,voxelLength=voxelLength,
                                pressureList=pressureList,
                                pressureCode=pressureCode,gamma=gamma)

#Save output    
tff.imsave(outputFileName,outputImage.astype(np.uint8))




