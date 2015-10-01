# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:15:13 2015

@author: ta240184
"""

#virtMat_dir="/home/ta240184/Documents/GitHub_Repo/Virtual_Materials"
#import sys
#sys.path.append(virtMat_dir)
#sys.path.append(virtMat_dir+"/Utilities")
#import tifffile as tff
#sys.path.append(virtMat_dir+"/PhysicalComputations")
#import PoreNetworkExtraction as pnex


import VirtualMaterials.Simulation.PoreNetworkExtraction as pnex
import VirtualMaterials.Utilities.Utilities  as utilities
#----------------------------------------------------------


inputFileName='/home/270.12-Modeling_PEMFC_Li/theseTristan/PSI_24BA/3DSamples/PSI_sampleDrainage_2540.tif'
outputFolder='/home/270.12-Modeling_PEMFC_Li/theseTristan/PSI_24BA/ResultsPNM/extractionResults/'

outputFileName=outputFolder+"PSI_sampleDrainage_pnmExtract_h8.mat"
print(outputFileName)

extractionResult=pnex.ExtractNetwork(image=inputFileName,phases={'void':False},seedMethod='hMaxima',seedParam=8)
pnex.SaveResults( outputFileName,extractionResult)

mTypeLink=utilities.BestMemoryType(extractionResult["imageLiens"].max())
mTypePore=utilities.BestMemoryType(extractionResult["imagePores"].max())

utilities.WriteTiff(extractionResult["imageLiens"].astype(mTypeLink),outputFileName+"_imageLiens.tif")    
utilities.WriteTiff(extractionResult["imagePores"].astype(mTypePore),outputFileName+"_imagePores.tif")
    