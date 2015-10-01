# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 15:18:49 2014

@author: TA240184
"""

import numpy as np
from scipy import ndimage
from skimage import morphology
import tifffile as tff
import pylab
import mahotas
from skimage import io as skimio

foldername='C:/Users/ta240184/Desktop/PSI_detection_conditions_limites/1770_pixels/'

img=tff.imread(foldername+'dry_wholeSample_medFilter4_1770.tif')

supportExterieur=tff.imread(foldername+'support_exterieur_1770.tif')
supportExterieur= (supportExterieur == 103)
img[ supportExterieur ] = 255

#ball=skimage.morphology.ball(radius)
#skimage.morphology.binary_opening(img, selem=ball)

membrane=skimio.MultiImage(foldername+"boundary_detected_1770.tif")
membrane=membrane.concatenate().astype(np.bool)
structuringElement = np.ones((3,3,3))
membraneLabel=ndimage.measurements.label(membrane, structure=structuringElement )[0]


#membraneLabel=tff.imread('C:/Users/ta240184/Desktop/PSI_detection_conditions_limites/membraneLabel.tif')
membraneSup = (membraneLabel == 30)
membraneInf = (membraneLabel == 50)

structuringElement=np.array([[[0, 0, 0],
        [0, 1, 0],
        [0, 0, 0]],

       [[0, 1, 0],
        [1, 1, 1],
        [0, 1, 0]],

       [[0, 0, 0],
        [0, 1, 0],
        [0, 0, 0]]], dtype=bool)

for i in range(2):
    membraneInf=mahotas.morph.dilate(membraneInf,structuringElement)
    membraneSup=mahotas.morph.dilate(membraneSup,structuringElement)

imgDrainage = img
imgDrainage[ membraneSup ] = 255
imgDrainage[ membraneInf ] = 0
imgDrainage[ supportExterieur ] = 255

tff.imsave('PSIimgDrainage_1770.tif',imgDrainage)
