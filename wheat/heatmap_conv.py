# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 11:34:31 2019

@author: X992029
"""

import cv2   
from matplotlib import pyplot as plt
import os
import sys
import glob
from PIL import Image
import numpy as np



def plot_image(img):
    """ Plot images
    Parameters
    ----------
    img : image (RGB)
        image
           
    """
    plt.figure()
    plt.axis("on")
    plt.imshow(img)
    plt.show()


def scale_number(unscaled, to_min, to_max, from_min, from_max):
    """ Scale number
    Parameters
    ----------
    unscaled: unscaled number
    to_min:
    to_max:
    from_min:
    from_max:
    """
    return (to_max-to_min)*(unscaled-from_min)/(from_max-from_min)+to_min


def extract_data_image(im, h, w):
    """ extract data from histogram digital image
    Keyword arguments:
    img -- image
    h -- number of rows in histogram
    w -- number of cols in histogram
    """
    im_h, im_w = im.shape
    w_num, h_num = im_w/w, im_h/h # 
    wi = np.arange(0, im_w, w_num)
    hi = np.arange(0, im_h, h_num)
    heatmap_array = np.zeros((h, w))
    for x in range(0, len(wi)): # len(wi)-1
        for y in range(0, len(hi)): # len(hi)-1
             
            crop_img = im[int(hi[y]+2): int(hi[y]+h_num), int(wi[x]+2): int(wi[x]+w_num)]
            hist = cv2.calcHist([crop_img],[0],None,[256],[0,256])
            n = np.where(hist == hist.max())[0][0]
            #n = round(np.mean(crop_img.flatten()))
            heatmap_array[y,x] = n
    return(heatmap_array)


def write_expressiondata(heatmap_array, t_min, t_max, family):
    """ export expression data to file
     Keyword arguments:
    heatmap_array - data extracted from histogram
    t_min -- min intensity in histogram
    t_max -- max intensity in histogram
    family -- gene family
    """
   

    minin = np.amin(heatmap_array)
    maxim = np.amax(heatmap_array)
    f = ''
    for x in range(0, heatmap_array.shape[0]): 
        f = f + "%s,%s," % (geneid[(x+1)][:-1], family)
        for y in range(0, heatmap_array.shape[1]): 
            heatmap_array[x,y]  = scale_number(heatmap_array[x, y], t_min, \
                         t_max, maxim, minin)
            f = f + "%5.2f, " % heatmap_array[x,y]
        f = f + "\n"
    return(f)



w = (10) # x bins
h = (91, 18, 38) # y bins
t_max = (8.73, 4.62, 5.22) # Max value in heatmap
t_min = 0 # Min value in heatmap
family_name = ('chl', 'wrky', 'myb')
f1 = open('data/sag_gene_data.csv', 'w')
gene_data = ''
stages = ('GS0','GS10', 'GS20', 'GS30',	'GS40',	'GS50',	'GS60',	'GS70',	\
              'GS80','GS90')
gene_data = "%s,%s,"% ('geneid', 'name')

for i in range(0, len(stages)):
    gene_data = gene_data + "%s," %stages[i]
gene_data = gene_data + "\n"
    
for i in range(0, len(family_name)):
    fname = 'data/' + '%s'%family_name[i] + '_genes.png'
    im = cv2.imread(fname,1) # Read image
    fname = 'data/' + '%s'%family_name[i] + '_geneid.csv'
    filename = open(fname, "r") # Read geneid
    geneid = filename.readlines() # Read geneid
    im = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)# Change to gray
    heatmap_array = extract_data_image(im, h[i], w)
    fname = 'data/' + '%s'%family_name[i] + '_gene_data.csv'
    gene_data = gene_data + write_expressiondata(heatmap_array, t_min, t_max[i], family_name[i]) 
    plt.imshow(heatmap_array, cmap='hot', interpolation='nearest')
    plt.show()

f1.write(gene_data)
f1.close()