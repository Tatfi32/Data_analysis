##########################################################################
#Libs
import scipy as sp
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda
import pandas as pd
import cv2
import os
from scipy import ndimage
from astropy.io import fits
from scipy.stats import pearsonr

###############################################################
#GOES soft X-ray flux in the 1–8 Å
goes = pd.read_csv('goes.txt', sep="  ", header=None)
goes.columns = ["Time (in sec)", "X-ray flux (in Wm^2)"]
Xray = goes["X-ray flux (in Wm^2)"]
timeGOES = goes["Time (in sec)"]/60
plt.figure()
plt.plot(timeGOES,Xray)
plt.yscale('log')

plt.xlabel('Time (min)')
plt.ylabel('X-ray flux (Wm−2)')
plt.title('GOES soft X-ray flux in the 1–8 Å')

print ('Xray Max ',np.max([Xray]))
print ('Xray Min ',np.min([Xray]))

##############################################################################
#AIA images
images = []
PixNum = []
KmPerPix = []
path = 'AIA_24_cal'
files= os.listdir(path)
#print(len(files))
for image_file in files:
    image_file =path + '/'+ image_file
    hdulist = fits.open(image_file, memmap=True)   
    image_data = fits.getdata(image_file, ext=0)

    
    #print(image_data.shape)
##############################################################################
#Print image

    #plt.figure()
    #plt.imshow(image_data)
    #plt.colorbar()
    
##############################################################################
#From fts file read info

    hdr = hdulist[0].header
    #print(repr(hdulist[0].header))
    #hdulist.info()
    RSUN_REF = hdulist[0].header['RSUN_REF'] #RSUN_REF - radius of Sun in m
    #print(RSUN_REF)
    DATE_OBS = hdulist[0].header['DATE-OBS']
    R_SUN = hdulist[0].header['R_SUN'] 
    #print(R_SUN)                            #R_SUN - Radius of the Sun’s image in pixels.
    hdulist.close() 
    MetrPerPix = (RSUN_REF/R_SUN)**2
    #print('km^2 per pix ',MetrPerPix/1000000) 
    KmPerPix.append(MetrPerPix/1000000)
##############################################################################
# print it
    
    syn_map = sunpy.map.Map(image_data,hdulist[0].header)
    #plt.show()
    #im = syn_map.plot()
    
#############################################################################
# Flare Square (median filter)
    image = image_data[600:850,1100:1400]
    image[image<4000] = 0
    image[image>=3000] = 255
    image = cv2.medianBlur(image,3) 
    PixNum.append(np.count_nonzero(image))
    Square =  np.array(PixNum)*MetrPerPix
    images.append(image)
    
#Total Square over time
TotalSquare = np.sum([Square])
print ('Total square (km^2) ', TotalSquare/1e+6)
#############################################################################
# Flare Square plot
plt.figure()
time = np.linspace(0, 59, len(Square))
plt.plot(time,Square)
plt.title('Flare square over time')
plt.ylabel('Flare area (m2)')
plt.xlabel('Time (min)')

# Diff X ray plot
difX = np.diff(Xray)
plt.figure()
time = np.linspace(0, 59, len(difX))
plt.plot(time,difX)
plt.yscale('log')
plt.title('Time derivative of the GOES over time')
plt.ylabel('Time derivative (W*m−2*min-1)')
plt.xlabel('Time (min)')



# X-ray flux  VS Flare Area interpolation

plt.figure()
plt.scatter(Square [5:45] ,difX[5:45])
polinomFirst = np.poly1d(np.polyfit(Square [5:20], difX[5:20], 1))
plt.plot(Square [5:20],polinomFirst(Square [5:20]),'g')
plt.ylim(0,0.00005)
plt.xlabel('Flare Area (m^2)')
plt.ylabel('X-ray flux derivative')
plt.title(' X-ray flux derivative VS Flare Area interpolation (time from 5 to 20 min)')

#############################################################################
#Calculates a Pearson correlation coefficient and the p-value for testing non-correlation
pearsonFirst = pearsonr(difX[5:20] , polinomFirst(Square[5:20]))
print  ('pearson[5:20] ' , pearsonFirst)




plt.figure()
plt.scatter(Square [5:45] ,difX[5:45])
polinomLast  = np.poly1d(np.polyfit(Square [30:45], difX[30:45], 1))
plt.plot(Square [30:45],polinomFirst(Square [30:45]),'r')
plt.ylim(0,0.00005)
plt.xlabel('Flare Area (m^2)')
plt.ylabel('X-ray flux derivative')
plt.title(' X-ray flux derivative  VS Flare Area interpolation (time from 30 to 45 min)')

#############################################################################
#Calculates a Pearson correlation coefficient and the p-value for testing non-correlation
pearsonLast = pearsonr(difX[30:45] , polinomLast(Square[30:45]))
print  ('pearson[25:45] ' , pearsonLast)


#Continuum
Con = fits.open('HMI_24_cal/HMI_continuum.fits')[0].data
Con = Con[600:850,1100:1400]
#Con  = np.flipud(Con)
Con = np.round(Con/239).astype('uint8')
#Magnetogram
Mag = fits.open('HMI_24_cal/HMI_magnetogram.fits')[0].data
Mag  = Mag [600:850,1100:1400]
#Mag = np.flipud(Mag)
Mag  = np.round((Mag +2815.1694)/19).astype('uint8')
#MaxFlare
MaxSquare = np.max([Square])
max_Square_index=np.where(Square == MaxSquare)
#print (max_Square_index)
MaxFlare=images[12]


#MinFlare
MinSquare = np.min([Square])
min_Square_index=np.where(Square == MinSquare)
#print (min_Square_index)
MinFlare=images[20] 


#contours on MaxFlare image
_, contours, _= cv2.findContours(MaxFlare.astype('uint8'), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) 
NewMaxFlare = cv2.drawContours(MaxFlare.astype('uint8'), contours, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on  AIA 1700 (at time of flare maximum) ')
plt.imshow(NewMaxFlare)

#contours on MinFlare image
_, contours2, _= cv2.findContours(MinFlare.astype('uint8'), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) 
NewMinFlare = cv2.drawContours(MinFlare.astype('uint8'), contours2, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on  AIA 1700 (at time of flare minimum) ')
plt.imshow(NewMinFlare)


area = []
for rows in range(MaxFlare.shape[0]):
    for columns in range(MaxFlare.shape[1]):
        for n in contours:
            if (cv2.pointPolygonTest(n, (rows,columns),measureDist=False)) == 1:
                area.append([rows,columns]) 

#contours on  Magnetogram (max)
NewMag = cv2.drawContours(Mag.astype('uint8'), contours, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on  Magnetogram (at time of flare maximum) ')
plt.imshow(NewMag)
#contours on Continuum
NewCon = cv2.drawContours(Con.astype('uint8'), contours, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on Continuum (at time of flare maximum) ')
plt.imshow(NewCon)


#contours on  Magnetogram (min)
NewMag2 = cv2.drawContours(Mag.astype('uint8'), contours2, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on  Magnetogram (at time of flare min) ')
plt.imshow(NewMag2)

#contours on Continuum (max)
NewCon = cv2.drawContours(Con.astype('uint8'), contours, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on Continuum (at time of flare maximum) ')
plt.imshow(NewCon)




#contours on Continuum (min)
NewCon2 = cv2.drawContours(Con.astype('uint8'), contours2, -1, (0, 0, 255),1)
plt.figure()
plt.title(' contours on Continuum (at time of flare min) ')
plt.imshow(NewCon2)

# magnetic flux underlying the flare ribbons, separately for the positive and negative polarity
Positive = 0
Negative = 0
EachKmPerPix = np.mean(KmPerPix)
Magnet = fits.open('HMI_24_cal/HMI_magnetogram.fits')[0].data
Magnet  = Magnet [600:850,1100:1400]
Magnet = np.flipud(Magnet)
for j in  area:
        if Magnet[j[0],j[1]]>0:
            Positive += Magnet[j[0],j[1]]
        elif Magnet[j[0],j[1]]<0:   
            Negative +=Magnet[j[0],j[1]]
print(Positive*10000000000*EachKmPerPix ,"Mx")
print(Negative*100000000000*EachKmPerPix ,"Mx")
