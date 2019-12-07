from PIL import Image, ImageDraw
import numpy as np
import cv2
import matplotlib.pyplot as plt
from PIL.TiffTags import TAGS
import pandas as pd
#import exifread


def ring(A,B,C):
    a = np.linalg.norm(C - B)
    b = np.linalg.norm(C - A)
    c = np.linalg.norm(B - A)
    s = (a + b + c) / 2
    R = a*b*c / 4 / np.sqrt(s * (s - a) * (s - b) * (s - c))
    b1 = a*a * (b*b + c*c - a*a)
    b2 = b*b * (a*a + c*c - b*b)
    b3 = c*c * (a*a + b*b - c*c)
    P = np.column_stack((A, B, C)).dot(np.hstack((b1, b2, b3)))
    P /= b1 + b2 + b3
    Radius = 71492 #km
    KmPerPix= Radius/R
    return (KmPerPix)



a_im = Image.open('heic0009a.tif')
#a_im.show(title='heic0009a.tif')
img=a_im.convert(mode='P')
a_array = np.array(img)
#contours of the disk
a_array = cv2.medianBlur(a_array,77)
a_data_array = cv2.Canny(a_array, threshold1=0, threshold2=255)
#plt.figure()
#plt.imshow(a_array)
#plt.figure()
#plt.imshow(a_data_array)
#Km per pix calculation
A = np.array([0, 590, 0])
B = np.array([1078, 798, 0])
C = np.array([2047,660 , 0])
KmPerPix1=ring(A,B,C)
#Aurora ellipse area
x, y =  a_im.size
eX, eY = 810,230 #Size of Bounding Box for ellipse
bbox =  (x/2-350 - eX/2, y/2-150 - eY/2, x/2 + eX/2, y/2 + eY/2)
draw = ImageDraw.Draw(a_im)
draw.ellipse(bbox, outline='red',width=2)
#a_im.show(title='heic0009a.tif')
del draw
Square = 3.14*810*230*KmPerPix1


b_im = Image.open('heic1613a.tif')
#b_im.show(title='heic1613a.tif')
b_array = np.array(b_im)
#b_im.show(title='heic0009a.tif')
img=b_im.convert(mode='P')
#contours of the disk
b_array = cv2.medianBlur(b_array,77)
b_data_array = cv2.Canny(b_array, threshold1=0, threshold2=255)
#plt.figure()
#plt.imshow(b_array)
#plt.figure()
#plt.imshow(b_data_array)
#Km per pix calculation
A2 = np.array([393, 86, 0])
B2 = np.array([394, 907, 0])
C2 = np.array([69, 495, 0])
KmPerPix2=ring(A2,B2,C2)
#Aurora ellipse area
x, y =  b_im.size
eX, eY = 250,180 #Size of Bounding Box for ellipse
draw = ImageDraw.Draw(b_im)
bbox =  (x/2 - eX/2, y-835 - eY/2, x/2 + eX/2, y-935 + eY/2)
#draw.ellipse((200,200,500,500),outline='red',width=2)
draw.ellipse(bbox, outline='red',width=2)
#b_im.show(title='heic0009a.tif')
del draw
Square2 = 3.14*250*180*KmPerPix2



c_im = Image.open('heic1721.tif')
#c_im.show(title='heic1721.tif')
c_array = np.array(c_im)
img=c_im.convert(mode='P')
c_data_array = np.array(img)
#contours of the disk
c_array = cv2.medianBlur(c_data_array, 45)
#c_data_array = cv2.Canny(c_array, threshold1=0, threshold2=255)
_, contours, _= cv2.findContours(c_data_array, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
cont=[]
area=[]
for contour in contours:
    M = cv2.moments(contour)
    if M['m00']>28 and M['m00']<313:
        cont.append(contour) 
        area.append(cv2.contourArea(contour))
#plot the contours      
c_array = cv2.drawContours(c_array, cont, -1, (0, 0, 255),1)
#plt.figure()
#plt.imshow(c_array)
#plt.figure()
#plt.imshow(c_data_array)
#Km per pix calculation
A3 = np.array([5, 200, 0])
B3 = np.array([488,100, 0])
C3 = np.array([240, 40, 0])
KmPerPix3=ring(A3,B3,C3)
#Aurora ellipse area
x, y =  c_im.size
eX, eY = 170,60 #Size of Bounding Box for ellipse
draw = ImageDraw.Draw(c_im)
bbox =  (x/2-10 - eX/2, y/2-70 - eY/2, x/2 + eX/2, y/2-70 + eY/2)
#draw.ellipse((200,200,500,500),outline='red',width=2)
draw.ellipse(bbox, outline='red',width=2)
#c_im.show(title='heic0009a.tif')
del draw
Square3 = 3.14*166*55*KmPerPix3



d_im = Image.open('heic011998.tif')
#d_im.show(title='heic011998.tif')
d_array = np.array(d_im)
img=d_im.convert(mode='P')
#contours of the disk
d_array = cv2.medianBlur(d_array,77)
d_data_array = cv2.Canny(d_array, threshold1=0, threshold2=255)
#plt.figure()
#plt.imshow(d_array)
#plt.figure()
#plt.imshow(d_data_array)
#Km per pix calculation
A4 = np.array([0, 129, 0])
B4 = np.array([325, 240, 0])
C4 = np.array([505, 398, 0])
KmPerPix4=ring(A4,B4,C4)
#Aurora ellipse area
x, y =  d_im.size
eX, eY = 240,100 #Size of Bounding Box for ellipse
draw = ImageDraw.Draw(d_im)
bbox =  (x/2-44 - eX/2, y/2-55 - eY/2, x/2-44 + eX/2, y/2-55 + eY/2)

draw.ellipse(bbox, outline='red',width=2)
#d_im.show(title='heic0009a.tif')
#del draw
Square4 = 3.14*240*100*KmPerPix4



path = r'Earth_Jupiter_OMNI_2000D182-2001D165_h.csv'
data = pd.read_csv(path)


data_347 = data[data['doy'] == 347].index
data_351 = data[data['doy'] == 351].index
our_data = data[data_347[0]:data_351[len(data_351)-1]]

our_data['br^2'] = our_data['br'].apply(lambda x: x**2)
our_data['bt^2'] = our_data['bt'].apply(lambda x: x**2)
our_data['bp^2'] = our_data['bp'].apply(lambda x: x**2)
our_data['Bm'] = our_data['br^2']+our_data['bt^2']+our_data['bp^2']
our_data['Bm'] = our_data['Bm'].apply(lambda x: x**0.5)




plt.figure()
plt.title ('Bm over time')
plt.ylabel('Bm [nT]')
plt.xlabel('Time (days of the year)') 
plt.plot(our_data['t'], our_data['Bm'])

our_data['Bm_min']= np.min(our_data['Bm'])


our_data['vr^2'] = our_data['vr'].apply(lambda x: x**2)
our_data['vt^2'] = our_data['vt'].apply(lambda x: x**2)
our_data['vp^2'] = our_data['vp'].apply(lambda x: x**2)
our_data['vm'] = our_data['vr^2']+our_data['vt^2']+our_data['vp^2']
our_data['vm'] = our_data['vm'].apply(lambda x: x**0.5)

plt.figure()
plt.title ('vm over time')
plt.ylabel('vm [km/s]')
plt.xlabel('time (days of the year)') 
plt.plot(our_data['t'], our_data['vm'],c = "black")




plt.figure()
plt.title ('Temperature over time')
plt.ylabel('Temperature')
plt.xlabel('Time (days of the year)') 
plt.plot(our_data['t'], our_data['T'],  c = "r")



plt.figure()
plt.title ('Density over time')
plt.ylabel('Denisity ')
plt.xlabel('Time (days of the year)') 
plt.plot(our_data['t'], our_data['rho'], c = "orange")




#contours
"""
_, contours, _= cv2.findContours(a_array, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
cont=[]
area=[]
for contour in contours:
    M = cv2.moments(contour)
    if M['m00']>15 and M['m00']<165330.9:
        cont.append(contour) 
        area.append(cv2.contourArea(contour))
#plot the contours      
a_array = cv2.drawContours(a_array, cont, -1, (0, 0, 255),1)
"""
#Image data 
"""
f = open('heic0009a.tif', 'rb')
tags = exifread.process_file(f)
for tag in tags.keys():
    if tag not in ('JPEGThumbnail', 'TIFFThumbnail', 'Filename', 'EXIF MakerNote'):
        #print ("Key: %s, value %s" % (tag, tags[tag]))    
"""
"""
with Image.open('heic0009a.tif') as img:
    meta_dict = {TAGS[key] : img.tag[key] for key in img.tag.keys()}
#print(meta_dict)
"""
