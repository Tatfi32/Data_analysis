import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def A_calc(i, m, AL):
    return (bAtab[m][0, i] + 
    bAtab[m][1, i]*np.log10(abs(AL)) +
    bAtab[m][2, i]*(np.log10(abs(AL))**2) + 
    bAtab[m][3, i]*(np.log10(abs(AL))**3))

def Theta_calc(m, t, AL):
    return (A_calc(0, m, AL) +
            A_calc(1, m, AL)*np.cos(np.pi/180*15*(t + alp(1, m, AL))) +
            A_calc(2, m, AL)*np.cos(np.pi/180*15*(2*t + alp(2, m, AL))) +
            A_calc(3, m, AL)*np.cos(np.pi/180*15*(3*t + alp(3, m, AL))))

def alp(i, m, AL):
    return (batab[m][0, i-1] + batab[m][1, i-1]*np.log10(abs(AL)) + 
    batab[m][2, i-1]*(np.log10(abs(AL))**2) + 
    batab[m][3, i-1]*(np.log10(abs(AL))**3))

#RANGE
AL_arr = [-10, -100, -200, -300, -400, -500, -600, -700, -800]
MLT_arr =[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
#A data array
bAtab = np.array([
        #first case
    [[-0.07, -10.06, -4.44, -3.77],
    [24.54, 19.83, 7.47, 7.90],
    [-12.53, -9.33, -3.01, -4.73],
    [2.15, 1.24, 0.25, 0.91]],
         #second case
    [[1.61, -9.59, -12.07, -6.56],
    [23.21, 17.78, 17.49, 11.44],
    [-10.97, -7.20, -7.96, -6.73],
    [2.03, 0.96, 1.15, 1.31]],
         #therd case
    [[3.44, -2.41, -0.74, -2.12],
    [29.77, 7.89, 3.94, 3.24],
    [-16.38, -4.32, -3.09, -1.67],
    [3.35, 0.87, 0.72, 0.31]]
])
#alpha data array
batab = np.array([
        #first case
    [[-6.61, 6.37, -4.48],
    [10.17, -1.10, 10.16],
    [-5.80, 0.34, -5.87],
    [1.19, -0.38, 0.98]],
         #second case
    [[-2.22, -23.98, -20.07],
    [1.50, 42.79, 36.67],
    [-0.58, -26.96, -24.20],
    [0.08, 5.56, 5.11]],
         #therd case
    [[-1.68, 8.69, 8.61],
    [-2.48, -20.73, -5.34],
    [1.58, 13.03, -1.36],
    [-0.28, -2.14, 0.76]]
])

case0 = [] # m = 0 - Poleward boundary
case1 = [] # m = 1 - Equatorward boundary
case2 = [] # m = 3 - Diffusive aurora boundary

for i in range(len(AL_arr)):
    for j in MLT_arr:
        case0.append(90 - Theta_calc(0, j, AL_arr[i]))
    plt.figure()   
    plt.xlabel('MLT')
    plt.ylabel('Latitude')
    plt.title(' AL=: %i' %AL_arr[i])
    plt.plot(MLT_arr, case0, label = "Poleward boundary")
    case0.clear()
    
    for l in MLT_arr:
        case1.append(90 - Theta_calc(1, l, AL_arr[i]))
     
    plt.xlabel('MLT')
    plt.ylabel('Latitude')
    plt.title('AL=: %i' %AL_arr[i])
    plt.plot(MLT_arr, case1, label = "Equatorward boundary")
    case1.clear()
    
    for z in MLT_arr:
        case2.append(90 - Theta_calc(2, z, AL_arr[i]))  
     
    plt.xlabel('MLT')
    plt.ylabel('Latitude')
    plt.title('AL=: %i' %AL_arr[i])
    plt.plot(MLT_arr, case2, label = "Diffusive aurora boundary")
    plt.legend()
    plt.show()
    case2.clear()
  
    
#Second part
path = r'omni.csv'
data = pd.read_csv(path)

#read table from the site with the real numper of AL for comparison purposes
day =data.iloc[:, 1]
hour = data.iloc[:, 2]
time = day*24 + hour
By = data.iloc[:, 3]
Bz = data.iloc[:, 4]
V = data.iloc[:, 5]
Kp = data.iloc[:, 6]
AL_data = data.iloc[:, 7]

#calculate theory values
alpha = 0.0044
BzBy_sqrt = np.sqrt(By**2+Bz**2)
BzBy_sqrt[BzBy_sqrt==0] = 1
theta = np.arccos(Bz/BzBy_sqrt)
ALkp = (18 - 12.3*Kp + 27.2*Kp**2 - 2*Kp**3)/300

#plot thetheory values and the real ones
plt.figure(figsize=(15,10))
plt.plot(time, ALkp, label="theory AL data", color= 'r')
plt.plot(time, AL_data, label="real AL data",color= 'b')
plt.xlabel('Time in hours')
plt.ylabel('AL values')
plt.legend()
plt.title('AL from the Kp reconstruction:')
plt.show()
print ('Mean square error =', np.average(np.sqrt(np.abs(ALkp**2 - AL_data**2))))


# Solar wind data reconstruction:
#calculate theoretical E
E = V*np.sqrt(By**2/2 + Bz**2)*(np.sin(theta/2)**4) + alpha * V**2 * np.sin(theta/2)**0.5
E_new = np.zeros(len(E)-2)
for i in range(len(E_new)):
    E_new[i] = (E[i] + 8*E[i+1] + 5*E[i+2])/14

lg = np.log10(E_new)
c = np.polyfit(lg,AL_data[2:], 4)
AL_new = c[0]*lg**4 + c[1]*lg**3 + c[2]*lg**2 + c[3]*lg + c[4]

#plot thetheory values and the real ones
plt.figure(figsize=(15,10))
plt.plot(time[2:], AL_data[2:], label="real AL data")
plt.plot(time[2:], AL_new, label="theory AL data")
plt.xlabel('Time in hours')
plt.title(' AL from solar wind data reconstruction:')
plt.ylabel('AL values')
plt.legend()
plt.show()
print ('Mean square error =', np.average(np.sqrt(np.abs(AL_new**2 - AL_data[2:]**2))))



