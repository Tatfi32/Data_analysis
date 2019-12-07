import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#function that returns dDst_dt
def Dst_prediction(dst,Q,tau):
    dDst_dt = (Q - dst/tau)
    return dDst_dt

def Dst_correction(Dst,b,Pdyn,c):
    return(Dst-b*Pdyn**0.5+c)

path = r'omni.csv'
data = pd.read_csv(path)

#read table from the site with the real numper of AL for comparison purposes
day =data.iloc[:, 1]
hour = data.iloc[:, 2]
time = day*24 + hour
Bz = data.iloc[:, 3]
V = data.iloc[:, 4]
Pdyn = data.iloc[:, 5]
Dst = data.iloc[:, 6]
#average of Pdyn
for i in range(len(Pdyn)):
    if Pdyn[i]>99:
        Pdyn[i]=0
#mean=Pdyn.mean()
for i in range(len(Pdyn)):
    if (Pdyn[i])<0.000001:
        Pdyn[i]=2.125355111444693
#print (Pdyn.describe())
        
#creating VBs  and Q data arrays
Q_AK1= np.zeros(len(Bz))  
Q_AK2= np.zeros(len(Bz))  
Q_UCB =np.zeros(len(Bz))  
VBs = np.zeros(len(Bz))
VBz= np.zeros(len(Bz))
for i in range(len(VBs)):
    VBz[i]=Bz[i]*V[i]*10**(-3)
    if VBz[i] < 0:
        VBs[i]=np.abs(VBz[i]) 
    Q_AK1[i]=(-1)*2.47*VBs[i]       #for AK1
    Q_AK2[i]=(-1)*4.4*(VBs[i]-0.5)  #for AK2
    if (VBs[i])>0.5:
        Q_UCB[i]=(-1)*4.32*VBs[i]-0.5*(Pdyn[i])**(1/3)  #for UCB
    else:
        Q_UCB[i]=0
#UCB
b=15.8
c=20  #nT
err_UCB= []
UCB= np.zeros(len(VBs))   
UCB[0] = Dst[0]
UCB_correction=np.zeros(len(VBs))   
for i in range(len(VBs)-1):    
    if (VBs[i])>4:
        tau=3
    else:
        tau=7.7      
    if (i%24)!=0:
        UCB[i+1]=odeint(Dst_prediction,UCB[i],Q_UCB[i], args=(tau,))
        UCB_correction[i] = Dst_correction(UCB[i],b,Pdyn[i],c)
        err_UCB.append((Dst[i] - UCB_correction[i])**2)
    else:
        UCB[i+1]=odeint(Dst_prediction,Dst[i],Q_UCB[i], args=(tau,))
        UCB_correction[i] = Dst_correction(UCB[i],b,Pdyn[i],c)

#plot results
plt.figure(figsize=(15,10))
#plt.plot(time,UCB,'r-',label='Theory without correction')
plt.plot(time,UCB_correction,'b-',label='UCB Theory with correction')
plt.plot(time,Dst,'g-',label='Real DST')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()

plt.plot(err_UCB)
plt.ylabel('Mean square error during the time UCB')
plt.show()
print ('Error UCB ', (sum(err_UCB)/len(err_UCB)))
print ('Mean square error without the peak ', (sum(err_UCB[3000:6000]) / len(err_UCB[3000:6000])))
#AK2
b=7.26
c=11  #nT
err_AK2= []
AK2= np.zeros(len(VBs))   
AK2[0] = Dst[0]
AK2_correction=np.zeros(len(VBs))   
for i in range(len(VBs)-1):    
    tau=2.4*(2.718261826**(9.74/(4.69+VBs[i])))  
    if (i%24)!=0:
        AK2[i+1]=odeint(Dst_prediction,AK2[i],Q_AK2[i], args=(tau,))
        AK2_correction[i] = Dst_correction(AK2[i],b,Pdyn[i],c)
        err_AK2.append((Dst[i] - AK2_correction[i])**2)
    else:
        AK2[i+1]=odeint(Dst_prediction,Dst[i],Q_AK2[i], args=(tau,))
        AK2_correction[i] = Dst_correction(AK2[i],b,Pdyn[i],c)
      
#plot results
plt.figure(figsize=(15,10))
#plt.plot(time,AK1,'r-',label='Theory without correction')
plt.plot(time,AK2_correction,'b-',label='AK2 Theory with correction')
plt.plot(time,Dst,'g-',label='Real DST')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()

plt.plot(err_AK2)
plt.ylabel('Mean square error during the time AK2')
plt.show()
print ('Mean square error during the time ', (sum(err_AK2) / len(err_AK2)))
print ('Mean square error without the peak ', (sum(err_AK2[3000:6000]) / len(err_AK2[3000:6000])))
#AK1
tau=17.0 #hour
b=8.74
c=11.05  #nT
err_AK1= []
AK1= np.zeros(len(VBs))   
AK1[0] = Dst[0]
AK1_correction=np.zeros(len(VBs))   
for i in range(len(VBs)-1):
    if (i%24)!=0:
        AK1[i+1]=odeint(Dst_prediction,AK1[i],Q_AK1[i], args=(tau,))
        AK1_correction[i] = Dst_correction(AK1[i],b,Pdyn[i],c)
        err_AK1.append((Dst[i] - AK1_correction[i])**2)
    else:
        AK1[i+1]=odeint(Dst_prediction,Dst[i],Q_AK1[i], args=(tau,))
        AK1_correction[i] = Dst_correction(AK1[i],b,Pdyn[i],c)
      
#plot results
plt.figure(figsize=(15,10))
#plt.plot(time,AK1,'r-',label='AK1 Theory without correction')
plt.plot(time,AK1_correction,'b-',label='AK1 Theory with correction')
plt.plot(time,Dst,'g-',label='Real DST')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()

plt.plot(err_AK1)
plt.ylabel('Mean square error during the time AK1')
plt.show()
print ('Mean square error during the time ', (sum(err_AK1) / len(err_AK1)))
print ('Mean square error without the peak ', (sum(err_AK1[3000:6000]) / len(err_AK1[3000:6000])))

"""
#plot results
plt.figure(figsize=(15,10))
plt.plot(time,AK1_correction,'r-',label='AK1 ')
plt.plot(time,AK2_correction,'b-',label='AK2 Theory with correction')
plt.plot(time,Dst,'g-',label='Real DST')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()
"""