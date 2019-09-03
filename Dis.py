import os
import glob
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import xml.dom.minidom


class point:
    def __init__(self,tag,r,proj):
        self.tag = tag
        self.r = [s for s in r]
        self.proj = proj


    def dis(self,other):
        return np.sum( (np.array(self.r)-np.array(other.r))**2 )**0.5
    def __lt__(self,other):
        return NotImplemented
    def __gt__(self,other):
        return NotImplemented
    def __le__(self,other):
        return NotImplemented
    def __ge__(self,other):
        return NotImplemented

class step_plist:
    def __init__(self,plist):
        self.plist = [s for s in plist]


def compute_msd(position, shifts):
    msds = np.zeros(len(shifts))
    stds = np.zeros(len(shifts))
    
    for i,shift in enumerate(shifts):
        diffs = position[:-shift if shift else None] - position[shift:]
        sqdlist = np.square(diffs)
        msds[i] = sqdlist.mean()
        stds[i] = sqdlist.std(ddof=1)
    return msds,stds


#---------------------------------------------------------------------------------------#
#			        Get input file
#---------------------------------------------------------------------------------------#

str = raw_input("Enter input file name: ")

str_time = str + '/SIGEPS'

str_COP = str + '/COP*.dat'
#---------------------------------------------------------------------------------------#
#				  Get time
#---------------------------------------------------------------------------------------#

tt = []    #in ns
tt.append(0.0)    #time zero

fcon = open(str_time)
line = fcon.readline()    #skip first line

line = fcon.readline()

while line : 
    line = line.strip()
    if line == "":
        line = fcon.readline()
        continue
    
    else:
        line = line.split()
        tt.append(float(line[1]))
        line = fcon.readline()

#---------------------------------------------------------------------------------------#
#				  Get initial position
#---------------------------------------------------------------------------------------#
fplist = []    #list to all plist for different steps

dom = xml.dom.minidom.parse('graph.xml')

root = dom.documentElement

itemlist = root.getElementsByTagName('object')

plist = []

for p_number in range(len(itemlist)):
    item = itemlist[p_number]
    plist.append(point(float(p_number),[float(item.getAttribute('cx')),float(item.getAttribute('cy')),float(item.getAttribute('cz'))],float(item.getAttribute('cz'))))

fplist.append(step_plist(plist))
#---------------------------------------------------------------------------------------#
#				  Get points
#---------------------------------------------------------------------------------------#

flist = glob.glob(str_COP)
flist.sort()

for ff in flist:
    fcon = open(ff)
    plist = []
    line = fcon.readline()    #skip first line
    line = fcon.readline()
    
    while line : 
        line = line.strip()
        if line == "":
            line = fcon.readline()
            continue

        else:
            line = line.split()
            plist.append(point(float(line[0]),[float(s) for s in line[1:4]],float(line[4])))
            line = fcon.readline()
    
    fplist.append(step_plist(plist))


dr = []
dr2 = []
dx = []
dy = []


for step in range(len(tt)):
    dr.append(fplist[step].plist[0].proj * 1.16 - fplist[0].plist[0].proj * 1.16)
    dx.append(fplist[step].plist[0].r[0] - fplist[0].plist[0].r[0])
    dy.append(fplist[step].plist[0].r[1] - fplist[0].plist[0].r[1])
    dr2.append((fplist[step].plist[0].proj - fplist[0].plist[0].proj)**2.0)

sum_z = 0.0
dt_mean = 0.0

for step in range(1,len(tt)):
    dt_mean += tt[step]-tt[step-1]
    sum_z += (tt[step]-tt[step-1])*dr[step]

dt_mean /= 1.0*(len(tt)-1.0)
print dt_mean

dr_cor = [dr[0]]

for step in range(1,len(tt)):
    dr_cor.append(dr[step] - 2.0*tt[step]/tt[-1]/tt[-1]*sum_z)






dshift = [10,15,20,50,60,80]


tau = [s*dt_mean for s in dshift]
msd_t,msd_t_std = compute_msd(np.array(dr),dshift)

msd_cor,msd_cor_std = compute_msd(np.array(dr_cor),dshift)

#---------------------------------------------------------------------------------------#
#				  New drift
#---------------------------------------------------------------------------------------#
msd_drift = []
msd_drift_std = []
confi_up = []
confi_down = []

plt.figure()
for shift_drift in dshift:
    dr_drift2 = []
    dr_drift_select=[]
    for ii,time in enumerate(tt[:-shift_drift]):
        int_g = 0.0
        for jj in range(shift_drift-1):
            int_g += dt_mean*(dr[ii+jj]- dr[ii])  
        dr_drift2.append( (dr[ii+shift_drift] - dr[ii] - 2.0/(shift_drift*1.0*dt_mean)*int_g )**2.0 )
    dr_drift2_mean = np.mean(dr_drift2)
    dr_drift2_std = np.std(dr_drift2, ddof=1)
    for xxx in dr_drift2:
        if xxx >  (dr_drift2_mean - 0.2 * dr_drift2_std) and xxx <  (dr_drift2_mean +  1.1 * dr_drift2_std):
            dr_drift_select.append(xxx)

    plt.plot([shift_drift*dt_mean for s in  dr_drift2],dr_drift2,'o')
    msd_drift.append( np.mean(dr_drift2) )
    msd_drift_std.append( np.std(dr_drift2, ddof=1))
    confi_up.append(dr_drift2_mean + dr_drift2_std*1.96/(len(dr_drift2)**0.5))
    confi_down.append(dr_drift2_mean - dr_drift2_std*1.96/(len(dr_drift2)**0.5))
plt.plot(tau,msd_drift)




for kk,ss in enumerate( msd_drift ):
    print msd_t[kk]/tau[kk]/20.0,ss/tau[kk]/20.0*3.0,confi_up[kk]/tau[kk]/20.0*3.0,confi_down[kk]/tau[kk]/20.0*3.0

plt.ylim(0,50)


#---------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------#
#				  Cut into n segment
#---------------------------------------------------------------------------------------#

n_segment = [s*1 for s in range(1,10)]
tau_segment = []
msd_segment = []
msd_segment_std = []


plt.figure()

for n in n_segment:
    dr_segment = []
    n_interval = ( len(tt)-1 )/n
    tau_segment.append(n_interval*dt_mean)
    for ii in range(n): 
        int_g = 0.0
        for jj in range(n_interval):
            int_g += dt_mean*(dr[ii*n_interval+jj] - dr[ii*n_interval] )  
        dr_segment.append( (dr[(ii+1)*n_interval-1] - dr[ii*n_interval] - 2.0/(n_interval*1.0*dt_mean)*int_g  )**2.0 )   
    plt.plot([n_interval*dt_mean for s in  dr_segment],dr_segment,'o')
    msd_segment.append( np.mean(dr_segment) )
    msd_segment_std.append( np.std(dr_segment, ddof=1))     

#for kk in range( len(msd_segment)-1 ):
#    print (msd_segment[kk+1]-msd_segment[kk])/(tau_segment[kk+1]-tau_segment[kk])/20.0*3.0
#---------------------------------------------------------------------------------------#

plt.figure()
plt.plot(tt,dr_cor)
#plt.xlim(0,10)

plt.figure()
plt.errorbar(tau,msd_t, yerr = msd_t_std)
plt.errorbar(tau,msd_drift, yerr = msd_drift_std)
#plt.errorbar(tau_segment,msd_segment, yerr = msd_segment_std)


plt.figure()
plt.plot(tau,msd_cor)

plt.show()



Cv_out = np.zeros(len(tt)*4)
Cv_out = Cv_out.reshape(len(tt),4)


for s in range(len(tt)):
    Cv_out[s][0] = tt[s]
    Cv_out[s][1] = dr[s]
    Cv_out[s][2] = tau[s] if s<len(dshift) else None
    Cv_out[s][3] = msd_cor[s] if s<len(dshift) else None


str = raw_input("Enter the file temperature (K): ")

str+= '.txt'

f = open(str,'w')
f.write(" #time(ns)  #displacement(A)      #tau(ns)           #MSD(A^2)\n")
for s in range(len(tt)):
    f.write('%f  %f  %f  %f\n'%(Cv_out[s][0],Cv_out[s][1],Cv_out[s][2],Cv_out[s][3]))















