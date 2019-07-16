import os
import glob
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import xml.dom.minidom


class point:
    def __init__(self,tag,r,proj,Eint):
        self.tag = tag
        self.r = [s for s in r]
        self.proj = proj
        self.Eint = Eint

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

    def average(self):
        x_average, y_average, z_average, p_average = 0.0, 0.0, 0.0, 0.0
        for s in self.plist:
            x_average += s.r[0] 
            y_average += s.r[1] 
            z_average += s.r[2] 
            p_average += s.proj
        x_average /= len(plist)
        y_average /= len(plist)
        z_average /= len(plist)
        p_average /= len(plist)
        return point(99,[x_average,y_average,z_average],p_average,0.0)


def compute_msd(position, shifts):
    msds = np.zeros(len(shifts))
    
    for i,shift in enumerate(shifts):
        diffs = position[:-shift if shift else None] - position[shift:]
        sqdlist = np.square(diffs)
        msds[i] = sqdlist.mean()
    return msds

#---------------------------------------------------------------------------------------#
#			   Get name of inputfile
#---------------------------------------------------------------------------------------#
inpu_name = raw_input('Please give the name of the input file (incluede_):')

res_in = 'res'+inpu_name

#---------------------------------------------------------------------------------------#
#				  Get time
#---------------------------------------------------------------------------------------#

tt = []    #in ns
tt.append(0.0)    #time zero

f_name = res_in + '/SIGEPS'
fcon = open(f_name)

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

'''
#---------------------------------------------------------------------------------------#
#				  Get initial position
#---------------------------------------------------------------------------------------#
fplist = []    #list to all plist for different steps

f_name = 'graph' + inpu_name + '.xml'

dom = xml.dom.minidom.parse(f_name)

root = dom.documentElement

itemlist = root.getElementsByTagName('object')

plist = []

for p_number in range(len(itemlist)):
    item = itemlist[p_number]
    plist.append(point(float(p_number),[float(item.getAttribute('cx')),float(item.getAttribute('cy')),float(item.getAttribute('cz'))],float(item.getAttribute('cz'))))



plist.append( step_plist(plist).average() )
fplist.append(step_plist(plist))
'''

#---------------------------------------------------------------------------------------#
#				  Get points
#---------------------------------------------------------------------------------------#
fplist = []    #list to all plist for different steps


f_name = res_in + '/COP*.dat'

flist = glob.glob(f_name)
flist.sort(key = str.lower)

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
            plist.append(point(int(line[0]),[float(s) for s in line[2:5]],float(line[5]),1.e6*1.e-30/1.6e-19*float(line[6])))
            line = fcon.readline()
    
    plist.append( step_plist(plist).average() )
    fplist.append(step_plist(plist))

    if ( len(plist) != len(fplist[0].plist) ) :
        print("The number of loops is chaning in file %s " %ff )


dr = np.zeros([len(tt),len(fplist[0].plist)])
dr2 = np.zeros(dr.shape)


for step in range(len(tt)):
    for ii in range(len(fplist[0].plist)):
      dr[step,ii] = fplist[step].plist[ii].proj - fplist[0].plist[ii].proj
      dr2[step,ii] = (fplist[step].plist[ii].proj - fplist[0].plist[ii].proj)**2.0


#---------------------------------------------------------------------------------------#
#	 Drift method in: S. L. Dudarev, C. R. Phys. 9, 409 (2009).
#---------------------------------------------------------------------------------------#
dt_mean = 0.0
sum_z = np.zeros(len(fplist[0].plist))

for step in range(1,len(tt)):
    dt_mean += tt[step]-tt[step-1]
    for ii in range(len(fplist[0].plist)):
        sum_z[ii] += (tt[step]-tt[step-1])*dr[step][ii]

dt_mean /= 1.0*(len(tt)-1.0)
print dt_mean


dr_cor = np.zeros(dr.shape)
dr_cor[0,:] = dr[0,:]

for step in range(1,len(tt)):
    for ii in range(len(fplist[0].plist)):
        dr_cor[step,ii] = dr[step,ii] - 2.0*tt[step]/tt[-1]/tt[-1]*sum_z[ii]     # R_drift(t) = R(t) - 2t/tau^2 * int_0^tau {R(t')dt'}  

 
#---------------------------------------------------------------------------------------#
#	                 Calculation of MSD cut in nt pieces
#---------------------------------------------------------------------------------------#

nt = 8.001

dshift = [10,30,50]
tau = np.array([s*dt_mean for s in dshift])

msd_cor_total = np.zeros([int(nt),len(fplist[0].plist),len(dshift)])
t_D = []
D_t = np.zeros([len(fplist[0].plist),int(nt)])

icompt = 1

for jj,kk in enumerate(tt):
    if kk >= tt[-1]/nt*1.0*icompt:
        msd_cor = np.zeros([len(fplist[0].plist),len(dshift)])  

   
        for ii in range(len(fplist[0].plist)):  # ii th loop
            msd_cor[ii,:] = compute_msd(dr[:jj,ii],dshift)
        for ii in range(len(fplist[0].plist)):  # ii th loop
            msd_cor_total[icompt-1,ii,:] = msd_cor[ii,:]
        for ii in range(len(fplist[0].plist)):  # ii th loop
            msd_cor[ii,:] /= (2.0*tau[:])       #A^2/ns  D = MSD/(2*tau)
        D_t[:,icompt-1] = np.mean(msd_cor,axis=1)[:]
        print msd_cor
        t_D.append(kk)
        icompt += 1





'''
#---------------------------------------------------------------------------------------#
#	                 Calculation of MSD total in 1 piece
#---------------------------------------------------------------------------------------#

dshift = [5,10,20,25,50,55,75,80,90,100,150]
tau = np.array([s*dt_mean for s in dshift])

msd_t = np.zeros([len(fplist[0].plist),len(dshift)])
msd_cor = np.zeros(msd_t.shape)

for ii in range(len(fplist[0].plist)):
    msd_t[ii,:] = compute_msd(dr[:,ii],dshift)
    msd_cor[ii,:] = compute_msd(dr_cor[:,ii],dshift)
'''

#---------------------------------------------------------------------------------------#
#	                     Single D for reference at 600K
#---------------------------------------------------------------------------------------#
# for simulation case
D0 = 128.25    #A^2/ns  D = MSD/(2*tau)
# for real case
D00 = D0 *3    #A^2/ns


ref_D0 = []
ref_D00 = []

for s in tau:
    ref_D0.append(D0*2.0*s)
    ref_D00.append(D00*2.0*s)

#---------------------------------------------------------------------------------------#
#	                    Calculating velocity
#---------------------------------------------------------------------------------------#
velocity_loop1 = []   #velocity of loop1

for i in range(len(tt)-1):
    velocity_loop1.append((fplist[i+1].plist[0].proj - fplist[i].plist[0].proj)/(tt[i+1]-tt[i]))  

mean_v = 0.0
for i,s in enumerate(tt):
    if s < 0.62:
        mean_v += velocity_loop1[i]*1.0


#---------------------------------------------------------------------------------------#
#	             Calculating distance between two loops
#---------------------------------------------------------------------------------------#
distance = []   #distance between two loops

for i,_ in enumerate(fplist):
    distance.append(abs(fplist[i].plist[1].proj - fplist[i].plist[0].proj))  


#---------------------------------------------------------------------------------------#
#	                     Plotting......
#---------------------------------------------------------------------------------------#
plt.figure()
plt.plot(tt,distance)

plt.figure()
plt.plot(tt[:-1],velocity_loop1)
plt.xlim(0,0.62)

plt.figure()
for ii in range(len(fplist[0].plist)):
    plt.plot(tt,dr[:,ii])

'''
for s in range(int(nt)):
    plt.figure()
    for ii in range(len(fplist[0].plist)):
        plt.plot(tau,msd_cor_total[s,ii,:])
    plt.plot(tau,ref_D0,'k--')
'''
plt.figure()
for ii in range(len(fplist[0].plist)):
    plt.plot(t_D,D_t[ii,:])

plt.plot(t_D,[D0 for s in t_D],'k--')


plt.figure()
plt.plot(tt,[s.plist[0].Eint for s in fplist])

plt.figure()
plt.plot(distance,[s.plist[0].Eint for s in fplist],'o')


'''
plt.figure()
for ii in range(len(fplist[0].plist)):
    plt.plot(tau,msd_cor[ii,:])
plt.plot(tau,ref_D0,'k--')



plt.figure()
for ii in range(len(fplist[0].plist)):
    plt.plot(tau,msd_t[ii,:])
plt.plot(tau,ref_D00,'k--')
'''
#---------------------------------------------------------------------------------------#
#	              Reference E_int
#---------------------------------------------------------------------------------------#
sys.path.append('/scratch/yli_scratch/m_loopp/doule/intE/')
import commu

plt.plot(commu.out1,commu.out2)    #ref_Eint
plt.plot(commu.out3,commu.out4)    #re_orient the normal of habit plan
plt.xlim(0,300)
plt.show()


#---------------------------------------------------------------------------------------#
#	                Output
#---------------------------------------------------------------------------------------#
flag = raw_input('Output? (y or n)')
if ( flag == 'n' or flag == 'N' ):
    sys.exit(0)


Cv_out = np.zeros(len(tt)*14)
Cv_out = Cv_out.reshape(len(tt),14)


for s in range(len(tt)):
    Cv_out[s][0] = tt[s]                     #time
    for ii in range(3):
        Cv_out[s][ii+1] = fplist[s].plist[ii].proj     #Z_cordinate
    Cv_out[s][4] = t_D[s] if s<len(t_D) else None      #t_cut for calculate D
    for ii in range(3):
        Cv_out[s][ii+5] = D_t[ii,s] if s<len(t_D) else None    #coeff_D
    Cv_out[s][8] = distance[s]                       #distance between 2 loop along z
    Cv_out[s][9] = fplist[s].plist[0].Eint           #intE from simulation
    Cv_out[s][10] = commu.out1[s] if s<len(commu.out1) else None            
    Cv_out[s][11] = commu.out2[s] if s<len(commu.out1) else None         #ref_intE
    Cv_out[s][12] = commu.out3[s] if s<len(commu.out3) else None            
    Cv_out[s][13] = commu.out4[s] if s<len(commu.out3) else None         #re_orient the normal of habit plan

outfile_name = raw_input("Enter the file name: ")

outfile_name+= '.txt'

f = open(outfile_name,'w')
f.write("#time(ns) #Z_cordinate_loop1(A) #Z_cordinate_loop2(A) #Z_cordinate_COM(A) #t_D(ns) #D_T_loop1(A^2/ns) #D_T_loop2(A^2/ns) #D_T_COM(A^2/ns) #distance_loop(A) #intE(eV) #zz(A) #intE_ref(eV) #zz(A) #intE_normal(eV)\n")
for s in range(len(tt)):
    for i in range(14):
        f.write('%f '%(Cv_out[s][i]))
    f.write(' \n')












