#!usr/bin/python
#*******************************************************************************
#	This file is used to calculate the Interaction energy	      
#*******************************************************************************

import os
import os.path
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.interpolate import spline
from pylab import *
import math 
		
#*******************************************************************************
#				Parameters				
#*******************************************************************************

xx = 120

miu = 0.393  #0.511   # 0.393   #  # eV/A^3
poison = 0.43  #0.29   # 0.43    #0.275
aa = 2.856 # A
radius = 45  #15   # 45 # A


burger_s = aa/2.0 * aa/2.0 * 3.0
area = radius*radius * 3.0/ 2.0 * 3.0**0.5

drift = 1.6  #1.6 for xx = 120 middle and 1.6 for final
factor = 1.28  # 1.25 for xx = 120 middle and 1.2 for final if take 0.29 v  # 1.3for xx = 120 middle and 1.16 for final if take 0.43 v   # factor = 0.0 for mu = 101MPa if take 0.275 v

#*******************************************************************************
#				Vectors				
#*******************************************************************************

v_x = [-1.0/(6.0**0.5),-1.0/(6.0**0.5),2.0/(6.0**0.5)]
v_y = [1.0/(2.0**0.5),-1.0/(2.0**0.5),0.0]
v_z = [1.0/(3.0**0.5),1.0/(3.0**0.5),1.0/(3.0**0.5)]

a = -0.2   #for 120x_120y from [-0.2,0.3,1.0] to [2.0,0.0,1.0]
b = 0.3
c = 1.0
n1 = [c*v_z[i] + b*v_y[i] + a*v_x[i] for i in range(3)]
n2 = [c*v_z[i] + b*v_y[i] + a*v_x[i] for i in range(3)]

n1_normal = 0.0
n2_normal = 0.0
for i in range(3):
    n1_normal += n1[i]**2.0
    n2_normal += n2[i]**2.0
n1_normal = n1_normal**0.5
n2_normal = n2_normal**0.5

n_1 = [s/n1_normal for s in n1]
n_2 = [s/n2_normal for s in n2]
print n_1


#*******************************************************************************
#			   Class of intE from simulation			
#*******************************************************************************

class point:
    def __init__(self,zcor,Eint):
        self.zcor = zcor
        self.Eint = Eint




#*******************************************************************************
#			   Theoritical Interaction Energy 			
#*******************************************************************************
def E_int(z,r):
    cos_theta = z/r*1.0
    E12 = miu*burger_s*area*area/4.0/math.pi/(1.0-poison) * (15.0*cos_theta**4.0 - 6.0*cos_theta**2.0 - 1.0) /r/r/r 
    return E12

def E_int_nr(v_r):
    r = 0.0
    for i in range(3):
        r += v_r[i]**2.0 
    r = r**0.5

    r_unit = []
    for i in range(3):
        r_unit.append(v_r[i]/r)
    
    n1n2 = 0.0
    n1r = 0.0
    n2r = 0.0
    b1r = 0.0
    b1b2 = 0.0
    b1n1 = 0.0
    b1n2 = 0.0
    for i in range(3):
        n1n2 += n_1[i]*n_2[i]
        n1r += n_1[i]*r_unit[i]
        n2r += n_1[i]*r_unit[i]
        b1r += v_z[i]*r_unit[i]
        b1b2 += v_z[i]*v_z[i]
        b1n1 += v_z[i]*n_1[i]
        b1n2 += v_z[i]*n_2[i]  
    b2n1 = b1n1
    b2n2 = b1n2 
    b2r = b1r

  
#    temp = 15.0 * n1r**2.0 * n2r**2.0 - (4.0*poison - 1.0) - 12.0*poison*n1n2*n1r*n2r - (1.0-2.0*poison)*(3.0*n1r**2.0 + 3.0*n2r**2.0 + 2.0*n1n2**2)  for prismatic
    temp = 15.0*b1r*b2r*n1r*n2r - 3.0*poison*(b1b2*n1r*n2r+b1n2*n1r*b2r+b2n1*n2r*b1r+n1n2*b1r*b2r) - 3.0*(1.0-2.0*poison)*b1n1*b2r*n2r-3.0*(1.0-2.0*poison)*b2n2*b1r*n1r - (1.0-2.0*poison)*(b1b2*n1n2+b1n2*b2n1)-(4.0*poison-1.0)*b1n1*b2n2
    E12 = miu*burger_s*area*area/4.0/math.pi/(1.0-poison) /r/r/r * temp
    return E12*factor


#---------------------------------------------------------------------------------------#
#				  Get points
#---------------------------------------------------------------------------------------#

f_name = 'res_xx' + str(int(xx)) + '_*/COP*.dat'

flist = glob.glob(f_name)
plist = []

for ff in flist:
    fcon = open(ff)
    line = fcon.readline()    #skip first line
    line = fcon.readline()
    line = fcon.readline()    # direct to 3rd line
    
    while line : 
        line = line.strip()
        if line == "":
            line = fcon.readline()
            continue

        else:
            line = line.split()
            plist.append(point(float(line[5]),float(line[6])))
            plist[-1].Eint =  plist[-1].Eint*1.e6*1.e-30/1.6e-19
            line = fcon.readline()
    

#---------------------------------------------------------------------------------------#
#		       Calculating relaxition time
#---------------------------------------------------------------------------------------#
F_plot = []

plist = sorted(plist,key=lambda point: point.zcor)
for i in range(1,len(plist)):
    F = -(plist[i].Eint - plist[i-1].Eint) / (plist[i].zcor - plist[i-1].zcor)
    F_plot.append(F)

Eint_cal = [s.Eint for s in plist]
z_min = Eint_cal.index(min(Eint_cal))

F_mean = 0.0
n_compt = 0.0
t_relax = 0.0
'''
for i in range(len(F_plot)):
    if ( plist[i].zcor > z_min and plist[i].zcor <= xx ):
        F_mean += F_plot[i]*1.0
        n_compt += 1.0

F_mean /= n_compt
t_relax = F_mean*1.6e-19/1.e6/1.e-30/0.08/6.0/4.5    #A/ns
print t_relax
t_relax = (xx-z_min)/2.0/t_relax #ns
print t_relax
'''
#*******************************************************************************
#				Other cases of Eint
#*******************************************************************************	

zz = np.linspace(0,500,5001)

E_plot = []



xx -= 6.0
for s in zz:
    rr = (xx*xx + s*s)**0.5
    E_plot.append(E_int(s,rr))



E_plot_nr = []

for s in zz:
    v_r = []
    for i in range(3):
        v_r.append(xx*v_x[i] + 0.0*v_y[i] + s*v_z[i])
    E_plot_nr.append(E_int_nr(v_r))



#---------------------------------------------------------------------------------------#
#		       Calculating life-time
#---------------------------------------------------------------------------------------#
rad_condition = 2.3    #1.0   1.8  2.25 2.3

#z_cal = [s.zcor for s in plist]
#Eint_cal = [s.Eint/rad_condition/rad_condition/rad_condition/rad_condition for s in plist]

z_cal = [s for s in zz]
Eint_cal = [s/rad_condition/rad_condition/rad_condition/rad_condition for s in E_plot_nr]

z_min,z_max = Eint_cal.index(min(Eint_cal)),Eint_cal.index(max(Eint_cal))
dh = z_cal[1] - z_cal[0]   #A

B_drag = 0.08    #MPa*ns
k_B = 1.38e-23   #J/K
kbt = k_B*200.0/1.6e-19    #eV
L_loop = radius * 6.0 /rad_condition   #A

print z_min*dh,z_max*dh

if xx == 50:
    dE2_min = (Eint_cal[1] - 2.0* Eint_cal[0] + Eint_cal[1])/dh/dh
    Eb = Eint_cal[z_max] - Eint_cal[0]
else:
    dE2_min = (Eint_cal[z_min+1] - 2.0* Eint_cal[z_min] + Eint_cal[z_min-1])/dh/dh
    Eb = Eint_cal[z_max] - Eint_cal[z_min]     #eV

dE2_max = (Eint_cal[z_max+1] - 2.0* Eint_cal[z_max] + Eint_cal[z_max-1])/dh/dh    #eV/A2  


tau = (abs(dE2_min*dE2_max))**0.5/L_loop*160.21766208*1000  # 1eV/A3 = 160.21766208*1000 MPa
tau = tau/(B_drag*1.e-9)/2.0/math.pi*math.exp(-Eb/kbt)    # 1/s
tau = 1.0/tau     #s
# tau = tau/3600.0/24.0    #day
print tau


#*******************************************************************************
#				Plotting
#*******************************************************************************	

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(zz,E_plot,linewidth=1.0, label = 'Generally used with approximation')
ax.plot([s.zcor for s in plist], [s.Eint for s in plist],'o',label = 'Implemented in simulation')
ax.plot(zz,E_plot_nr,linewidth=1.0)
ax.set_xlabel('Displacement of loop 2 along its glide cylinder (Angstrom)')
ax.set_ylabel('E_int (eV)')



F_plot.append(F_plot[-1])   #eV/A
ax2 = ax.twinx()
ax2.plot([s.zcor for s in plist],[s*1.6e-19*1.e24 for s in F_plot],'k',linewidth=1.0)
ax2.set_xlabel('ZZ')
ax2.set_ylabel('F_int')

plt.show()


'''		
flag = raw_input('Output? (y or n)')
if ( flag == 'n' or flag == 'N' ):
    sys.exit(0)


Cv_out = np.zeros(len(zz)*2)
Cv_out = Cv_out.reshape(len(zz),2)

out3 = [s.zcor for s in plist]
out4 = [s.Eint for s in plist]

for s in range(len(zz)):
    Cv_out[s][0] = zz[s]          
    Cv_out[s][1] = E_plot_nr[s]        #ref_intE


outfile_name = raw_input("Enter the file name: ")

outfile_name+= '.txt'

f = open(outfile_name,'w')
f.write("#zz(A) #intE_nr(eV)\n")
for s in range(len(zz)):
    for i in range(2):
        f.write('%f '%(Cv_out[s][i]))
    f.write(' \n')

'''





