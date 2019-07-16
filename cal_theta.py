import os
import glob
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import xml.dom.minidom


class point:
    def __init__(self,tag,r):
        self.tag = tag
        self.r = [s for s in r]


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

def vec(p0,p1):
    vect = []
    for i in range(3):
          vect.append(p0.r[i] - p1.r[i])
    return vect

def prod_vec(v1,v2):
    vect = []
    vect.append(v1[1]*v2[2] - v1[2]*v2[1])
    vect.append(- v1[0]*v2[2] + v1[2]*v2[0])
    vect.append(v1[0]*v2[1] - v1[1]*v2[0])
    return vect

def theta(v1,v2):
    product_s = 0.0
    nor_v1 = 0.0
    nor_v2 = 0.0
    for i in range(3):
        nor_v1 += v1[i]*v1[i]
        nor_v2 += v2[i]*v2[i]
        product_s += v1[i]*v2[i]

    product_s = product_s/ (nor_v1**0.5) / (nor_v2**0.5)


    return math.acos(product_s)


#---------------------------------------------------------------------------------------#
#				  Get Angle
#---------------------------------------------------------------------------------------#
f_name = 'res_xx120_zz120/DDD*.vtk'

flist = glob.glob(f_name)
flist.sort(key = str.lower)

tilt_angle = []
tilt_angle2 = []

for ff in flist:

    plist= []

    tag_vec = []

    fcon = open(ff)

    line = fcon.readline()

    while line : 


        line = line.strip().split()
    
        if line == []:
            line = fcon.readline()
            continue


        elif line[0] == 'POINTS':
            num_point = int(line[1])
            for i in range(num_point):
                line = fcon.readline()
                line = line.strip().split()
                plist.append(point(i,[float(line[0]), float(line[1]), float(line[2])]))
            line = fcon.readline() 
            continue

        elif line[0] == 'LINES':
            num_tag = int(line[1])
            for i in range(num_tag):
                line = fcon.readline()
                line = line.strip().split()
                tag_vec.append(int(line[1]))
            line = fcon.readline()
            break
        
        else:
            line = fcon.readline()
            continue
    


    tag0 = tag_vec[0]
    tag1 = tag_vec[1]
    tag2 = tag_vec[2]

    tag3 = tag_vec[3]
    tag4 = tag_vec[4]


    vect1 = vec(plist[tag0],plist[tag1])
    vect2 = vec(plist[tag0],plist[tag2])

    vect3 = vec(plist[tag1],plist[tag3])
    vect4 = vec(plist[tag1],plist[tag4])   

    n1 = [0,0,1]

    n2 = prod_vec(vect1,vect2)
    n3 = prod_vec(vect3,vect4)


    thetav1v2 = theta(n1,n2)/math.pi*180.0
    if thetav1v2 > 90.:
        thetav1v2 = 180. - thetav1v2

    tilt_angle.append(thetav1v2)

    thetav3v4 = theta(n1,n3)/math.pi*180.0
    if thetav3v4 > 90.:
        thetav3v4 = 180. - thetav3v4

    tilt_angle2.append(thetav3v4)

plt.plot([s for s in range(len(tilt_angle))],tilt_angle)
plt.plot([s for s in range(len(tilt_angle2))],tilt_angle2)
plt.xlim([0,1000])
plt.show()






