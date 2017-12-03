import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import sys

dat_path = "./data/"
plot_path = "./plots/"
allfiles = listdir(dat_path)
allfiles_sorted = sorted(allfiles)

nx = int(sys.argv[1])
ny = int(sys.argv[2])
prop = int(sys.argv[3])

if prop>5 or prop<1 :
    prop = 5

def choose_prop(x):
    return {
        1 : 'u_velocity',
        2 : 'v_velocity',
        3 : 'pressure',
        4 : 'density',
        5 : 'Temperature'
    }.get(x,'Unknown')

proprty = choose_prop(prop)

frame_no=0
for fil_name in allfiles_sorted:
    fname = dat_path + fil_name
    data = np.loadtxt(fname)
    Z = np.reshape(data[:,prop+1],(ny,nx))
    X = np.reshape(data[:,0],(ny,nx))
    Y = np.reshape(data[:,1],(ny,nx))
    plt.title(proprty,loc='left')
    plt.title(fil_name)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.contour(X,Y,Z,30)
    plt.colorbar()
    plt.savefig( plot_path + proprty + str(frame_no).zfill(6) +'.png')
    plt.clf()
    frame_no = frame_no +1
    if(frame_no % 10==0):
        sys.stdout.write('\rPlotting file '+str(frame_no)+'/'+str( len(allfiles) ) )
        sys.stdout.flush()
print("\nDone Plotting\n!!!!!!---------#########################---------!!!!!!")
