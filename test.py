import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import colors
from matplotlib import cm
import math


infile = open("test.txt",'r')
maxX = 50
maxY = 50

pop = np.zeros((maxX,maxY))

def xy2i(x,y,mx,my):
	return x*my+y

def i2xy(i,mx,my):
	return (i/my,i%my)

def wrap_around(x,w):
	return ((x%w)+w)%w

def disperse(x,y,nx,ny,mx,my):
	dx = x+nx
	dy = y+ny
	newX = wrap_around(dx,mx)
	newY = wrap_around(dy,my)
	return xy2i(newX,newY,mx,my)

def distDisperse(x,y,mx,my):
	a = np.random.uniform(0,2*math.pi)
	r = 4
	dX = math.floor(r*math.cos(a)+x+0.5)
	dY = math.floor(r*math.sin(a)+y+0.5)
	newX = wrap_around(dX,mx)
	newY = wrap_around(dY,my)
	pop[newX][newY] += 1
'''
for i in xrange(100000):
	distDisperse(10,10,maxX,maxY)
'''

for line in infile:
	xy = line.strip().split('\t')
	x = float(xy[0])
	y = float(xy[1])
	pop[x][y] += 1



fig,ax = plt.subplots(figsize=(5,5))
cax = ax.imshow(pop, interpolation='nearest',cmap=cm.Greys)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off',labelsize=6)
plt.tick_params(axis='y',which='both',right='off',left='off',labelleft='off',labelsize=6)
plt.savefig("test.png", dpi=1000)

plt.show()




