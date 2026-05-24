import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
import numpy as np
import math
import gmpy2

sqrt = math.sqrt
floor = math.floor
ceil = math.ceil

def points_old(n):
    # Create a sphere
    r = 1
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    Sx = r*sin(phi)*cos(theta)
    Sy = r*sin(phi)*sin(theta)
    Sz = r*cos(phi)

    # Get coordinates of sums of squares
    xx = np.array([])
    yy = np.array([])
    zz = np.array([])
    
    for x in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
        for y in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
            rem = n - x**2 - y**2
            if gmpy2.is_square(rem):
                if rem != 0:
                    z = - sqrt(rem)
                    xx = np.append(xx, x/sqrt(n))
                    yy = np.append(yy, y/sqrt(n))
                    zz = np.append(zz, z/sqrt(n))
                z = sqrt(rem)
                xx = np.append(xx, x/sqrt(n))
                yy = np.append(yy, y/sqrt(n))
                zz = np.append(zz, z/sqrt(n))

    #Set colours and render
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    ax.scatter(xx,yy,zz,color="k",s=3)
    ax.plot_surface(Sx, Sy, Sz,  rstride=1, cstride=1, color='r', alpha=1, linewidth=0)

    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("auto")
    plt.tight_layout()
    return plt, xx, yy, zz


def points(n):
    # Create a sphere
    r = 1.0
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:100j, 0:2 * pi:100j]
    Sx = r*sin(phi)*cos(theta)
    Sy = r*sin(phi)*sin(theta)
    Sz = r*cos(phi)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    mlab.clf()

    # Get coordinates of sums of squares
    xx = np.array([])
    yy = np.array([])
    zz = np.array([])
    
    for x in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
        for y in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
            rem = n - x**2 - y**2
            if gmpy2.is_square(rem):
                if rem != 0:
                    z = - sqrt(rem)
                    xx = np.append(xx, x/sqrt(n))
                    yy = np.append(yy, y/sqrt(n))
                    zz = np.append(zz, z/sqrt(n))
                z = sqrt(rem)
                xx = np.append(xx, x/sqrt(n))
                yy = np.append(yy, y/sqrt(n))
                zz = np.append(zz, z/sqrt(n))

    # Create figure
    mlab.mesh(Sx, Sy, Sz, color = (1,0.0,0.0))
    mlab.points3d(xx, yy, zz, scale_factor=0.05)

    print(xx.size)
    mlab.show()


def count(n):
    res = 0
    for x in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
        for y in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
            rem = n - x**2 - y**2
            if gmpy2.is_square(rem):
                if rem == 0:
                    res = res + 1
                else:
                    res = res + 2

    return res


def largest(a, b):
    index = a
    current = 0
    for n in range(a, b+1):
        if count(n) > current:
            index = n
            current = count(n)

    return index, current


def all_points(high, low=1, residues=[1], m=1):
    # Create a sphere
    r = 1.0
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:100j, 0:2 * pi:100j]
    Sx = r*sin(phi)*cos(theta)
    Sy = r*sin(phi)*sin(theta)
    Sz = r*cos(phi)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    mlab.clf()

    # Get coordinates of sums of squares
    xx = np.array([])
    yy = np.array([])
    zz = np.array([])

    for n in range(low, high+1):
        isres = False
        for res in residues:
            if (n%m) == (res%m):
                isres = True
                break

        if not isres and residues != []:
            continue
        
        for x in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
            for y in range(int(floor(-sqrt(n))), int(floor(sqrt(n))) + 1):
                rem = n - x**2 - y**2
                if gmpy2.is_square(rem):
                    if rem != 0:
                        z = - sqrt(rem)
                        xx = np.append(xx, x/sqrt(n))
                        yy = np.append(yy, y/sqrt(n))
                        zz = np.append(zz, z/sqrt(n))
                    z = sqrt(rem)
                    xx = np.append(xx, x/sqrt(n))
                    yy = np.append(yy, y/sqrt(n))
                    zz = np.append(zz, z/sqrt(n))

    # Create figure
    mlab.mesh(Sx, Sy, Sz, color = (1,0.0,0.0))
    mlab.points3d(xx, yy, zz, scale_factor=0.02)

    print(xx.size)
    mlab.show()
        
        
                    
                    
    
