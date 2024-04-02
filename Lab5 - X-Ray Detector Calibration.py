import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.ndimage.morphology import binary_hit_or_miss

#%% ######### INIT ##########

# Path
folder = r"C:/Users/eduar/Documents/Faculdade/Mestrado IST/2º ano/Space Engineering/Laboratories/5. X-Ray Detectors"

## Integral image
imageIntegral = np.zeros((256,256))

## Energy
energy = np.zeros((256,256))

## Single Pixel Event Pattern
singlePixelEvent =      [[0,0,0],
                        [0,1,0],
                        [0,0,0]]

#%% ######### LOAD ##########

## Load config 
a = np.loadtxt(os.path.join(folder,"config","calib_a.txt"))
b = np.loadtxt(os.path.join(folder,"config","calib_b.txt"))
c = np.loadtxt(os.path.join(folder,"config","calib_c.txt"))
t = np.loadtxt(os.path.join(folder,"config","calib_t.txt"))

## Find data
dataList = glob(os.path.join(folder,"Data","*.txt"))

## Load data
image = np.loadtxt(dataList[0])
# image = np.loadtxt(data)

## get non-zero element
idx = np.nonzero(image)

## Calibration
imageCalib = c

# Get list of energies from frame
energy = np.concatenate((energy,imageCalib))

# integral image 
imageIntegral += np.where(image > 0, 1, 0)

# Calculate the average of the entire array
a_ = np.mean(a)
b_ = np.mean(b)
c_ = np.mean(c)
t_ = np.mean(t)

def energy_x(t, a, y, b, c):
    x = (t*a + y - b + np.sqrt((b + t*a - y)**2 + 4*a*c)) / (2*a)
    return x

X = energy_x(t_, a_, y_, b_, c_)
Y = 

#%% ######### PRINT / PLOT ##########

plt.figure()
plt.imshow(imageIntegral, cmap='gray_r')
plt.xlabel("x (px)")
plt.ylabel("y (px)")
plt.colorbar(label="Number of counts (--)")


plt.figure()
plt.grid(True, which = 'both', alpha=0.3)
plt.grid(b=True, which='major', linestyle='-')
plt.grid(b=True, which='minor', linestyle='--')
plt.minorticks_on()

plt.plot(X,Y,".")
plt.ylabel("Number of counts (--)")
plt.xlabel("Energy (keV)")

plt.show()


# print("Maximum peak energy: {} keV".format(bin_edges[np.argmax(hist)]))
