from matplotlib import pyplot as plt
import numpy as np
import matplotlib.ticker
from dateutil import parser
import csv
import os
import pandas as pd
import math
import PySimpleGUI as sg

plt.rcParams.update({'font.size': 9}) #to make it more readable
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

xData = []
yData = []

i = 0

dataPlot = open("plotData.txt", "r")

for line in dataPlot:
    if any(c.isalpha() for c in line): continue

    myIntegers = [int(x) for x in line.split()]    
   
    if i % 2 == 0:
        yData.append(myIntegers)

    else:
        xData.append(myIntegers)

    i += 1


xDataJacubi = np.array(xData[0])
xDataGauss = np.array(xData[1])
xDataLU = np.array(xData[2])

yDataJacubi = np.array(yData[0])
yDataGauss = np.array(yData[1])
yDataLU = np.array(yData[2])

plt.xlabel('Sizes')
plt.ylabel('Times [ms]')
plt.title('Method comparison')

plt.plot(xDataJacubi, yDataJacubi, "-b", label="jacubi")

plt.plot(xDataGauss, yDataGauss, "-r", label = "gauss")

plt.plot(xDataLU, yDataLU, "-y", label = "LU decomposition")

plt.legend(loc="upper left")

plt.savefig('methodComparison.png')

