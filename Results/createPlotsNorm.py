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

dataPlot = open("plotDataNorm.txt", "r")

for line in dataPlot:
    if any(c.isalpha() for c in line):
       if not any(c == '.' for c in line):
           continue

    myIntegers = [float(x) for x in line.split()]    
   
    if i % 2 == 0:
        yData.append(myIntegers)

    else:
        xData.append(myIntegers)

    i += 1


xDataJacubi = np.array(xData[0])
xDataGauss = np.array(xData[1])

yDataJacubi = np.array(yData[0])
yDataGauss = np.array(yData[1])

plt.xlabel('Iterations')
plt.ylabel('Norms')
plt.title('Norms for iterations')

plt.semilogy(xDataJacubi, yDataJacubi, "-b", label="Jacobi")

plt.semilogy(xDataGauss, yDataGauss, "-r", label = "Gauss")

plt.legend(loc="upper left")

plt.savefig('norms.png')

