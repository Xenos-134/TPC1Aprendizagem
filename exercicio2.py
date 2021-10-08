
from typing import ClassVar
import numpy
from numpy.core.arrayprint import IntegerFormat
from numpy.core.defchararray import array

# importing the required module
import matplotlib.pyplot as plt


clumpThickness = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
cellSizeUniformity = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
cellShapeUniformity = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
marginalAdhesion = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
singleEpiCellSize = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
bareNuclei = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
blandChromatin = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
normalNucleoi = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
mitoses = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}



def calculateClumptThickness(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[0])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[0])-1]+=1

def calculateCellSizeUniformity(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[1])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[1])-1]+=1

def calculateCellShapeUniformity(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1


def calculateMarginalAdhesion(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1


def calculateSingleEpiCellSize(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1


def calculateBareNuclei(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1


def calculateBlandChromatin(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1


def calculateMitoses(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[2])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[2])-1]+=1

f = open("breast.w.arff", "r")
fw = open("filteredData", "a")
for line in f:
    if(line.find("?")==-1 & line.find("@")==-1):
        fw.write(line)
       
        #if(line.split(",")[9][:-1]=="benign"):clumpThickness["B"][int(line.split(",")[0])-1]+=1
        #else: clumpThickness["M"][int(line.split(",")[0])-1]+=1
        calculateClumptThickness(line)




print("Our clump Thickness Bennign", clumpThickness)

#plt.plot(clumpThicknessCountM)
#plt.plot(clumpThicknessCountB)

width = 0.35  # the width of the bars
x = numpy.arange(10)

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, clumpThickness["B"], width, label='Benign')
rects2 = ax.bar(x + width/2, clumpThickness["M"], width, label='Malign')




# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of cases')
ax.set_title('Severity')
ax.set_xticks(x)
ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10])
ax.legend()

fig.tight_layout()

plt.show()