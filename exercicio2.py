
from re import T
from typing import ClassVar
import numpy
from numpy.core.arrayprint import IntegerFormat
from numpy.core.defchararray import array

#Training algorithms
from  sklearn import neighbors
from sklearn import tree


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

trainingData = []
dataLabels = []


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
        clumpThickness["B"][int(line.split(",")[3])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[3])-1]+=1


def calculateSingleEpiCellSize(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[4])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[4])-1]+=1


def calculateBareNuclei(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[5])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[5])-1]+=1


def calculateBlandChromatin(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[6])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[6])-1]+=1


def calculateMitoses(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[7])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[7])-1]+=1


#Open input file and filters lines by its content
#Appends to list of frequencies
f = open("breast.w.arff", "r")
fw = open("filteredData", "a")
for line in f:
    if(line.find("?")==-1 & line.find("@")==-1):
        #Depois fazer write de sample organizado
        fw.write(line)
        splitLine = line.split(",")
        type = 1 if (splitLine[9][:-1]=="benign") else 0
        dataLabels.append(type)
        trainingData.append([
            splitLine[0],
            splitLine[1],
            splitLine[2],
            splitLine[3],
            splitLine[4],
            splitLine[5],
            splitLine[6],
            splitLine[7],
            splitLine[8],
        ])
        calculateClumptThickness(line)


print("Our clump Thickness Bennign", clumpThickness)



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

#plt.show()
print(dataLabels)

clf = neighbors.KNeighborsClassifier(n_neighbors=3)
clf.fit(trainingData, dataLabels)
print(clf.predict([[5,1,1,1,2,1,3,1,1]]))
print(clf.predict_proba([[5,1,1,5,2,1,4,1,10]]))