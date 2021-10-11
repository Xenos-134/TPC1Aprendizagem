
from re import T
from typing import ClassVar
import numpy
from numpy.core.arrayprint import IntegerFormat
from numpy.core.defchararray import array
from scipy.sparse import data

#Cross Validation 
from sklearn.model_selection import cross_val_score, StratifiedKFold



#Training algorithms
from  sklearn import neighbors
from sklearn import tree
import sklearn
from sklearn.naive_bayes import GaussianNB


# importing the required module
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')
from sklearn.utils import shuffle


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
        #fw.write(line)
        splitLine = line.split(",")
        type = 1 if (splitLine[9][:-1]=="benign") else 0
        dataLabels.append(type)
        trainingData.append([
            float(splitLine[0]),
            float(splitLine[1]),
            float(splitLine[2]),
            float(splitLine[3]),
            float(splitLine[4]),
            float(splitLine[5]),
            float(splitLine[6]),
            float(splitLine[7]),
            float(splitLine[8]),
        ])
        calculateClumptThickness(line)


#print("Our clump Thickness Bennign", clumpThickness)

x = numpy.linspace(0, 2 * numpy.pi, 400)
y = numpy.sin(x ** 2)

xz = numpy.arange(10)

fig, axs = plt.subplots(3, 3)
width = 0.35  # the width of the bars

#rects1 = ax.bar(xz - width/2, clumpThickness["B"], width, label='Benign')
#rects2 = ax.bar(xz + width/2, clumpThickness["M"], width, label='Malign')

axs[0,0].bar(xz - width/2, clumpThickness["B"], width, label='Benign')
axs[0,0].bar(xz + width/2, clumpThickness["M"], width, label='Malign')

axs[0,0].set_ylabel('Number of cases')
axs[0,0].set_title('Severity')
axs[0,0].set_xticks(xz)
axs[0,0].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[0,0].legend()

fig.tight_layout()

axs[0, 1].plot(x, y, 'tab:orange')
axs[0, 1].set_title('Axis [0, 1]')
axs[1, 0].plot(x, -y, 'tab:green')
axs[1, 0].set_title('Axis [1, 0]')
axs[1, 1].plot(x, -y, 'tab:red')
axs[1, 1].set_title('Axis [1, 1]')
axs[0, 2].plot(x, y)
axs[0, 2].set_title('Axis [0, 0]')
axs[2, 0].plot(x, y)
axs[2, 0].set_title('Axis [0, 0]')
axs[2, 2].plot(x, y)
axs[2, 2].set_title('Axis [0, 0]')
axs[1, 2].plot(x, y)
axs[1, 2].set_title('Axis [0, 0]')
axs[2, 1].plot(x, y)
axs[2, 1].set_title('Axis [0, 0]')

for ax in axs.flat:
    ax.set(xlabel='Number of Cases', ylabel='Class')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()


plt.show()

clf = neighbors.KNeighborsClassifier(n_neighbors=3, metric="euclidean")
clf.fit(trainingData, dataLabels)

#       >>kNN individual
#print(clf.predict([[5,1,1,1,2,1,3,1,1]]))
#print(clf.predict_proba([[5,1,1,5,2,1,4,1,10]]))


##Cross Fold Validation with 10 groups
knn_cv = neighbors.KNeighborsClassifier(n_neighbors=3)
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=107)
cv_scores = cross_val_score(knn_cv, trainingData, dataLabels, cv=cv)
print("Lista de Cv Scores :kNN",cv_scores)
print("Medium Cv: kNN", sum(cv_scores/10))

gnb = GaussianNB()
gaussianPred = cross_val_score(gnb, trainingData, dataLabels, cv=cv)
print("Gaussian score", sum(gaussianPred)/10)