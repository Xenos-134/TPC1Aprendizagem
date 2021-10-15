
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
from sklearn.naive_bayes import MultinomialNB


# importing modules for plots
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')
#from sklearn.utils import shuffle



#################################################
#          Initialization of emty arrays for    #
# each parameter                                #
################################################
clumpThickness = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
cellSizeUniformity = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
cellShapeUniformity = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
marginalAdhesion = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
singleEpiCellSize = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
bareNuclei = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
blandChromatin = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
normalNucleoi = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}
mitoses = {"B":[0,0,0,0,0,0,0,0,0,0], "M": [0,0,0,0,0,0,0,0,0,0]}


#################################################
#       trainingData -Array that recives all    #
# filtered data from doc                        #
################################################
trainingData = []
dataLabels = []


#################################################
#   Functions that append data to resp array    #
################################################

def calculateClumptThickness(line):
    if(line.split(",")[9][:-1]=="benign"):
        clumpThickness["B"][int(line.split(",")[0])-1]+=1
        return
    clumpThickness["M"][int(line.split(",")[0])-1]+=1

def calculateCellSizeUniformity(line):
    if(line.split(",")[9][:-1]=="benign"):
        cellSizeUniformity["B"][int(line.split(",")[1])-1]+=1
        return
    cellSizeUniformity["M"][int(line.split(",")[1])-1]+=1

def calculateCellShapeUniformity(line):
    if(line.split(",")[9][:-1]=="benign"):
        cellShapeUniformity["B"][int(line.split(",")[2])-1]+=1
        return
    cellShapeUniformity["M"][int(line.split(",")[2])-1]+=1


def calculateMarginalAdhesion(line):
    if(line.split(",")[9][:-1]=="benign"):
        marginalAdhesion["B"][int(line.split(",")[3])-1]+=1
        return
    marginalAdhesion["M"][int(line.split(",")[3])-1]+=1


def calculateSingleEpiCellSize(line):
    if(line.split(",")[9][:-1]=="benign"):
        singleEpiCellSize["B"][int(line.split(",")[4])-1]+=1
        return
    singleEpiCellSize["M"][int(line.split(",")[4])-1]+=1


def calculateBareNuclei(line):
    if(line.split(",")[9][:-1]=="benign"):
        bareNuclei["B"][int(line.split(",")[5])-1]+=1
        return
    bareNuclei["M"][int(line.split(",")[5])-1]+=1


def calculateBlandChromatin(line):
    if(line.split(",")[9][:-1]=="benign"):
        blandChromatin["B"][int(line.split(",")[6])-1]+=1
        return
    blandChromatin["M"][int(line.split(",")[6])-1]+=1


def calculateNormalNucleoi(line):
    if(line.split(",")[9][:-1]=="benign"):
        normalNucleoi["B"][int(line.split(",")[7])-1]+=1
        return
    normalNucleoi["M"][int(line.split(",")[7])-1]+=1

def calculateMitoses(line):
    if(line.split(",")[9][:-1]=="benign"):
        mitoses["B"][int(line.split(",")[8])-1]+=1
        return
    mitoses["M"][int(line.split(",")[8])-1]+=1



#Open input file and filters lines by its content
#Appends to list of frequencies
f = open("breast.w.arff", "r")
fw = open("filteredData", "a")
for line in f:
    if(line.find("?")==-1 & line.find("@")==-1):
        #fw.write(line)
        splitLine = line.split(",")
        type = 1 if (splitLine[9][:-1]=="benign") else 0
        dataLabels.append(type)
        #appends to global arrray - used in kNN and Bayanesian prediction
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

        #appends data to parameter array - used in plot
        calculateClumptThickness(line)
        calculateCellSizeUniformity(line)
        calculateCellShapeUniformity(line)
        calculateMarginalAdhesion(line)
        calculateSingleEpiCellSize(line)
        calculateBareNuclei(line)
        calculateBlandChromatin(line)
        calculateNormalNucleoi(line)
        calculateMitoses(line)
        


    #spacing used for plot column
x = numpy.linspace(0, 2 * numpy.pi, 400)
y = numpy.sin(x ** 2)

xz = numpy.arange(10)

fig, axs = plt.subplots(3, 3)
width = 0.35  # the width of the bars


#################################################
#          Clump Thickness                     #
################################################
axs[0,0].bar(xz - width/2, clumpThickness["B"], width, label='Benign')
axs[0,0].bar(xz + width/2, clumpThickness["M"], width, label='Malign')
axs[0,0].set_xticks(xz)
axs[0,0].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[0,0].legend()
axs[0, 0].set_title('Clump Thickness')



#################################################
#          Cell Size Uniformity                #
################################################
axs[0,1].bar(xz - width/2, cellSizeUniformity["B"], width, label='Benign')
axs[0,1].bar(xz + width/2, cellSizeUniformity["M"], width, label='Malign')
axs[0,1].set_xticks(xz)
axs[0,1].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[0,1].legend()
axs[0,1].set_title('Cell Size Uniformity')


#################################################
#          Cell Shape  Uniformity              #
################################################
axs[0,2].bar(xz - width/2, cellShapeUniformity["B"], width, label='Benign')
axs[0,2].bar(xz + width/2, cellShapeUniformity["M"], width, label='Malign')
axs[0,2].set_xticks(xz)
axs[0,2].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[0,2].legend()
axs[0,2].set_title('Cell Shape Uniformity')


#################################################
#          Cell Marginal Adhe sion             #
################################################
axs[1,0].bar(xz - width/2, marginalAdhesion["B"], width, label='Benign')
axs[1,0].bar(xz + width/2, marginalAdhesion["M"], width, label='Malign')
axs[1,0].set_xticks(xz)
axs[1,0].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[1,0].legend()
axs[1,0].set_title('Cell Marginal Adhesion')


#################################################
#          Cell Epicell size                   #
################################################
axs[1,1].bar(xz - width/2, singleEpiCellSize["B"], width, label='Benign')
axs[1,1].bar(xz + width/2, singleEpiCellSize["M"], width, label='Malign')
axs[1,1].set_xticks(xz)
axs[1,1].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[1,1].legend()
axs[1,1].set_title('Cell Single Epicell Size')


#################################################
#          Cell Bare Nuclei                    #
################################################
axs[1,2].bar(xz - width/2, bareNuclei["B"], width, label='Benign')
axs[1,2].bar(xz + width/2,bareNuclei["M"], width, label='Malign')
axs[1,2].set_xticks(xz)
axs[1,2].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[1,2].legend()
axs[1,2].set_title('Cell Bare Nuclei')


#################################################
#          Cell Bland Chromatin                 #
################################################
axs[2,0].bar(xz - width/2, blandChromatin["B"], width, label='Benign')
axs[2,0].bar(xz + width/2,blandChromatin["M"], width, label='Malign')
axs[2,0].set_xticks(xz)
axs[2,0].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[2,0].legend()
axs[2,0].set_title('Cell Bland Chromatin')


#################################################
#          Cell Normal Nucleoi                 #
################################################
axs[2,1].bar(xz - width/2, normalNucleoi["B"], width, label='Benign')
axs[2,1].bar(xz + width/2, normalNucleoi["M"], width, label='Malign')
axs[2,1].set_xticks(xz)
axs[2,1].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[2,1].legend()
axs[2,1].set_title('Normal Nucleoi')

#################################################
#          Cell Normal Nucleoi                 #
################################################
axs[2,2].bar(xz - width/2, mitoses["B"], width, label='Benign')
axs[2,2].bar(xz + width/2, mitoses["M"], width, label='Malign')
axs[2,2].set_xticks(xz)
axs[2,2].set_xticklabels([1,2,3,4,5,6,7,8,9,10])
axs[2,2].legend()
axs[2,2].set_title('Mitoses')


for ax in axs.flat:
    ax.set(xlabel='Severity', ylabel='Number of cases')

fig.tight_layout()
plt.show()

clf = neighbors.KNeighborsClassifier(n_neighbors=3)
clf.fit(trainingData, dataLabels)

#       >>kNN individual
#print(clf.predict([[5,1,1,1,2,1,3,1,1]]))
#print(clf.predict_proba([[5,1,1,5,2,1,4,1,10]]))


#Cross Fold Validation with 10 groups
#################################################
#                  QUESTION 6                   #
################################################
print("*************Question number 6*************")
for i in [3, 5, 7]:
    knn_cv = neighbors.KNeighborsClassifier(n_neighbors=i, metric="euclidean")
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=107)
    cv_scores = cross_val_score(knn_cv, trainingData, dataLabels, cv=cv)
    print("Medium CV: kNN for n_neighbors = " + str(i) + ": " + str(sum(cv_scores/10)))

print("Best number of neighbors is 5")

#################################################
#                  QUESTION 7                   #
################################################
print("*************Question number 7*************")

knn_cv = neighbors.KNeighborsClassifier(n_neighbors=3, metric="euclidean")
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=107)
cv_scores = cross_val_score(knn_cv, trainingData, dataLabels, cv=cv)
print("Medium CV: kNN", sum(cv_scores/10))
gnb = MultinomialNB()
gaussianPred = cross_val_score(gnb, trainingData, dataLabels, cv=cv)
print("Gaussian score", sum(gaussianPred)/10)

print("-->We conclude that KNN is better for this particular dataset")

#################################################
#                  QUESTION 8                  #
################################################
print("*************Question number 8*************")

print("-->Naive Bayes assumes that all predictors (or features) are independent, rarely happening in real life.")
print("-->This algorithm faces the 'zero-frequency problem' where it assigns zero probability to a categorical variable whose category in the test data set wasn't available in the training dataset.")