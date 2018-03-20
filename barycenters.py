'''
    File name: barycenters2.py
    Author: Natalia Villegas Franco
    Date created: 9/30/2017
    Date last modified: 1/28/2018
    Python Version: 3.6.2
'''

import pandas as pd
import os
import numpy as np
import gmplot
import pickle
import subprocess
import re
import datetime
import calendar
from time import time
import sys

global t_init 

t_init = time()

## Things to set:
# - Set date range
y1=2016
y2=2016
m1=1
m2=12
# - Set the method.
method='exact'
#method='2aprox'
# - Group of months
gm=1
# - Tolerance taking the final police (zj>tolerance)
tolerance=0.0

# ----------------------------------- FUNCTIONS -----------------------------------

## Function which sum a number of months "months" to a datetime "date"
def add_months(date,months):
    month = date.month - 1 + months # Sum directly the months
    year = int(date.year + month / 12 ) # Compute number of years
    month = month % 12 + 1 # Compute number of months
    return datetime.date(year,month,1) # Return as a datetime

## Function which defines all the possible combinations to compute xj from xik
def combinations(months,number,it,result,variables):
    if (number < len(months)): # It is not in the last month: save the combinations and recursive calls
        for x in range(len(months[number])): 
            # Where the combinations of indexes are generated
            it[number]=x; 
            # Recursive function (generate all the combinations deeply)
            combinations(months,number+1,it,result,variables) # The first one called is the last one done!
    else: # The last month: no more recurseive calls and compute the result in this recursive call
        xik=[]
        xj=[]
        for idx,var in enumerate(variables): # In that case: lat and lon
            # Mean of all the combinations saved in the 'it' indeces array, and save in result array
            xik.append(np.array([months[i][var][it[i]] for i in range(len(it))]))
            xj.append(np.mean(xik[idx]))
            result[idx].append(xj[idx]) 
        result[2].append([it[i] for i in range(len(it))]);
        result[3].append(xik)
        #result[3].append(sum((xik[0]-xj[0])**2+(xik[0]-xj[0])**2))

## Function which calculates cikj that is the distances between xj and xik for all the possibilities
def cdistance(xik,xj,cikj):
    for i in range(len(xik)):
        cikj.append([])
        for k in range(len(xik[i])): 
            cikj[i].append([])
            for j in range(len(xj)):
                # cikj is saved in that indexes orden for AMPL reasons
                cikj[i][k].append((xik[i][k][0]-xj[j][0])**2+(xik[i][k][1]-xj[j][1])**2)

## Function which prints into the terminal and into the result.txt the given string p
def prints(p):
    print(p)
    fr.write(p+'\n')

# ----------------------------------- MAIN PROGRAM -----------------------------------

## Create the range of dates
# Define the start and the end datetimes
dts=datetime.datetime(y1,m1,1)
dte=datetime.datetime(y2,m2,1)
# Number of months between
nmonths=(dte.year - dts.year) * 12 + dte.month - dts.month +1
# Create the dates range
rdates=[ add_months(dts,m) for m in range(nmonths) ]

## Open the file where we are going to print the results.
fr=open('r'+str(method)+'_'+str(m1)+'-'+str(y1)+'_'+str(m2)+'-'+str(y2)+'_'+str(gm)+'.txt', 'w')

prints('Number of months: '+str(nmonths)+'\n')
prints('Taking the data from murders.csv...');

## Take murders information
df=pd.read_csv('murder.csv', parse_dates=['FIRST_OCCU'])
df=df.set_index(['FIRST_OCCU']); # set the index to be the dates

# Calculates and prints the elapsed time
t_final = time()
prints('Elapsed time: '+str(t_final-t_init)+'s\n')

## Generate data list "imonths" where there are all the data in the range dates.
imonths=[df.loc['%d-%d' % (d.year,d.month)] for d in rdates]
# Square approx coordinates of Denver
lat1=39.605200
lat2=39.929262
lon1=-105.105060
lon2=-104.581836
# Remove spurious data (not in the coordinates of Denver)
for index,m in enumerate(imonths):
    imonths[index]=m[(m['GEO_LON']>lon1) & (m['GEO_LON']<lon2) & (m['GEO_LAT']>lat1) & (m['GEO_LAT']<lat2)]
# Agruping the months for a number of groups specified in gm
dfs=[]
months=[]
for i,m in enumerate(imonths):
    dfs.append(m)
    if i%gm == gm-1: # This is the last element of each group
        months.append(pd.concat(dfs)) # Concatenates all the data frames in dfs
        dfs=[]

prints('Generating the variables for the model...')
prints('Number of the groups: '+str(len(months)))

## Generate xik in the model from the data
# Taking the longitudes
dlon=[[months[i]['GEO_LON'][k] for k in range(len(months[i]['GEO_LON']))] for i in range(len(months))]
# Taking the latitudes
dlat=[[months[i]['GEO_LAT'][k] for k in range(len(months[i]['GEO_LAT']))] for i in range(len(months))]
# Generating xik from latitudes and longitudes
xik=[]
for n in range(len(dlon)):
    xik.append(np.column_stack((dlon[n],dlat[n])))

prints('Number of murders: '+str(sum([len(xik[i]) for i in range(len(xik))])))

## Generate xj in the model from the data... 
if method=='2aprox':
    # ...using 2-aproxim implementation.
    prints('With 2-aprox!')
    # xj are directly the coordinates of xik
    xj=[x for sublist in xik for x in sublist]
else:
    # ...using all the combinations.
    prints('Exact solution!')
    variables=['GEO_LON','GEO_LAT']
    result=[[] for i in range(len(variables)+2)]
    # xj are all the combinations from xik and saved in result array
    combinations(months,0,[None]*(len(months)),result,variables)
    xj=np.column_stack((result[0],result[1]))
    co=result[2];
    xc=result[3];
    lambdas=np.array([(len(m)+0.)/sum([len(m) for m in months]) for m in months])
    ch=[np.dot(lambdas,(x[0]-xj[i][0])**2+(x[1]-xj[i][1])**2) for i,x in enumerate(xc)]

## Generate cikj that is the distances between xj and xik for all the possibilities    
cikj=[]
cdistance(xik,xj,cikj)

# result[2] is the indices of the xik for each xj
NN=len(co);
MM=np.array(co[-1])+1;

coeff=[];
for i in range(len(MM)):
    coeff.append(np.zeros((MM[i],NN)))
    for j in range(len(co)):
        coeff[i][int(co[j][i])][j]=1

# Calculates and prints the elapsed time
t_new = time()
prints('Elapsed time: '+str(t_new-t_final)+'s\n')

prints('Creating barycenters.dat...')

## Create barycenters.dat filling all the data
# Generate some more data needed
# Number of months (N)
N=len(months) 
# Number of possible police location (SO)
SO=len(xj) 
# Number of murders on each month labeling each month (PI)
l=[]
for i in range(1,N+1):
    l.append(i)
    l.append(len(months[i-1]))
PI=' '.join(map(str,l)) 
# Range between 1 and the max number of murders in one month (D)
ma=max([len(mon) for mon in months])
D=' '.join(map(str,range(1,ma+1)))
# Weight of each murder in each month (d)
td=[]
wmurd=[]
for i in range(N):
    td.append([i+1])
    for j in range(len(months[i])):
        td[i].append(1./len(months[i]))
        wmurd.append(1./len(months[i]));
    for j in range(len(months[i]),ma):
        td[i].append('.')
print(td)
d='\n'.join(map(str,[' '.join(map(str,k)) for k in td]))
# Weight of each month: # murders on a month divided by total # of murders in all the months (lambda: 'lam')
l=[]
for i in range(1,N+1):
    l.append(i)
    l.append((len(months[i-1])+0.)/sum([len(m) for m in months]))   
lam=' '.join(map(str,l))
# Range between 1 and the number of possible police location (so)
so=' '.join(map(str,range(1,SO+1)))
# New ch instead of c
cha=[]
for i in range(1,len(ch)+1):
    cha.append(i)
    cha.append(ch[i-1])  
chs=' '.join(map(str,cha))

prints('Lambdas: '+str(lam))

# Write all the information into barycenters.dat
f=open('barycenters.dat', 'w')
f.write('param N:= '+str(N)+';\n')
f.write('param SO:= '+str(SO)+';\n')
f.write('param PI:= '+str(PI)+';\n')
f.write('param d: '+str(D)+' :=\n')
f.write(str(d)+' ;\n')
f.write('param lambda := '+str(lam)+' ;\n')
f.write('param c\n')
for k in range(ma):
    f.write(': '+str(so)+'\n')
    f.write(': '+str(' '.join(map(str,[k+1] * SO)))+' :=\n')
    nn=np.where([k in ele for ele in [range(x) for x in [len(cikj[l]) for l in range(len(cikj))]]])[0]
    for i in nn:
        f.write(str(i+1)+' '+str(' '.join(map(str,cikj[i][k])))+'\n')
f.write(';\n')
f.write('param coef\n')
for k in range(ma):
    f.write(': '+str(so)+'\n')
    f.write(': '+str(' '.join(map(str,[k+1] * SO)))+' :=\n')
    nn=np.where([k in ele for ele in [range(x) for x in [len(coeff[l]) for l in range(len(coeff))]]])[0]
    for i in nn:
        f.write(str(i+1)+' '+str(' '.join(map(str,coeff[i][k])))+'\n')
f.write(';\n')
f.write('param ch := '+str(chs)+' ;\n')
f.close()

# Calculates and prints the elapsed time
t_final = time()
prints('Elapsed time: '+str(t_final-t_new)+'s\n')


#####################################
#xiklon=[[x[0] for x in sublist ] for sublist in xik]
#xiklat=[[x[1] for x in sublist ] for sublist in xik]

#with open('gmap.pickle', 'rb') as f: 
#    gmap = pickle.load(f)
#for i in range(len(xiklat)):
#    marker=i % 12 + 1 # number of the marker, from 1 to 12
#    gmap.scatter(xiklat[i],xiklon[i],'#month%d' % marker, marker=True) # scatter plot over the map
#for i in range(len(xiklat)):
#    gmap.scatter(xiklat[i],xiklon[i],'#month8', marker=True)
#gmap.draw('map.html')
## Change the path of the markers to a local path (pointing into a python folder, sometimes does not find them)
#with open('map.html', "r") as sources:
#    lines = sources.readlines()
#with open('map.html', "w") as sources:
#    for line in lines:
#        sources.write(re.sub('C:\\\\Anaconda2\\\\lib\\\\site-packages\\\\gmplot\\\\markers','./markers',line))

#####################################

#import matplotlib.pyplot as plt
#from sklearn import cluster, datasets, mixture
#from sklearn.neighbors import kneighbors_graph
#from sklearn.preprocessing import StandardScaler
#from itertools import cycle, islice

#X=[x for sublist in xik for x in sublist]
#Xlon=[x[0] for sublist in xik for x in sublist]
#Xlat=[x[1] for sublist in xik for x in sublist]

    # normalize dataset for easier parameter selection
    #   X = StandardScaler().fit_transform(X)

#spectral = cluster.SpectralClustering(
#        n_clusters=8, eigen_solver='arpack',
#        affinity="nearest_neighbors")

#spectral.fit(X)

#y = spectral.labels_.astype(np.int)

#print(y)

#colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a',
#                                             '#f781bf', '#a65628', '#984ea3',
#                                             '#999999', '#e41a1c', '#dede00']),
#                                      int(max(y) + 1))))

#plt.scatter(Xlon,Xlat, s=10, color=colors[y])
#plt.show()

#sys.exit()


prints('Solving the model using AMPL...')
## Execute AMPL through barycenters.run using barycenters.mod and the new barycenters.dat
command="ampl barycenters2.run"
process=subprocess.Popen(command.split(), stdout=subprocess.PIPE)
out=process.communicate()

# Print the output from AMPL
prints(str(out))

# Calculates and prints the elapsed time
t_new = time()
prints('Elapsed time: '+str(t_new-t_final)+'s\n')

prints('Taking the outputs from AMPL...')

## Use AMPL output
# Take data from AMPL outputs (transport.csv and output.csv)
if method=='2aprox':
    # read the file transport.csv
    df2=pd.read_csv('transport.csv',header=None)
    print(df2)
    # compute all the |Pi|
    PIv=[ len(m) for m in months ]
    trans=[]
    xjik=[]
    count=0
    # save all the transport variable in trans array
    for i in range(N):
        trans.append([])
        for j in range(SO):
            trans[i].append([])
            for k in range(PIv[i]):
                trans[i][j].append(df2[0][count])
                if trans[i][j][k] > tolerance: # select only when larger than tolerance
                    # create a list with all the xik for each xj with zj larger than tolerance
                    lis=np.array([ x[0] for x in xjik ])
                    if not any(lis==j):
                        xjik.append([j,[xik[i][k]]])
                    else:
                        xjik[int(list(lis).index(j))][1].append(xik[i][k])
                count+=1
    # save the final xj in result variable
    result=[ x[1] for x in xjik ]
else:  
    # read 'output.csv' into data frame 'df'
    df=pd.read_csv('output.csv',header=None)
    # indexes where larger than tolerance (we discard the others)
    index=df[df[1]>tolerance]
    # the previous command keeps the original index numbers, then reset the index numbers from 0
    index=index.reset_index()
    # take the final xj using the indexes from the first column of the csv (index[0])
    result=[ xj[int(x)] for x in np.array([index[0][i] for i in range(len(index[0]))])-1 ]
    # saving the weights of each xj (zj) which are in the second column of the csv (index[1])
    weights=np.array([index[1][i] for i in range(len(index[1]))])
    prints(str(weights))
    prints('Sum of the weights: '+str(sum(weights)))

# Calculates and prints the elapsed time
t_final=time()
prints('Elapsed time: '+str(t_final-t_new)+'s\n')

prints('Generating map.html...')

## Open Denver map
# - Save a map of Denver in 'gmap'
#gmap = gmplot.GoogleMapPlotter.from_geocode("Denver");
# - Save 'gmap' map in the file 'gmap.pickle'
#with open('gmap.pickle', 'wb') as f:  
#    pickle.dump(gmap, f)
# - Getting back the map in 'gmap.pickle' into the variable 'gmap'
with open('gmap.pickle', 'rb') as f: 
    gmap = pickle.load(f)

## Take lons and lats of the xik and result
xiklon=[[x[0] for x in sublist ] for sublist in xik]
xiklat=[[x[1] for x in sublist ] for sublist in xik]
if method=='2aprox':
    lons=[[ x[0] for x in sub ] for sub in result]
    xjlon=[np.mean(x) for x in lons]
    lats=[[ x[1] for x in sub ] for sub in result]
    xjlat=[np.mean(x) for x in lats]
else:
    xjlon=[x[0] for x in result]
    xjlat=[x[1] for x in result]

## Calculate the new objective function for 2aprox method (after changing the position)
if method=='2aprox':
    index=[ x[0] for x in xjik ]
    xjn=xj
    for ix,ind in enumerate(index):
        xjn[ind]=[xjlon[ix],xjlat[ix]]
    cresult=[]
    cdistance(xik,xjn,cresult)
    # computing the objective function before changing the position of the police
    objfunc=0
    for i in range(N):
        suma=0
        # labdas
        lamb=(len(months[i])+0.)/sum([len(m) for m in months])
        for j in range(len(xj)):
            for k in range(len(months[i])):
                # sum all the distances by transport variable
                suma+=cikj[i][k][j]*trans[i][j][k]
        # sum over all the labmda by the previous sum
        objfunc+=lamb*suma
    # computing the objective function after changing the position of the police
    objfunc2=0
    for i in range(N):
        suma=0
        # labdas
        lamb=(len(months[i])+0.)/sum([len(m) for m in months])
        for j in range(len(xjn)):
            for k in range(len(months[i])):
                # sum all the distances by transport variable
                suma+=cresult[i][k][j]*trans[i][j][k]
        # sum over all the labmda by the previous sum
        objfunc2+=lamb*suma

    prints('Objective function without modification: '+str(objfunc))
    prints('Objective function with modification: '+str(objfunc2))

## Print xik and result in the map
# - print all the xik using month classification by colors
for i in range(len(xiklat)):
    marker=i % 12 + 1 # number of the marker, from 1 to 12
    gmap.scatter(xiklat[i],xiklon[i],'#month%d' % marker, marker=True) # scatter plot over the map
# - print all the xik using a red marker
#for i in range(len(xiklat)):
#    gmap.scatter(xiklat[i],xiklon[i],'#month8', marker=True)
# - print all the resulting xj using a blue marker
gmap.scatter(xjlat,xjlon,'#police', marker=True)

xiklons=[x for sub in xiklon for x in sub]
xiklats=[x for sub in xiklat for x in sub]

murders = pd.DataFrame(
    {'Lon': xiklons,
     'Lat': xiklats,
     'Weights': wmurd
    })

police = pd.DataFrame(
    {'Lon': xjlon,
     'Lat': xjlat,
     'Weights': weights
    })

murders.to_csv('out_murders.csv');
police.to_csv('out_police.csv');

## Save the map
gmap.draw(str(method)+'map_'+str(m1)+'-'+str(y1)+'_'+str(m2)+'-'+str(y2)+'_'+str(gm)+'.html')

## Change the path of the markers to a local path (pointing into a python folder, sometimes does not find them)
with open(str(method)+'map_'+str(m1)+'-'+str(y1)+'_'+str(m2)+'-'+str(y2)+'_'+str(gm)+'.html', "r") as sources:
    lines = sources.readlines()
with open(str(method)+'map_'+str(m1)+'-'+str(y1)+'_'+str(m2)+'-'+str(y2)+'_'+str(gm)+'.html', "w") as sources:
    for line in lines:
        sources.write(re.sub('C:\\\\Anaconda2\\\\lib\\\\site-packages\\\\gmplot\\\\markers','./markers',line))

# Calculates and prints the elapsed time
t_new=time()
prints('Elapsed time: '+str(t_new-t_final)+'s\n')

# Calculates and prints the total time
prints('Total time: '+str(t_new-t_init)+'s')
