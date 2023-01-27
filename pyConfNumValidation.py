"""

Code to validate an OpenFOAM model of the
laboratory confluence

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

from sklearn.linear_model import LinearRegression

import pivToolbox.pivToolbox as piv
import pyFoamTools.pyFoamTools as foam

#rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':'8'})
rc('font',**{'family':'serif','serif':['Times New Roman'],'size':'7'})
rc('text', usetex=True)


"""
-------------------------------------------------------------------------------
--------------------Section #1 - Data import and treatment---------------------
-------------------------------------------------------------------------------
"""



"""
Step #1
Specify paths, case parameters, etc. 
"""


inst_PIV_path=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\Instantaneous\Case_8\Case_8_0.06s_inst_data"
mean_PIV_path=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\PlanarPIV_at_4D\Case_8_"

openFOAM_path=r'X:\LabConfluence_results\Hi2\lab1mm_pp'
cfdMeanTime='161'
cfdProbeTime='126'
figureOutputPath=openFOAM_path+r'\postProcessing\Figures'

try:
    os.makedirs(figureOutputPath)
    print('Created Figure folder.')
         
except FileExistsError:
    print('Figure folder already exists. Not overwriting.')
    pass



"""
Step #2

Import instantaneous PIV data output from DaVis (i.e. multiple files for each time)

CAUTION: DaVis always exports 'u' as the left-right component in the field of view, and 'v' the bottom-up component)
At the 4D plane, 'u' from DaVis is actually the lateral (v) component and 'v' from DaVis is actually the vertical (w) component
remember, 'negative' lateral velocities are towards the left

"""

vInstPIV_paths=inst_PIV_path+"\\Inst_u"
vInstPIV_files=piv.get_PIV_files(vInstPIV_paths)

wInstPIV_paths=inst_PIV_path+"\\Inst_v"
wInstPIV_files=piv.get_PIV_files(wInstPIV_paths)

label=['Case_8_']
times=[]

for file in vInstPIV_files:
    times.append(os.path.splitext(os.path.basename(file))[0])

instPIVData={}
instPIVData['v']=piv.load_instPIV_data(vInstPIV_files,times,'v',ref=True,shift=[0,25.46,0]) 
instPIVData['w']=piv.load_instPIV_data(wInstPIV_files,times,'w',ref=True,shift=[0,25.46,0]) 


"""
Step #3

Import mean PIV data output from DaVis (i.e. a single files for field exported for each experiment)
"""


"""Provide data labels and load data to pandas dataframes"""
meanPIV_fieldKeys=['B00001','B00002','B00005','B00006','B00007','B00008'] #u,v,TKE
meanPIV_fieldNames=['u','v','tke','XY','XX','YY']
meanPIV_limitFolderTo='Smoothed_9x9'
meanPIV_shotLabels=['Shot_1','Shot_2']
meanPIV_planeLabels=['Case_8_']

meanPIV_files=piv.get_PIV_files(mean_PIV_path,limitFolderTo=meanPIV_limitFolderTo)

meanPIV_data=piv.load_PIV_dataframes(meanPIV_files,meanPIV_planeLabels,meanPIV_shotLabels,meanPIV_fieldNames,meanPIV_fieldKeys)

meanPIV_treatedData={}
meanPIV_treatedData['Vector']=piv.treatMeanVecData(meanPIV_data,[0,25.46,0],[-110,110,0,80],2.5,avg=True,norm=False)
meanPIV_treatedData['Scalar']=piv.treatScalarData(meanPIV_data,['tke','XY','XX','YY'],[0,25.46,0],[-110,110,0,80],2.5,avg=True,norm=False)



meanPIV_fullFOV=meanPIV_treatedData['Vector']['Case_8_']['avg']


"""
Step #4 - Make a grid of points to sample on 4D plane
"""

vArrays=[]
wArrays=[]

for time in times:
    vArrays.append(instPIVData['v'][time]['v'])
    wArrays.append(instPIVData['w'][time]['w'])

M, N = 12, 5
x = np.linspace(-90,90,M+1)
y = np.linspace(10,60,N+1)
X,Y = np.meshgrid(x,y)

pivPoints = np.vstack([X.ravel(), Y.ravel()]).transpose()


"""
Step #5 - Extract mean data from instantaneous PIV data on grid 
"""

vMean=[]
wMean=[]

for count, point in enumerate(pivPoints):
    print('Extracting mean data at point %s' %point)
    vMean.append(piv.extractInstPointData(vArrays,point).mean())
    wMean.append(piv.extractInstPointData(wArrays,point).mean())

mean_vecPIV=pd.DataFrame(pivPoints,columns=['y','z'])
mean_vecPIV['v']=vMean
mean_vecPIV['w']=wMean


"""
Step #6 - Extract mean CFD data on grid points

Extract data from a single OpenFOAM timestep (i.e. for mean values or turbulent stats)

-creates set of points to sample from OpenFOAM data
-creates a 'probes' file for use with postProcess function
-use the probes file with: mpirun -np XX postProcess -func probes -parallel -time 'XX'

"""

#write a probe point file for OpenFOAM to extract data with
M, N = 12, 5
y = np.linspace(0.09,-0.09,M+1)
z = np.linspace(0.01,0.06,N+1)
Y,Z = np.meshgrid(y,z)
cfdPoints = np.vstack([Y.ravel(), Z.ravel()]).transpose()

pointsToExtract=pd.DataFrame(cfdPoints,columns=['y','z'])
pointsToExtract['x']=0.28 #fixed value of x at 4D = 28 cm
pointsToExtract = pointsToExtract[['x', 'y', 'z']]

foam.writeProbesFileFromArray(np.asarray(pointsToExtract),'UMean',openFOAM_path+'\system',"probes")

#Now run: mpirun -np XX postProcess -func probes -parallel -time 'XX' from case directory
#this will create a file 'openFOAM_path/postProcessing/probes/0/UMean' containing the sampled data for use below

"""
Step #7 - Import mean probe data from OpenFOAM
"""

#create list of point names
ptNames=[]
for count, point in enumerate(cfdPoints):
    ptNames.append('P%d' %count)
    
cfd_probePath=openFOAM_path+r'\postProcessing\probes\%s\UMean' %cfdMeanTime
treatedPath=openFOAM_path+r'\postProcessing\probes\%s\treated\U' %cfdMeanTime
try:
    os.makedirs(openFOAM_path+r'\postProcessing\probes\%s\treated' %cfdMeanTime)
    standinFile = open(treatedPath, "w")
    standinFile.close()
         
except FileExistsError:
    standinFile = open(treatedPath, "w")
    standinFile.close()
    pass

mean_vecCFD=foam.readProbeFile(cfd_probePath, treatedPath, 'vector', ['u','v','w'], ptNames).pointData
mean_vecCFD=pd.concat(mean_vecCFD).reset_index()
mean_vecCFD=mean_vecCFD.drop(['level_0', 'level_1','Time'], axis=1)    # clean up
mean_vecCFD = pd.concat([mean_vecCFD, pointsToExtract], axis=1)        # add x, y, z to dataframe

mean_vecCFD['v']=mean_vecCFD['v']*-1                                   #to make leftwards flow negative, like in the PIV!!
mean_vecCFD['y']=mean_vecCFD['y']*-1                                   #to invert lateral axis


"""
Step #8 - Import instantaneous probe data from OpenFOAM
"""

#instantaneous data collected during runtime
cfd_probePath=openFOAM_path+r'\postProcessing\probesU\%s\U' % cfdProbeTime 
treatedPath=openFOAM_path+r'\postProcessing\probesU\%s\treated\U' % cfdProbeTime

try:
    os.makedirs(openFOAM_path+r'\postProcessing\probesU\%s\treated'  % cfdProbeTime )
    standinFile = open(treatedPath, "w")
    standinFile.close()
except FileExistsError:
    standinFile = open(treatedPath, "w")
    standinFile.close()
    pass

pointNames=['p0','p1','p2','p3','p4','p5']
cfd_instProbeData=foam.readProbeFile(cfd_probePath, treatedPath, 'vector', ['u','v','w'], pointNames)


"""
Step # 9: Extract and import line data from OpenFOAM case
"""

def writeSingleGraphOpenFOAMfile(fileName,directory,start,end,fields,nPoints):
        
    f= open(directory+'//'+fileName,"w+")
        
    f.write("start (%s);\n" % start)
    f.write("end (%s);\n" % end)
    f.write("fields (%s);\n" %fields)

    f.write("interpolationScheme cellPoint;\n")
    f.write("setFormat raw;\n")
    
    f.write("setConfig\n")
    f.write("{\n")
    f.write("type uniform;\n")
    f.write("axis distance;\n")
    f.write("nPoints %d;\n" %nPoints)
    f.write("}\n")

    f.write('#includeEtc "caseDicts/postProcessing/graphs/graph.cfg"\n')
    f.close()  
    print('Successfully wrote array to file.')

def writeSingleGraphShellScript(outputDirectory,outputFileName,nbrProcs,lineNames,time):
    
    f= open(outputDirectory+'//'+outputFileName,"w+")
    f.write("#!/bin/bash\n")
    for name in lineNames:
        f.write("mpirun -np %d postProcess -func %s -parallel -time '%s'\n" % (nbrProcs, name, time))
    f.close()  
    print('Successfully wrote bash shell script for postProcess.')



lineNames=['line_0','line_neg035','line_neg07','line_pos035','line_pos07']
lineLocations=[['0.28 0.0 0.0','0.28 0.0 0.07'],
               ['0.28 0.035 0.0','0.28 0.035 0.07'],
               ['0.28 0.07 0.0','0.28 0.07 0.07'],
               ['0.28 -0.035 0.0','0.28 -0.035 0.07'],
               ['0.28 -0.07 0.0','0.28 -0.07 0.07']
               ]

for count, line in enumerate(lineNames):
    print(line)
    nmbPoints= 14
    fields='alpha.denseMean UMean turbulenceProperties.RMean UPrime2Mean'
    writeSingleGraphOpenFOAMfile(line,openFOAM_path+'\system',lineLocations[count][0],lineLocations[count][1],fields,nmbPoints)

writeSingleGraphShellScript(openFOAM_path, 'postProcessLines.sh', 64, lineNames, cfdMeanTime)
#Afterwards, run the shell script from the OpenFOAM case file. 
#You might need to run dos2unix on the file to overcome a ^M bad interpreter error


#load line data

cfdLineData={}
for line in lineNames:
    file=openFOAM_path+r'\postProcessing\%s\%s\%s.xy' %(line,cfdMeanTime,line)
    cfdLineData[line]=pd.read_table(file,delim_whitespace=True,names=['z','uMean','vMean','wMean','resXX','resXY','resXZ','resYY','resYZ','resZZ', 'sgsXX','sgsXY','sgsXZ','sgsYY','sgsYZ','sgsZZ'])
    cfdLineData[line]['vMean']=cfdLineData[line]['vMean']*-1 #to match coordinate orientation of PIV (i.e. negative lateral flow to the left)
    cfdLineData[line]['tke2C']=0.5*(cfdLineData[line]['resYY']+cfdLineData[line]['resZZ'] + cfdLineData[line]['sgsYY']+cfdLineData[line]['sgsZZ'])
    cfdLineData[line]['tke3C']=0.5*(cfdLineData[line]['resYY']+cfdLineData[line]['resYY']+cfdLineData[line]['resZZ'] + cfdLineData[line]['sgsXX']+cfdLineData[line]['sgsYY']+cfdLineData[line]['sgsZZ'])





"""
-------------------------------------------------------------------------------
----------------------Section #2 - Figure production---------------------------
-------------------------------------------------------------------------------
"""




"""
Grid plots of v and w between CFD and PIV
"""

#Lateral component
fig, ax = plt.subplots(2,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)
im=ax[0].scatter(mean_vecPIV['y']*0.001,mean_vecPIV['z']*0.001, c=mean_vecPIV['v'],cmap='seismic')
im=ax[1].scatter(mean_vecCFD['y'], mean_vecCFD['z'],c=mean_vecCFD['v'],cmap='seismic')
ax[0].set_xlim([-0.1,0.1])   
ax[0].set_ylim([0.0,0.07]) 
ax[1].set_xlim([-0.1,0.1])   
ax[1].set_ylim([0.0,0.07])                   
ax[0].set_xlabel('$y$ (m)', fontsize=10)
ax[1].set_xlabel('$y$ (m)', fontsize=10)
ax[0].set_ylabel('$v$ (m)', fontsize=10)
ax[1].set_ylabel('$v$ (m)', fontsize=10)
ax[0].text(0.02,0.9,"a)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[0].transAxes)
ax[1].text(0.02,0.9,"b)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[1].transAxes)
ax[0].text(0.05,-0.25,"$\mathbf{PIV}$", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[0].transAxes)
ax[1].text(0.05,-0.25,"$\mathbf{CFD}$", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[1].transAxes)
cbar=fig.colorbar(im, ax=ax.ravel().tolist(),shrink=0.5)
cbar.set_label('$\overline{v}$ (m/s)', rotation=90)
plt.savefig(figureOutputPath+'\Gridv_compare.pdf', format='pdf', dpi=600)


#Vertical component
fig, ax = plt.subplots(2,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)
im=ax[0].scatter(mean_vecPIV['y']*0.001,mean_vecPIV['z']*0.001, c=mean_vecPIV['w'],cmap='seismic')
im=ax[1].scatter(mean_vecCFD['y'], mean_vecCFD['z'],c=mean_vecCFD['w'],cmap='seismic')
ax[0].set_xlim([-0.1,0.1])   
ax[0].set_ylim([0.0,0.07]) 
ax[1].set_xlim([-0.1,0.1])   
ax[1].set_ylim([0.0,0.07])                   
ax[0].set_xlabel('$y$ (m)', fontsize=10)
ax[1].set_xlabel('$y$ (m)', fontsize=10)
ax[0].set_ylabel('$w$ (m)', fontsize=10)
ax[1].set_ylabel('$w$ (m)', fontsize=10)
ax[0].text(0.02,0.9,"a)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[0].transAxes)
ax[1].text(0.02,0.9,"b)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[1].transAxes)
ax[0].text(0.05,-0.25,"$\mathbf{PIV}$", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[0].transAxes)
ax[1].text(0.05,-0.25,"$\mathbf{CFD}$", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[1].transAxes)
cbar=fig.colorbar(im, ax=ax.ravel().tolist(),shrink=0.5)
cbar.set_label('$\overline{w}$ (m/s)', rotation=90)
plt.savefig(figureOutputPath+'\Gridw_compare.pdf', format='pdf', dpi=600)


"""
Compare mean PIV and CFD - RMSE, scatter & linear regression
"""


avg_v_CFD=mean_vecCFD['v'].mean()
avg_v_PIV=mean_vecPIV['v'].mean()

avg_w_CFD=mean_vecCFD['w'].mean()
avg_w_PIV=mean_vecPIV['w'].mean()

#RMSE of lateral and vertical components
RMSE_v=np.sqrt(((mean_vecCFD['v']-mean_vecPIV['v'])**2).mean())
print('RMSE_v is %.2f%% of mean v' %(RMSE_v/abs(avg_v_PIV)*100))
RMSE_w=np.sqrt(((mean_vecCFD['w']-mean_vecPIV['w'])**2).mean())
print('RMSE_w is %.2f%% of mean w' %(RMSE_w/abs(avg_w_PIV)*100))

#regression of v
modelv=LinearRegression().fit(np.asarray(mean_vecPIV['v']).reshape((-1, 1)),mean_vecCFD['v'])
scorev=modelv.score(np.asarray(mean_vecPIV['v']).reshape((-1, 1)),mean_vecCFD['v'])
print('The coefficient of determination for v is: %f' %scorev)
print(f"intercept: {modelv.intercept_}")
print(f"slope: {modelv.coef_}")

#regression of w
modelw=LinearRegression().fit(np.asarray(mean_vecPIV['w']).reshape((-1, 1)),mean_vecCFD['w'])
scorew=modelw.score(np.asarray(mean_vecPIV['w']).reshape((-1, 1)),mean_vecCFD['w'])
print('The coefficient of determination for w is: %f' %scorew)
print(f"intercept: {modelw.intercept_}")
print(f"slope: {modelw.coef_}")

#plot for article
fig, ax = plt.subplots(1,2,frameon=False,constrained_layout=True)
fig.set_figwidth(5.5,forward=True)
fig.set_figheight(2.4,forward=True)

ax[0].plot(mean_vecPIV['v'], modelv.coef_*mean_vecPIV['v']+modelv.intercept_,zorder=0,c='gray')
ax[1].plot(mean_vecPIV['w'], modelw.coef_*mean_vecPIV['w']+modelw.intercept_,zorder=0,c='gray')

ax[0].scatter(mean_vecPIV['v'],mean_vecCFD['v'],s=20,facecolor='none',edgecolor='black',zorder=1)
ax[1].scatter(mean_vecPIV['w'],mean_vecCFD['w'],s=20,facecolor='none',edgecolor='black',zorder=1)
ax[0].set_xlim([-0.05,0.05])   
ax[0].set_ylim([-0.05,0.05])                   
ax[1].set_xlim([-0.03,0.03])  
ax[1].set_ylim([-0.03,0.03])  

ax[0].set_xlabel('$v_{piv}$ (m/s)', fontsize=10)
ax[1].set_xlabel('$w_{piv}$ (m/s)', fontsize=10)
ax[0].set_ylabel('$v_{cfd}$ (m/s)', fontsize=10)
ax[1].set_ylabel('$w_{cfd}$ (m/s)', fontsize=10)

ax[0].text(0.05,0.95,"a)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[0].transAxes)
ax[0].text(0.02,-0.02,"y = %.2fx + %.2e" %(modelv.coef_,modelv.intercept_), fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')
ax[0].text(0.02,-0.026,"R$^{2}$ = %.2f" %(scorev), fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')

ax[1].text(0.05,0.95,"b)", fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center',transform=ax[1].transAxes)
ax[1].text(0.01,-0.02,"y = %.2fx + %.2e" %(modelw.coef_,modelw.intercept_), fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')
ax[1].text(0.01,-0.023,"R$^{2}$ = %.2f" %(scorew), fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')

plt.savefig(figureOutputPath+'\Scatter_compare.pdf', format='pdf', dpi=600)



"""
Line plot compairisons between CFD and PIV
-Use: "mpirun -np 64 postProcess -func singleGraph -parallel -time 'xx'" to extract UMean on a line
-Repeat for the 5 lines (i.e. y positions of -0.07, -0.035, 0, 0.035, 0.07 )
"""
#requires postProcess line data from OpenFOAM case! 


#extract lines from mean PIV flow field data
linePIV_neg07=piv.extractProfile(meanPIV_fullFOV['db'], 'x', 'y', -70)
linePIV_neg035=piv.extractProfile(meanPIV_fullFOV['db'], 'x', 'y', -35)
linePIV_0=piv.extractProfile(meanPIV_fullFOV['db'], 'x', 'y', 0)
linePIV_pos035=piv.extractProfile(meanPIV_fullFOV['db'], 'x', 'y', 35)
linePIV_pos07=piv.extractProfile(meanPIV_fullFOV['db'], 'x', 'y', 70)

pivLines=[linePIV_neg07,linePIV_neg035,linePIV_0,linePIV_pos035,linePIV_pos07]

linePIV_neg07.columns=['y','z','v','w']
linePIV_neg035.columns=['y','z','v','w']
linePIV_0.columns=['y','z','v','w']
linePIV_pos035.columns=['y','z','v','w']
linePIV_pos07.columns=['y','z','v','w']


#Create plots for publication

#lateral components
fig, ax = plt.subplots(1,5,frameon=False,constrained_layout=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(3.0,forward=True)

ax[0].scatter(cfdLineData['line_neg07']['vMean'],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='b')
ax[1].scatter(cfdLineData['line_neg035']['vMean'],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='b')
ax[2].scatter(cfdLineData['line_0']['vMean'],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='b')
ax[3].scatter(cfdLineData['line_pos035']['vMean'],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='b')
lineCFD=ax[4].scatter(cfdLineData['line_pos07']['vMean'],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='b',label='CFD')

ax[0].scatter(linePIV_neg07['v'],linePIV_neg07['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[1].scatter(linePIV_neg035['v'],linePIV_neg035['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[2].scatter(linePIV_0['v'],linePIV_0['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[3].scatter(linePIV_pos035['v'],linePIV_pos035['z']*0.001,s=20,facecolor='none',edgecolor='black')
linePIV=ax[4].scatter(linePIV_pos07['v'],linePIV_pos07['z']*0.001,s=20,facecolor='none',edgecolor='black',label='PIV')

ax[4].legend(handles=[linePIV, lineCFD])

for i in range(len(ax)):
    ax[i].set_xlim([-0.04,0.04])   
    ax[i].set_ylim([0,0.07])            
         
ax[0].title.set_text('$y$=-0.070 m')
ax[1].title.set_text('$y$=-0.035 m')
ax[2].title.set_text('$y$=-0.000 m')
ax[3].title.set_text('$y$= 0.035 m')
ax[4].title.set_text('$y$= 0.070 m')

ax[0].set_xlabel('$w$ (m/s)', fontsize=10)
ax[1].set_xlabel('$w$ (m/s)', fontsize=10)
ax[2].set_xlabel('$w$ (m/s)', fontsize=10)
ax[3].set_xlabel('$w$ (m/s)', fontsize=10)
ax[4].set_xlabel('$w$ (m/s)', fontsize=10)
ax[0].set_ylabel('$z$ (m)', fontsize=10)

plt.savefig(figureOutputPath+'\Profiles_vcompare.pdf', format='pdf', dpi=600)


#vertical components
fig, ax = plt.subplots(1,5,frameon=False,constrained_layout=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(3.0,forward=True)

ax[0].scatter(cfdLineData['line_neg07']['wMean'],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='b')
ax[1].scatter(cfdLineData['line_neg035']['wMean'],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='b')
ax[2].scatter(cfdLineData['line_0']['wMean'],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='b')
ax[3].scatter(cfdLineData['line_pos035']['wMean'],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='b')
lineCFD=ax[4].scatter(cfdLineData['line_pos07']['wMean'],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='b',label='CFD')

ax[0].scatter(linePIV_neg07['w'],linePIV_neg07['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[1].scatter(linePIV_neg035['w'],linePIV_neg035['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[2].scatter(linePIV_0['w'],linePIV_0['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[3].scatter(linePIV_pos035['w'],linePIV_pos035['z']*0.001,s=20,facecolor='none',edgecolor='black')
linePIV=ax[4].scatter(linePIV_pos07['w'],linePIV_pos07['z']*0.001,s=20,facecolor='none',edgecolor='black',label='PIV')

ax[4].legend(handles=[linePIV, lineCFD])


for i in range(len(ax)):
    
    ax[i].set_xlim([-0.015,0.015])   
    ax[i].set_ylim([0,0.07])            
         
ax[0].title.set_text('$y$=-0.070 m')
ax[1].title.set_text('$y$=-0.035 m')
ax[2].title.set_text('$y$=-0.000 m')
ax[3].title.set_text('$y$= 0.035 m')
ax[4].title.set_text('$y$= 0.070 m')

ax[0].set_xlabel('$w$ (m/s)', fontsize=10)
ax[1].set_xlabel('$w$ (m/s)', fontsize=10)
ax[2].set_xlabel('$w$ (m/s)', fontsize=10)
ax[3].set_xlabel('$w$ (m/s)', fontsize=10)
ax[4].set_xlabel('$w$ (m/s)', fontsize=10)
ax[0].set_ylabel('$z$ (m)', fontsize=10)

plt.savefig(figureOutputPath+'\Profiles_wcompare.pdf', format='pdf', dpi=600)




"""
Compare time series of instantaneous velocities
"""

pivProbe0=piv.extractInstPointData(vArrays, [0, 10])
pivProbe1=piv.extractInstPointData(vArrays, [0, 35])
pivProbe2=piv.extractInstPointData(vArrays, [0, 60])
pivProbe3=piv.extractInstPointData(vArrays, [-60, 10])
pivProbe4=piv.extractInstPointData(vArrays, [-60, 35])
pivProbe5=piv.extractInstPointData(vArrays, [-60, 60])


def treatProbeData(df):
    cfdProbe=np.asarray(df['v'])*-1
    # cfdProbe=cfdProbe[80000:119630]
    cfdProbe=cfdProbe[0::30]#30 is for a 0.002 timestep, it matches the 0.06 s sampling from PIV
    return cfdProbe

cfdProbe0=treatProbeData(cfd_instProbeData.pointData['p0'])
cfdProbe1=treatProbeData(cfd_instProbeData.pointData['p1'])
cfdProbe2=treatProbeData(cfd_instProbeData.pointData['p2'])
cfdProbe3=treatProbeData(cfd_instProbeData.pointData['p3'])
cfdProbe4=treatProbeData(cfd_instProbeData.pointData['p4'])
cfdProbe5=treatProbeData(cfd_instProbeData.pointData['p5'])


fig, ax = plt.subplots(6,1,frameon=False,constrained_layout=True,sharex=True,sharey=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(8,forward=True)
lineCFD0=ax[0].plot(cfdProbe0, c='black',label='CFD',lw=1)
lineCFD1=ax[1].plot(cfdProbe1, c='black',label='CFD',lw=1)
lineCFD2=ax[2].plot(cfdProbe2, c='black',label='CFD',lw=1)
lineCFD3=ax[3].plot(cfdProbe3, c='black',label='CFD',lw=1)
lineCFD4=ax[4].plot(cfdProbe4, c='black',label='CFD',lw=1)
lineCFD5=ax[5].plot(cfdProbe5, c='black',label='CFD',lw=1)
linePIV0=ax[0].plot(pivProbe0, c='blue',label='PIV',lw=1)
linePIV1=ax[1].plot(pivProbe1, c='blue',label='PIV',lw=1)
linePIV2=ax[2].plot(pivProbe2, c='blue',label='PIV',lw=1)
linePIV3=ax[3].plot(pivProbe3, c='blue',label='PIV',lw=1)
linePIV4=ax[4].plot(pivProbe4, c='blue',label='PIV',lw=1)
linePIV5=ax[5].plot(pivProbe5, c='blue',label='PIV',lw=1)

ax[5].xaxis.set_ticks([0,100,200,300,400,500,600,700])
ax[5].set_xticklabels(['0','6','12','18','24','30','36','42'])
ax[0].title.set_text('point: y=0.0, z=0.010 ')
ax[1].title.set_text('point: y=0.0, z=0.035 ')
ax[2].title.set_text('point: y=0.0, z=0.060 ')
ax[3].title.set_text('point: y=-0.06, z=0.010 ')
ax[4].title.set_text('point: y=-0.06, z=0.035 ')
ax[5].title.set_text('point: y=-0.06, z=0.060 ')

ax[5].set_xlim([0,600])  
ax[5].set_xlabel('$t$ (s)', fontsize=10)

for axes in ax:
    axes.set_ylabel('$v$ (m/s)', fontsize=10)

plt.savefig(figureOutputPath+'\TimesSeries.pdf', format='pdf', dpi=600)






"""
-------------------------------------------------------------------------------
--------------------Comparison with turbulent statistics-----------------------
-------------------------------------------------------------------------------
"""


"""
Reynolds XY compairison plot 
"""

PIVscalar='XY'
cfdScalar1='resXY'
cfdScalar2='sgsXY'

#extract lines from mean PIV flow field data
PIV=meanPIV_treatedData['Scalar']['Case_8_'][PIVscalar]['avg']['db']
linePIV_neg07=piv.extractProfile(PIV, 'x', 'y', -70)
linePIV_neg035=piv.extractProfile(PIV, 'x', 'y', -35)
linePIV_0=piv.extractProfile(PIV, 'x', 'y', 0)
linePIV_pos035=piv.extractProfile(PIV, 'x', 'y', 35)
linePIV_pos07=piv.extractProfile(PIV, 'x', 'y', 70)

# pivLines=[linePIV_neg07,linePIV_neg035,linePIV_0,linePIV_pos035,linePIV_pos07]

linePIV_neg07.columns=['y','z',PIVscalar]
linePIV_neg035.columns=['y','z',PIVscalar]
linePIV_0.columns=['y','z',PIVscalar]
linePIV_pos035.columns=['y','z',PIVscalar]
linePIV_pos07.columns=['y','z',PIVscalar]


#Create line plots for publication

#lateral components
fig, ax = plt.subplots(1,5,frameon=False,constrained_layout=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(3.0,forward=True)
ax[0].scatter(cfdLineData['line_neg07'][cfdScalar1]+cfdLineData['line_neg07'][cfdScalar2],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='b')
ax[1].scatter(cfdLineData['line_neg035'][cfdScalar1]+cfdLineData['line_neg035'][cfdScalar2],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='b')
ax[2].scatter(cfdLineData['line_0'][cfdScalar1]+cfdLineData['line_0'][cfdScalar2],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='b')
ax[3].scatter(cfdLineData['line_pos035'][cfdScalar1]+cfdLineData['line_pos035'][cfdScalar2],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='b')
lineCFD=ax[4].scatter(cfdLineData['line_pos07'][cfdScalar1]+cfdLineData['line_pos07'][cfdScalar2],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='b',label='CFD')
ax[0].scatter(linePIV_neg07[PIVscalar],linePIV_neg07['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[1].scatter(linePIV_neg035[PIVscalar],linePIV_neg035['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[2].scatter(linePIV_0[PIVscalar],linePIV_0['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[3].scatter(linePIV_pos035[PIVscalar],linePIV_pos035['z']*0.001,s=20,facecolor='none',edgecolor='black')
linePIV=ax[4].scatter(linePIV_pos07[PIVscalar],linePIV_pos07['z']*0.001,s=20,facecolor='none',edgecolor='black',label='PIV')

ax[4].legend(handles=[linePIV, lineCFD])

for i in range(len(ax)):
    #ax[i].set_xlim([-0.04,0.04])   
    ax[i].set_ylim([0,0.07])            
         

ax[0].title.set_text('$y$=-0.070 m')
ax[1].title.set_text('$y$=-0.035 m')
ax[2].title.set_text('$y$=-0.000 m')
ax[3].title.set_text('$y$= 0.035 m')
ax[4].title.set_text('$y$= 0.070 m')

# ax[0].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[1].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[2].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[3].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[4].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[0].set_ylabel('$z$ (m)', fontsize=10)

plt.savefig(figureOutputPath+'\Profiles_XYcompare.pdf', format='pdf', dpi=600)


"""
Turbulent kinetic energy comparison plot
"""


PIVscalar='tke'
cfdScalar1='tke2C'

#extract lines from mean PIV flow field data
PIV=meanPIV_treatedData['Scalar']['Case_8_'][PIVscalar]['avg']['db']
linePIV_neg07=piv.extractProfile(PIV, 'x', 'y', -70)
linePIV_neg035=piv.extractProfile(PIV, 'x', 'y', -35)
linePIV_0=piv.extractProfile(PIV, 'x', 'y', 0)
linePIV_pos035=piv.extractProfile(PIV, 'x', 'y', 35)
linePIV_pos07=piv.extractProfile(PIV, 'x', 'y', 70)

# pivLines=[linePIV_neg07,linePIV_neg035,linePIV_0,linePIV_pos035,linePIV_pos07]

linePIV_neg07.columns=['y','z',cfdScalar1]
linePIV_neg035.columns=['y','z',cfdScalar1]
linePIV_0.columns=['y','z',cfdScalar1]
linePIV_pos035.columns=['y','z',cfdScalar1]
linePIV_pos07.columns=['y','z',cfdScalar1]


#Create line plots for publication

#lateral components
fig, ax = plt.subplots(1,5,frameon=False,constrained_layout=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(3.0,forward=True)
ax[0].scatter(cfdLineData['line_neg07'][cfdScalar1],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='b')
ax[1].scatter(cfdLineData['line_neg035'][cfdScalar1],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='b')
ax[2].scatter(cfdLineData['line_0'][cfdScalar1],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='b')
ax[3].scatter(cfdLineData['line_pos035'][cfdScalar1],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='b')
lineCFD=ax[4].scatter(cfdLineData['line_pos07'][cfdScalar1],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='b',label='CFD')

ax[0].scatter(linePIV_neg07[cfdScalar1],linePIV_neg07['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[1].scatter(linePIV_neg035[cfdScalar1],linePIV_neg035['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[2].scatter(linePIV_0[cfdScalar1],linePIV_0['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[3].scatter(linePIV_pos035[cfdScalar1],linePIV_pos035['z']*0.001,s=20,facecolor='none',edgecolor='black')
linePIV=ax[4].scatter(linePIV_pos07[cfdScalar1],linePIV_pos07['z']*0.001,s=20,facecolor='none',edgecolor='black',label='PIV')

ax[4].legend(handles=[linePIV, lineCFD])

for i in range(len(ax)):
    #ax[i].set_xlim([-0.04,0.04])   
    ax[i].set_ylim([0,0.07])            
         

ax[0].title.set_text('$y$=-0.070 m')
ax[1].title.set_text('$y$=-0.035 m')
ax[2].title.set_text('$y$=-0.000 m')
ax[3].title.set_text('$y$= 0.035 m')
ax[4].title.set_text('$y$= 0.070 m')

# ax[0].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[1].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[2].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[3].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[4].set_xlabel('$w$ (m/s)', fontsize=10)
# ax[0].set_ylabel('$z$ (m)', fontsize=10)

plt.savefig(figureOutputPath+'\Profiles_TKEcompare.pdf', format='pdf', dpi=600)



"""
3D versus 2D Turbulent kinetic energy comparison plot
"""


PIVscalar='tke'
cfdScalar1='tke2C'
cfdScalar2='tke3C'

#extract lines from mean PIV flow field data
PIV=meanPIV_treatedData['Scalar']['Case_8_'][PIVscalar]['avg']['db']
linePIV_neg07=piv.extractProfile(PIV, 'x', 'y', -70)
linePIV_neg035=piv.extractProfile(PIV, 'x', 'y', -35)
linePIV_0=piv.extractProfile(PIV, 'x', 'y', 0)
linePIV_pos035=piv.extractProfile(PIV, 'x', 'y', 35)
linePIV_pos07=piv.extractProfile(PIV, 'x', 'y', 70)

# pivLines=[linePIV_neg07,linePIV_neg035,linePIV_0,linePIV_pos035,linePIV_pos07]

linePIV_neg07.columns=['y','z',cfdScalar1]
linePIV_neg035.columns=['y','z',cfdScalar1]
linePIV_0.columns=['y','z',cfdScalar1]
linePIV_pos035.columns=['y','z',cfdScalar1]
linePIV_pos07.columns=['y','z',cfdScalar1]


#Create line plots for publication

#lateral components
fig, ax = plt.subplots(1,5,frameon=False,constrained_layout=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(3.0,forward=True)
ax[0].scatter(cfdLineData['line_neg07'][cfdScalar1],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='b')
ax[1].scatter(cfdLineData['line_neg035'][cfdScalar1],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='b')
ax[2].scatter(cfdLineData['line_0'][cfdScalar1],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='b')
ax[3].scatter(cfdLineData['line_pos035'][cfdScalar1],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='b')
lineCFD1=ax[4].scatter(cfdLineData['line_pos07'][cfdScalar1],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='b',label='CFD_2Ctke')

ax[0].scatter(cfdLineData['line_neg07'][cfdScalar2],cfdLineData['line_neg07']['z'],s=20,facecolor='none',edgecolor='r')
ax[1].scatter(cfdLineData['line_neg035'][cfdScalar2],cfdLineData['line_neg035']['z'],s=20,facecolor='none',edgecolor='r')
ax[2].scatter(cfdLineData['line_0'][cfdScalar2],cfdLineData['line_0']['z'],s=20,facecolor='none',edgecolor='r')
ax[3].scatter(cfdLineData['line_pos035'][cfdScalar2],cfdLineData['line_pos035']['z'],s=20,facecolor='none',edgecolor='r')
lineCFD2=ax[4].scatter(cfdLineData['line_pos07'][cfdScalar2],cfdLineData['line_pos07']['z'],s=20,facecolor='none',edgecolor='r',label='CFD_3Ctke')


ax[0].scatter(linePIV_neg07[cfdScalar1],linePIV_neg07['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[1].scatter(linePIV_neg035[cfdScalar1],linePIV_neg035['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[2].scatter(linePIV_0[cfdScalar1],linePIV_0['z']*0.001,s=20,facecolor='none',edgecolor='black')
ax[3].scatter(linePIV_pos035[cfdScalar1],linePIV_pos035['z']*0.001,s=20,facecolor='none',edgecolor='black')
linePIV=ax[4].scatter(linePIV_pos07[cfdScalar1],linePIV_pos07['z']*0.001,s=20,facecolor='none',edgecolor='black',label='PIV_2Ctke')

ax[4].legend(handles=[linePIV, lineCFD1,lineCFD2])

for i in range(len(ax)):
    #ax[i].set_xlim([-0.04,0.04])   
    ax[i].set_ylim([0,0.07])            
         

ax[0].title.set_text('$y$=-0.070 m')
ax[1].title.set_text('$y$=-0.035 m')
ax[2].title.set_text('$y$=-0.000 m')
ax[3].title.set_text('$y$= 0.035 m')
ax[4].title.set_text('$y$= 0.070 m')

ax[0].set_xlabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)
ax[1].set_xlabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)
ax[2].set_xlabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)
ax[3].set_xlabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)
ax[4].set_xlabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)
ax[0].set_ylabel('$tke$ (m$^{2}$/s$^{2}$)', fontsize=10)

ax[0].set_ylabel('z (m)', fontsize=10)

plt.savefig(figureOutputPath+'\Profiles_3Dversus2DTKEcompare.pdf', format='pdf', dpi=600)
















"""
Testing zone - all this code might not run well ...
"""



"""
Spectral analysis of PIV data at 4D 
"""


#inst_PIV_path=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\Instantaneous\Case_8\Case_8_0.06s_inst_data"
from threading import Thread

#paths to lateral velocity time-series profiles
eq1_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_1_full_Hz_shot1_lat"
eq2_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_2_full_Hz_shot1_lat"
lo1_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_4_full_Hz_shot1_lat"
lo2_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_5_full_Hz_shot1_lat"
hi1_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_7_full_Hz_shot1_lat"
hi2_lat_path=r"C:\Users\Jason\Desktop\SpectralAnalysis\Lateral_timeSeries\Case_8_full_Hz_shot1_lat"

lateralPaths=[eq1_lat_path,eq2_lat_path,lo1_lat_path,lo2_lat_path,hi1_lat_path,hi2_lat_path]
names=['eq1_lat','eq2_lat','lo1_lat','lo2_lat','hi1_lat','hi2_lat']


#instLatData=[{} for x in lateralPaths]

instLatData={}
vArrays={}
    
def import_inst_lat_PIV_data(case,names,result,array,index):
    
    lateralDataFiles=piv.get_PIV_files(case)

    times=[]    
    for file in lateralDataFiles:
        times.append(os.path.splitext(os.path.basename(file))[0])

    instLatData[names[index]]=piv.load_instPIV_data(lateralDataFiles,times,'v',ref=True,shift=[0,25.46,0]) 
    
    array=[]
    for time in times:
        array.append(instLatData[names[index]][time]['v'])
    
    vArrays[names[index]]=array
    


threads = []

for i in range(len(lateralPaths)):
    
    process = Thread(target=import_inst_lat_PIV_data, args=[lateralPaths[i],names,instLatData,vArrays,i])
    process.start()
    threads.append(process)


instLatData['eq2_lat']['B00001']
vArrays['eq1_lat']





location=[0,35]


pivProbe0=piv.extractInstPointData(vArrays['eq1_lat'], location)
pivProbe1=piv.extractInstPointData(vArrays['eq2_lat'], location)
pivProbe2=piv.extractInstPointData(vArrays['lo1_lat'], location)
pivProbe3=piv.extractInstPointData(vArrays['lo2_lat'], location)
pivProbe4=piv.extractInstPointData(vArrays['hi1_lat'], location)
pivProbe5=piv.extractInstPointData(vArrays['hi2_lat'], location)

#np.savetxt('hi2_y0_z035_66.66Hz_w.csv', pivProbe2, fmt='%.7f', delimiter=',')



fig, ax = plt.subplots(6,1,frameon=False,constrained_layout=True,sharex=True,sharey=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(8,forward=True)

linePIV0=ax[0].plot(pivProbe0, c='blue',label='PIV',lw=1)
linePIV1=ax[1].plot(pivProbe1, c='blue',label='PIV',lw=1)
linePIV1=ax[2].plot(pivProbe2, c='blue',label='PIV',lw=1)
linePIV1=ax[3].plot(pivProbe3, c='blue',label='PIV',lw=1)
linePIV1=ax[4].plot(pivProbe4, c='blue',label='PIV',lw=1)
linePIV1=ax[5].plot(pivProbe5, c='blue',label='PIV',lw=1)


ax[5].xaxis.set_ticks([0,100,200,300,400,500,600,700])
ax[5].set_xticklabels(['0','6','12','18','24','30','36','42'])
ax[0].title.set_text('point: y=0.0, z=0.010 ')
ax[1].title.set_text('point: y=0.0, z=0.035 ')
ax[2].title.set_text('point: y=0.0, z=0.060 ')
ax[3].title.set_text('point: y=-0.06, z=0.010 ')
ax[4].title.set_text('point: y=-0.06, z=0.035 ')
ax[5].title.set_text('point: y=-0.06, z=0.060 ')



for axes in ax:
    axes.set_ylabel('$v$ (m/s)', fontsize=10)

plt.savefig(figureOutputPath+'\TimesSeries.pdf', format='pdf', dpi=600)



from scipy import signal


fig, ax = plt.subplots(6,1,frameon=False,constrained_layout=True,sharex=True,sharey=True)
fig.set_figwidth(8,forward=True)
fig.set_figheight(8,forward=True)



freq=66.6667

colors=['black','blue','darkgray']

nperseg=1024

f0, pxx_den0 = signal.welch(x=pivProbe0, fs=freq,nperseg=nperseg,detrend='constant')
f1, pxx_den1 = signal.welch(x=pivProbe1, fs=freq,nperseg=nperseg,detrend='constant')
f2, pxx_den2 = signal.welch(x=pivProbe2, fs=freq,nperseg=nperseg,detrend='constant')
f3, pxx_den3 = signal.welch(x=pivProbe3, fs=freq,nperseg=nperseg,detrend='constant')
f4, pxx_den4 = signal.welch(x=pivProbe4, fs=freq,nperseg=nperseg,detrend='constant')
f5, pxx_den5 = signal.welch(x=pivProbe5, fs=freq,nperseg=nperseg,detrend='constant')

ax[0].plot(f0,pxx_den0)
ax[1].plot(f1,pxx_den1)
ax[2].plot(f2,pxx_den2)
ax[3].plot(f3,pxx_den3)
ax[4].plot(f4,pxx_den4)
ax[5].plot(f5,pxx_den5)


ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_ylim([0.00000001,0.0004])
ax[0].set_xlim([0.0,5])


ax.set_xlim([0.1,10])  
ax[5].set_xlabel('$t$ (s)', fontsize=10)



