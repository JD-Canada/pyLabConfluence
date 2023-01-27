
"""
Code to answer Reviewer 2 of JFM article
- Performs circulation analysis
- LIF and circulation correlations
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sci
from matplotlib import rc

import pivToolbox.pivToolbox as piv
from scipy import ndimage
import re


rc('font',**{'family':'serif','serif':['Times New Roman'],'size':'7'})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amssymb}')



"""
Import vorticity and swirl strength data exported from DaVis
"""


vortPaths=r"Y:\LaboratoryConfluence\CirculationAnalysis\Hi2_vorticity"
swirlPaths=r"Y:\LaboratoryConfluence\CirculationAnalysis\Hi2_swirl"


#width and height of .dat file export (available in .dat file header)
I=164
J=104

vortFiles=piv.get_PIV_files(vortPaths)
swirlFiles=piv.get_PIV_files(swirlPaths)

# Needed to properly sequence images of files on my Ubuntu server

vortFiles.sort(key=lambda f: int(re.sub('\D', '', f)))
swirlFiles.sort(key=lambda f: int(re.sub('\D', '', f)))

times=[]
for file in vortFiles:
    times.append(os.path.splitext(os.path.basename(file))[0])


Data={}
Data['vort']=piv.load_instPIV_data(vortFiles,times,'vort')
Data['swirl']=piv.load_instPIV_data(swirlFiles,times,'swirl')



"""
Import binarized LIF data exported from DaVis
"""

lifPaths=r"Y:\LaboratoryConfluence\LIF of cases\Hi2_toExtractLIFtimeseries"

lifFiles=piv.get_PIV_files(lifPaths)

#This is needed to get properly sequenced images when the files are stored on my Ubuntu server
#If the same files are on my Windows computer, get_PIV_files returns a correctly organized list
lifFiles.sort(key=lambda f: int(re.sub('\D', '', f)))

timesLIF=[]
for file in lifFiles:
    timesLIF.append(os.path.splitext(os.path.basename(file))[0])

#Data={}
Data['LIF']=piv.load_instPIV_data(lifFiles,timesLIF,'LIF')





"""
Turn the Tecplot .dat imported data into a 2D array for processing with numpy

-Pixel lengths are available in header of .dat file
-long and short refer to the long and short sides of the image
-'I' and 'J' in the header correspond to the long and short side of the image respectively (must be specified above)

"""


for item in Data['vort']:
    Data['vort'][item]['array']=piv.convertDatTo2Dnumpy(Data['vort'][item]['vort'], I,J)
    
    #crop to avoid shear induced vorticity near bed
    Data['vort'][item]['array']=Data['vort'][item]['array'][0:160,45:97]


for item in Data['swirl']:
    Data['swirl'][item]['array']=piv.convertDatTo2Dnumpy(Data['swirl'][item]['swirl'], I,J)
    
    #crop to avoid shear induced vorticity near bed
    Data['swirl'][item]['array']=Data['swirl'][item]['array'][0:160,45:97]





"""
Extract time series from LIF data
"""

LIFsums=[]
for item in Data['LIF']:
    Data['LIF'][item]['array']=piv.convertDatTo2Dnumpy(Data['LIF'][item]['LIF'], 133,85)
    
    #crop to avoid shear induced vorticity near bed
    Data['LIF'][item]['array']=Data['LIF'][item]['array'][53:91,28:70] #from -5 cm to 1cm (3 cm on each side of -2 cm)

    LIFsums.append(Data['LIF'][item]['array'].sum())

LIFtimeSeries=np.asarray(LIFsums)
LIFtimeSeries=LIFtimeSeries/1596 #divide by number of pixels in sampling region







 
# fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
# fig.set_figwidth(6,forward=True)
# fig.set_figheight(3,forward=True)

# for item in Data['LIF']:
#     ax.imshow(np.rot90(np.fliplr(Data['LIF'][item]['array'])))
#     print(item)
#     plt.pause(0.1)



"""
Circulation analysis
"""

circData={}
circData['clock']={}
circData['clock']['circ']=[]
circData['clock']['circSum']=[]

circData['clock']['img']=[]
circData['clock']['area']=[]

circData['anti']={}
circData['anti']['circ']=[]
circData['anti']['circSum']=[]

circData['anti']['img']=[]
circData['anti']['area']=[]

circData['clock']['filt_c']=[]
circData['anti']['filt_antic']=[]

for item in Data['vort']:
    
    
    image=Data['vort'][item]['array']
    
    c, antic, img, anti_img, area, antiarea, circSum, anticircSum, filt_c, filt_antic = piv.getCirculation(Data['vort'][item]['array'], Data['swirl'][item]['array'], 1, 0.00127322**2,1e-4)
    
    
    circData['clock']['circ'].append(c)
    circData['anti']['circ'].append(antic)
    
    circData['clock']['circSum'].append(circSum)
    circData['anti']['circSum'].append(anticircSum)

    circData['clock']['img'].append(img)
    circData['anti']['img'].append(anti_img)

    circData['clock']['area'].append(area)
    circData['anti']['area'].append(antiarea) 

    circData['clock']['filt_c'].append(filt_c)
    circData['anti']['filt_antic'].append(filt_antic) 

timeSeriesClock_nonDim=np.asarray(circData['clock']['circSum'])*200.08 #200.02 = ((0.07/0.0714)/0.07**2) which non-dimensionalizes the m2/s units of circulation 
#timeSeriesFiltClock_nonDim=np.asarray(circData['clock']['filt_c'])*200.08 #200.02 = ((0.07/0.0714)/0.07**2) which non-dimensionalizes the m2/s units of circulation 

timeSeriesAnti_nonDim=np.asarray(circData['anti']['circSum'])*200.08








"""
Export circulation time series data for later use (avoid having to rerun this code multiple times)
"""


np.savetxt("Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_1_clock.txt", timeSeriesClock_nonDim, fmt='%f')
np.savetxt("Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_1_anti.txt", timeSeriesAnti_nonDim, fmt='%f')





"""
Import lateral velocity time series data for the Hi2 case at various points
"""

#brings in u and v Davis data (CAUTION: remember u is always the left-right component, and v is always the bottom-up component in DaVis exports)
#therefore on the 4D PIV plane in the confluence u in DaVis is actually the lateral (v) component, and v in DaVis is therefore the vertical (w) component
vPaths=r"Y:\LaboratoryConfluence\CirculationAnalysis\Hi2_lateral"
vFiles=piv.get_PIV_files(vPaths)

vFiles.sort(key=lambda f: int(re.sub('\D', '', f)))


label=['Case_8_']
times=[]

for file in vFiles:
    times.append(os.path.splitext(os.path.basename(file))[0])

Data['v']=piv.load_instPIV_data(vFiles,times,'v',ref=True,shift=[0,25.46,0]) 


"""
Average lateral velocity over the same sampling volume as the LIF
"""


vMean=[]
for item in Data['v']:
    Data['v'][item]['array']=piv.convertDatTo2Dnumpy(Data['v'][item]['v'][:,:-1], 164,104) #[:,:-1] is to remove an extra column that appears for some reason ...
    
    #crop to avoid shear induced vorticity near bed
    Data['v'][item]['top']=Data['v'][item]['array'][50:100,41:69] #from -5 cm to 1cm (3 cm on each side of -2 cm)
    Data['v'][item]['bottom']=Data['v'][item]['array'][50:100,69:97]
   
vMeanTop=[]
vMeanBottom=[]

for item in Data['v']:
   
    vMeanTop.append(Data['v'][item]['top'].mean())
    vMeanBottom.append(Data['v'][item]['bottom'].mean())
    
    
    
latTop=np.asarray(vMeanTop)
latBottom=np.asarray(vMeanBottom)
#latTimeSeries=latTimeSeries/1596 #divide by number of pixels in sampling region

latTop = (latTop - latTop.mean())/latTop.std()
latBottom = (latBottom- latBottom.mean())/latBottom.std()


plt.imshow(Data['v']['B00001']['array'])

plt.imshow(Data['v']['B00001']['bottom'])

print('hi')

vArrays=[]

for time in times:
    vArrays.append(Data['v'][time]['v'])
        

point_top=[-0,50]   
point_bottom=[-0,20]   

top=piv.extractInstPointData(vArrays,point_top)
bottom=piv.extractInstPointData(vArrays,point_bottom)


topNor=(top-top.mean())/top.std()
bottomNor=(bottom-bottom.mean())/bottom.std()

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

topNor_smoothed=moving_average(topNor,5)
bottomNor_smoothed=moving_average(bottomNor,5)



from scipy import stats


"""
Pearson correlations between the various time series
"""

#bottom v versus LIF
stats.pearsonr(bottomNor_smoothed,LIFtimeSeries[:-5] ) #had to remove a few elements

#top v versus LIF
stats.pearsonr(latTimeSeries,LIFtimeSeries[:-1] ) #had to remove a few elements

#spatially averaged lateral versus clockwise circulation
stats.pearsonr(latTimeSeries,timeSeriesClock_nonDim ) #had to remove a few elements

stats.pearsonr(latTimeSeries,latTimeSeries) #had to remove a few elements

"""
Time series plots for Hi2 case
"""
fig, ax = plt.subplots(3,1,frameon=False,constrained_layout=True)
fig.set_figwidth(5,forward=True)
fig.set_figheight(3,forward=True)

ax[0].plot(LIFtimeSeries,color='black',lw=0.5, label='$sum(I)$')

ax[0].vlines(x=[145,342,785,1103,1312,1592,1819,2041], ymin=0, ymax=30, colors='gray',ls='--',lw=0.5)

ax[0].set_ylim([0,1])
ax[0].set_xlim([0,2720])
ax[0].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[0].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[0].set_xlabel('$\\tilde{t}$')
ax[0].set_ylabel('sum(I)')


# ax[1].plot(latTop,color='red',lw=0.5, label='$sum(I)$')
ax[1].plot(latBottom,color='blue',lw=0.5, label='$sum(I)$')

ax[1].vlines(x=[145,342,785,1103,1312,1592,1819,2041], ymin=-4, ymax=4, colors='gray',ls='--',lw=0.5)

ax[1].set_ylim([-3,2])
ax[1].set_xlim([0,2720])
ax[1].yaxis.set_ticks([-4,-2,0,2,4])
ax[1].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[1].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[1].set_xlabel('$\\tilde{t}$')
ax[1].set_ylabel('$\\tilde{v} / \\sigma_{\\tilde{v}}$')



ax[2].plot(-timeSeriesClock_nonDim,color='blue',lw=0.5,label='$|sum(\\tilde{\\Gamma}_{\curvearrowright})|$')
ax[2].plot(timeSeriesAnti_nonDim,color='red',lw=0.5,label='$sum(\\tilde{\\Gamma}_{\curvearrowleft})$')
ax[2].vlines(x=[145,342,785,1103,1312,1592,1819,2041], ymin=0, ymax=3, colors='gray',ls='--',lw=0.5)

ax[2].set_ylim([0,3.2])
ax[2].set_xlim([0,2720])
ax[2].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[2].yaxis.set_ticks([0,1,2,3])
ax[2].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[2].set_xlabel('$\\tilde{t}$')
ax[2].set_ylabel('$|sum(\\tilde{\\Gamma})|$')



#ax[1].hlines(y=[2.8], xmin=145, xmax=175, colors='gray',ls='-',lw=0.5)
#[1].text(240,2.8,'lag', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')



#ax[2].set_ylabel('$|sum(\\tilde{\\Gamma})|$')


ax[0].text(50,0.85,'a)', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')
ax[1].text(50,2.7,'b)', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')
ax[2].text(50,3.2,'c)', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')


# ax[1].annotate("lag", xy=(145, 2.7), xytext=(175, 2.7 ),arrowprops=dict(arrowstyle="<-"))

#ax[1].text(500,2.4,'$|sum(\\tilde{\\Gamma}_{\circlearrowright})|$', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')
#ax[1].text(500,0.4,'$sum(\\tilde{\\Gamma}_{\circlearrowleft})$', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')

#ax[0].legend(frameon=False,bbox_to_anchor = [0.97, 0.2])
ax[1].legend(frameon=False)




plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_Hi1LIFGamma.pdf', format='pdf', dpi=600)








"""
Lines plots of total circulation
"""

fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)

ax.plot(-timeSeriesClock)
ax.plot(timeSeriesAnti)



ax.set_ylim([-15000,15000])
ax.set_xlim([0,60])



"""
Animation of labelled blobs
"""

fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)

for i in range(2767):
    ax.cla()
    plt.imshow(circData['anti']['img'][i]) 

    # ax.plot(circData['clock']['circ'][i])
    
    # ax.plot(circData['anti']['circ'][i][(circData['anti']['area'][i]<200)])
    
    #n, bins, patches = ax.hist(x=circData['clock']['area'][i], bins='auto', color='#0504aa',alpha=0.7, rwidth=0.85)
    
    print(i)
    plt.pause(0.1)
   

"""
Focus on total clockwise circulation, but only in the SOV zone
"""

sovData={}
sovData['clock']={}
sovData['clock']['circ']=[]
sovData['clock']['circSum']=[]

sovData['clock']['img']=[]
sovData['clock']['area']=[]

sovData['anti']={}
sovData['anti']['circ']=[]
sovData['anti']['circSum']=[]

sovData['anti']['img']=[]
sovData['anti']['area']=[]


for item in Data['vort']:
    
   
    
    c, antic, l, antil, area, antiarea, circSum, anticircSum = getCirculation(Data['vort'][item]['array'][0:50,:], 2, 1.27322,50)
    
    
    sovData['clock']['circ'].append(c)
    sovData['anti']['circ'].append(antic)
    
    sovData['clock']['circSum'].append(circSum)
    sovData['anti']['circSum'].append(anticircSum)

    sovData['clock']['img'].append(l)
    sovData['anti']['img'].append(antil)

    sovData['clock']['area'].append(area)
    sovData['anti']['area'].append(antiarea) 


fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)

for i in range(2767):
    ax.cla()
    plt.imshow(sovData['clock']['img'][i]) 

    # ax.plot(circData['clock']['circ'][i])
    
    # ax.plot(circData['anti']['circ'][i][(circData['anti']['area'][i]<200)])
    
    #n, bins, patches = ax.hist(x=circData['clock']['area'][i], bins='auto', color='#0504aa',alpha=0.7, rwidth=0.85)
    
    
    plt.pause(0.1)



#plt.imshow(Data['vort']['B00001']['array'])#[0:60,:])#


#plt.imshow(Data['vort'][item]['array'][0:60,:])#


sovCirculation=np.asarray(sovData['clock']['circSum'])

fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)

ax.plot((sovCirculation-sovCirculation.mean())/sovCirculation.mean())
ax.plot((LIFtimeSeries-LIFtimeSeries.mean())/LIFtimeSeries.mean())




"""
Cross-correlation analysis of total circulation and LIF time series

This provides a good Python explanation for autocorrelation and other methods:
    https://scicoding.com/cross-correlation-in-python-3-popular-packages/
    https://scicoding.com/4-ways-of-calculating-autocorrelation-in-python/
"""

import numpy


mean_LIF = np.mean(LIFtimeSeries)
var_LIF = np.var(LIFtimeSeries)
nor_LIFts = (LIFtimeSeries - mean_LIF)

mean_anti = np.mean(timeSeriesAnti_nonDim)
var_anti = np.var(timeSeriesAnti_nonDim)
nor_anti = timeSeriesAnti_nonDim - mean_anti

mean_clock = np.mean(timeSeriesClock_nonDim)
var_clock = np.var(timeSeriesClock_nonDim)
nor_clock = timeSeriesClock_nonDim - mean_clock


corr = np.correlate(nor_LIFts, nor_LIFts, 'full')#[len(timeSeriesClock)-1:] 
corr = corr/ (len(nor_LIFts)*np.std(nor_LIFts) * np.std(nor_LIFts)) # Normalization


#LIF and clockwise
#This indicates a lag of about 41 timesteps (0.205 s) between peak LIF and peak clockwise circulation
corr1 = np.correlate(nor_LIFts, nor_clock, 'full')#[len(timeSeriesClock)-1:] 
corr1 = corr1/ (len(nor_LIFts)*np.std(nor_clock) * np.std(nor_LIFts)) # Normalization


corr2 = np.correlate(nor_clock, nor_clock, 'full')#[len(timeSeriesClock)-1:] 
corr2 = corr2/ (len(nor_clock)*np.std(nor_clock) * np.std(nor_clock)) # Normalization


corr3 = np.correlate(nor_anti, nor_anti, 'full')#[len(timeSeriesClock)-1:] 
corr3 = corr3/ (len(nor_anti)*np.std(nor_anti) * np.std(nor_anti)) # Normalization

len(nor_LIFts)

np.std(nor_clock)
np.std(nor_LIFts)

#Anticlockwise circulation is also lagged by 41 timesteps
corr2 = np.correlate(nor_LIFts, nor_anti, 'same')#[len(timeSeriesClock)-1:] 
corr2 = corr2/ (len(nor_LIFts) * np.std(nor_anti) * np.std(nor_LIFts)) # Normalization

#Shows nor_clock and nor_anti are most correlated at 0 lag
corr3 = np.correlate(nor_clock, nor_anti, 'full')#[len(timeSeriesClock)-1:] 
corr3 = corr3/ (len(nor_clock) * np.std(nor_clock) * np.std(nor_anti)) # Normalization

fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(5,forward=True)
fig.set_figheight(3,forward=True)

ax.plot(corr,lw=0.5,color='black',label='$R_{II}$')
ax.plot(corr2,lw=0.5,color='darkgray',label='$R_{\\tilde{\\Gamma}_{\\circlearrowright}\\tilde{\\Gamma}_{\\circlearrowright}}$',ls='dashdot')
ax.plot(corr1,lw=0.5,color='blue',label='$R_{I\\tilde{\\Gamma}_{\\circlearrowright}}$')


ax.set_ylim([-0.412,1])
ax.set_xlim([2333,3201])

ax.vlines(x=[2767], ymin=-0.412, ymax=1, colors='gray',ls='--',lw=0.5)
ax.hlines(y=[0], xmin=2087, xmax=3447, colors='gray',ls='--',lw=0.5)



ax.xaxis.set_ticks([2333,2420,2507,2593,2680,2767,2854,2941,3027,3114,3201])
ax.set_xticklabels(['-6.5','-5.2','-3.9','-2.6','-1.3','0.0','1.3','2.6','3.9','5.2','6.5'])
ax.set_xlabel('$\\tilde{t}$')
ax.set_ylabel('$R_{ii}$')


ax.legend(frameon=False)

plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_Correlations.pdf', format='pdf', dpi=600)


# acorr = np.correlate(nor_LIFts, nor_clock, 'full')

"""
Power spectral density analysis:
    
    Good link describing PSD in Python:
        https://scicoding.com/calculating-power-spectral-density-in-python/
"""



import scipy.signal


#(f, S) = scipy.signal.periodogram(nor_LIFts, 66.666, scaling='density')

(f_lif, S_lif)= scipy.signal.welch(nor_LIFts, 66.6666, nperseg=2*1024)
(f_clock, S_clock)= scipy.signal.welch(nor_clock, 66.6666, nperseg=2*1024)
(f_anti, S_anti)= scipy.signal.welch(nor_anti, 66.6666, nperseg=2*1024)



fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(6,forward=True)
fig.set_figheight(3,forward=True)

ax.plot(f_lif,S_lif)
ax.plot(f_clock,S_clock)
ax.plot(f_anti,S_anti)


ax.set_ylim([10e-2,10e8])
ax.set_xlim([0,1])
ax.set_yscale('log')
ax.set_xscale('log')



"""
Pearson correlation analysis
"""

    
pear1=scipy.stats.pearsonr(-nor_clock/mean_clock, nor_anti/mean_anti)#0.45 correlation, highly significant

    





"""
Circulation timeseries figures
"""
Eq1_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_1_clock.txt")
Eq1_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_1_anti.txt")

Eq1_clock.mean()
StreaEq1_anti.mean()

Lo1_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Lo_1_clock.txt")
Lo1_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Lo_1_anti.txt")

Lo1_clock.mean()
Lo1_anti.mean()

Hi1_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Hi_1_clock.txt")
Hi1_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Hi_1_anti.txt")


Eq2_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_2_clock.txt")
Eq2_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Eq_2_anti.txt")

Eq2_clock.mean()
Eq2_anti.mean()

Lo2_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Lo_2_clock.txt")
Lo2_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Lo_2_anti.txt")

Lo2_clock.mean()
Lo2_anti.mean()

Hi2_clock = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Hi_2_clock.txt")
Hi2_anti = np.loadtxt(r"Y:\LaboratoryConfluence\CirculationAnalysis\CirculationData\Hi_2_anti.txt")

Hi2_clock.mean()

Hi2_anti.mean()

"""
Six panel circulation time series figure
"""


fig, ax = plt.subplots(3,2,frameon=False,constrained_layout=True)
fig.set_figwidth(5,forward=True)
fig.set_figheight(4,forward=True)


lw = 0.25
i=0
j=0

ax[i][j].plot(-Eq1_clock,color='blue',lw=lw)
ax[i][j].plot(Eq1_anti,color='red',lw=lw)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
#ax[i][j].set_xlabel('$\\tilde{t}$')


i=0
j=1

ax[i][j].plot(-Eq2_clock,color='blue',lw=lw,label='clockwise')
ax[i][j].plot(Eq2_anti,color='red',lw=lw, label='anticlockwise')

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
#ax[i][j].set_xlabel('$\\tilde{t}$')



i=1
j=0

ax[i][j].plot(-Lo1_clock,color='blue',lw=lw)
ax[i][j].plot(Lo1_anti,color='red',lw=lw)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
#ax[i][j].set_xlabel('$\\tilde{t}$')


i=1
j=1

ax[i][j].plot(-Lo2_clock,color='blue',lw=lw)
ax[i][j].plot(Lo2_anti,color='red',lw=lw)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
#ax[i][j].set_xlabel('$\\tilde{t}$')


i=2
j=0

ax[i][j].plot(-Hi1_clock,color='blue',lw=lw)
ax[i][j].plot(Hi1_anti,color='red',lw=lw)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[i][j].set_xlabel('$\\tilde{t}$')


i=2
j=1

ax[i][j].plot(-Hi2_clock,color='blue',lw=lw)
ax[i][j].plot(Hi2_anti,color='red',lw=lw)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[i][j].set_xlabel('$\\tilde{t}$')

coords=[[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]]


xRemoveTicks=[[0,0],[0,1],[1,0],[1,1]]
yRemoveTicks=[[0,1],[1,1],[2,1]]


for i in xRemoveTicks:
    ax[i[0]][i[1]].xaxis.set_major_locator(plt.NullLocator())

for i in yRemoveTicks:
    ax[i[0]][i[1]].yaxis.set_major_locator(plt.NullLocator())

#Fd=['n.d.','n.d.','n.d.','4.72','2.62','1.31','3.34','1.85','0.93','2.72','1.51','0.73']
#cases=['1','2','4','5','7','8']

letters=['a) $eq_1$','b) $eq_2$','c) $lo_1$','d) $lo_2$','e) $hi_1$','f) $hi_2$']

for index, i in enumerate(coords):

    #ax[i[0]][i[1]].add_patch(Rectangle((75,0), 25, 12, facecolor="white",zorder=2))
    ax[i[0]][i[1]].text(220,2.75,"%s" % letters[index], fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')

plt.setp(ax[-1, :], xlabel='$\\tilde{t}$')
plt.setp(ax[:, 0], ylabel='$|sum(\\tilde{\\Gamma})|$')

ax[0][1].legend(frameon=False)


plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_CircTimeSeries.pdf', format='pdf', dpi=600)


#ax[0].set_ylabel('sum(I)')


ax[0].vlines(x=[145,1103,1312,1592,1819,2041], ymin=0, ymax=30, colors='gray',ls='--',lw=0.5)

ax[1].plot(-timeSeriesClock_nonDim,color='blue',lw=0.5,label='$|sum(\\tilde{\\Gamma}_{\circlearrowright})|$')
ax[1].plot(timeSeriesAnti_nonDim,color='red',lw=0.5,label='$sum(\\tilde{\\Gamma}_{\circlearrowleft})$')

ax[1].set_ylim([0,3.2])
ax[1].set_xlim([0,2720])
ax[1].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[1].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
ax[1].set_xlabel('$\\tilde{t}$')
ax[1].set_ylabel('$sum(\\tilde{\\Gamma})$')

ax[1].vlines(x=[145,1103,1312,1592,1819,2041], ymin=0, ymax=30, colors='gray',ls='--',lw=0.5)


ax[0].text(50,0.9,'a)', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')
ax[1].text(50,2.7,'b)', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')

#ax[1].text(500,2.4,'$|sum(\\tilde{\\Gamma}_{\circlearrowright})|$', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')
#ax[1].text(500,0.4,'$sum(\\tilde{\\Gamma}_{\circlearrowleft})$', fontsize=6, color='black', weight='bold', zorder=3 , ha='center', va='center')

#ax[0].legend(frameon=False,bbox_to_anchor = [0.97, 0.2])
ax[1].legend(frameon=False)




plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_Hi1LIFGamma.pdf', format='pdf', dpi=600)






Eq1_clock_mean = Eq1_clock.mean()
Eq1_anti_mean = Eq1_anti.mean()

Lo1_clock_mean = Lo1_clock.mean()
Lo1_anti_mean = Lo1_anti.mean()

Hi1_clock_mean = Hi1_clock.mean()
Hi1_anti_mean = Hi1_anti.mean()


Eq2_clock_mean = Eq2_clock.mean()
Eq2_anti_mean = Eq2_anti.mean()

Lo2_clock_mean = Lo2_clock.mean()
Lo2_anti_mean = Lo2_anti.mean()

Hi2_clock_mean = Hi2_clock.mean()
Hi2_anti_mean = Hi2_anti.mean()



Lo1_clock_mean_nor = Lo1_clock.mean()/Eq1_clock_mean
Lo1_anti_mean_nor = Lo1_anti.mean()/Eq1_clock_mean

Hi1_clock_mean_nor = Hi1_clock.mean()/Eq1_clock_mean
Hi1_anti_mean_nor = Hi1_anti.mean()/Eq1_clock_mean


Eq2_clock_mean = Eq2_clock.mean()
Eq2_anti_mean = Eq2_anti.mean()

Lo2_clock_mean_nor = Lo2_clock.mean()/Eq2_clock_mean
Lo2_anti_mean_nor = Lo2_anti.mean()/Eq2_clock_mean

Hi2_clock_mean_nor = Hi2_clock.mean()/Eq2_clock_mean
Hi2_anti_mean_nor = Hi2_anti.mean()/Eq2_clock_mean



equal=[1,Lo1_clock_mean_nor,Hi1_clock_mean_nor]
unequal=[1,Lo2_clock_mean_nor,Hi2_clock_mean_nor]


fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
fig.set_figwidth(5,forward=True)
fig.set_figheight(4,forward=True)


density=[0.00,0.33,0.66]

ax.scatter(density,[Eq1_clock_mean, Lo1_clock_mean, Hi1_clock_mean],color='gray')
ax.scatter(density,[Eq2_clock_mean, Lo2_clock_mean, Hi2_clock_mean],color='black')

ax.scatter(density,unequal,color='red')

ax.errorbar(density,equal, yerr=)

ax[i][j].set_ylim([0,3])
ax[i][j].set_xlim([0,2720])
ax[i][j].xaxis.set_ticks([0,340,680,1020,1360,1700,2040,2380,2720])
ax[i][j].set_xticklabels(['0','5','10','15','20','25','30','35','40'])
#ax[i][j].set_xlabel('$\\tilde{t}$')


























# """
# Detect vortices using OpenCV

# Ideally, you save the data from each timestep into a dataframe or something in this processing step
# """

# for item in Data['sw']:
    
#     image=Data['sw'][item]['img']
#     color = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
#     (T, threshInv) = cv2.threshold(image, 20, 255, cv2.THRESH_BINARY)
    
    
#     contours, _ = cv2.findContours(threshInv, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#     print(len(contours), "objects were found in this image.")
    

#     cntarea=[]
#     xcntcoord=[]
#     ycntcoord=[]
#     filteredContours=[]

#     #loop through contours
#     for c in contours:
        
        
#         if cv2.contourArea(c) < 10:          
#             continue                          
        
#         filteredContours.append(c)
#         cntarea.append(cv2.contourArea(c)) 
#         M = cv2.moments(c)                 
#         cx = int(M['m10']/M['m00'])        
#         cy = int(M['m01']/M['m00'])        
#         xcntcoord.append(cx)
#         ycntcoord.append(cy)


#     #remove small contours
#     # contourList=list(contours)
#     # def by_size(contourList, size):
#     #     return [contour for contour in contourList if len(contour) > size]
#     # filteredContours=by_size(contourList,3)
    
#     cv2.drawContours(color, filteredContours, -1, (0, 255, 0), 1)
    
#     cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
#     cv2.resizeWindow("Resized_Window", 1000, 700)
    
#     color=cv2.flip(color,0)
#     color=cv2.rotate(color,cv2.ROTATE_90_CLOCKWISE)
    
#     cv2.imshow("Resized_Window", color)    
    
    
    
    
#     key = cv2.waitKey(2) & 0xFF
#     if key == ord("q"):
#         break
#     print('shown')
# #test=convertDatToOpenCVImage(Data['sw']['B00001']['sw'], 164,104).astype(np.uint8)

# len(contourList[10])






















