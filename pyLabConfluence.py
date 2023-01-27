"""
Code to treat Figures 1 to 4 in the JFM article (time averaged data)
and Figure 6 (instantaneous vectors over LIF background)

Code was tested and fully functional on January 18th, 2023
"""

import pivToolbox.pivToolbox as piv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle
import os


rc('font',**{'family':'serif','serif':['Times New Roman'],'size':'7'})
rc('text', usetex=True)



""" Step 1 - Bring in mean vector field data


    **Important - only export scalar data from DaVis (i.e. lateral and vertical velocity components in separate .dat files**

"""

#modify as needed
path=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\PlanarPIV_at_4D"


files=piv.get_PIV_files(path)

"""Provide data labels and load data to pandas dataframes"""
fieldKeys=['B00001','B00002','B00003','B00004','B00005','B00006','sw','vor'] #u,w,TKE, Rvw,sw,vor
fieldNames=['u','v','stdv','stdw','tke','rss','sw','vor']
shot_labels=['Shot_1','Shot_2']

#These are the names of the individual folders in path, you might only need one
plane_labels=['Case_1_',
              'Case_2_',
              'Case_4_',
              'Case_5_',
              'Case_7_',
              'Case_8_']



#main dictionary holding untreated data (i.e. output from DaVis)
data=piv.load_PIV_dataframes(files,plane_labels,shot_labels,fieldNames,fieldKeys)



"""
Interpolate vector data to a regular grid 
"""
treatedData={}
treatedData['Vector']=piv.treatMeanVecData(data,[0,25.46,0],[-110,110,0,80],2.5,avg=True,norm=False)


"""
Interpolate scalar data to a regular grid 
"""
scalarNames=['stdv','stdw','tke','rss','sw','vor']
treatedData['Scalar']=piv.treatScalarData(data,scalarNames,[0,25.46,0],[-110,110,0,80],2.5,avg=True)


"""
Define some parameters needed to non-dimensionalize data
"""

velocityScale=0.071 #bulk velocity of cold channel
lengthScale=0.07    #depth in the mixing interface







"""
Step 2 - Make mean flow field figures
"""


def commonSubPlotElements():
    
    """
    Function to place lettering and common elements on the subplots of the various figures
    """
    
    plt.setp(ax, xticks=[-98, -49,0,49,98], xticklabels=['-1.40', '-0.70', '0.00','0.70','1.40'],yticks=[0, 35, 70],yticklabels=['0', '0.5', '1'])
    
    for index, i in enumerate(coords):
        ax[i[0]][i[1]].set_xlim([-98,98])            
        ax[i[0]][i[1]].set_ylim([0,70])    
    
    xRemoveTicks=[[0,0],[0,1],[1,0],[1,1]]
    yRemoveTicks=[[0,1],[1,1],[2,1]]
    
    
    for i in xRemoveTicks:
        ax[i[0]][i[1]].xaxis.set_major_locator(plt.NullLocator())
    
    for i in yRemoveTicks:
        ax[i[0]][i[1]].yaxis.set_major_locator(plt.NullLocator())
    
    
    letters=['a) $eq_1$','b) $eq_2$','c) $lo_1$','d) $lo_2$','e) $hi_1$','f) $hi_2$']
    
    for index, i in enumerate(coords):
    
        ax[i[0]][i[1]].add_patch(Rectangle((65,0), 33, 12, facecolor="white",zorder=2))
        ax[i[0]][i[1]].text(82,5,"%s" % letters[index], fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')
    
    plt.setp(ax[-1, :], xlabel='$\\tilde{y}$')
    plt.setp(ax[:, 0], ylabel='$\\tilde{z}$')





"""
Vector magnitude and streamline figure (Fig. 2 in article)
"""

coords=[[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]]
cases=['Case_1_','Case_2_','Case_4_','Case_5_','Case_7_','Case_8_']

fig, ax = plt.subplots(3,2,frameon=False,constrained_layout=True)

fig.set_figwidth(5.0,forward=True)
fig.set_figheight(4,forward=True)



for index, case in enumerate(cases):
    
    
    if case in treatedData['Vector']:
        print(case)

        u=treatedData['Vector'][case]['avg']['u_grid']/velocityScale
        v=treatedData['Vector'][case]['avg']['v_grid']/velocityScale     
        color=np.sqrt(np.multiply(u,u)+np.multiply(v,v))

        X=treatedData['Vector'][case]['avg']['X']
        Y=treatedData['Vector'][case]['avg']['Y']  
        
        import matplotlib.colors
        norm = matplotlib.colors.Normalize(vmin=0,vmax=0.65,clip=False)
        
        levels = np.linspace(0.0, 0.65, 8)
        CS=ax[coords[index][0]][coords[index][1]].contourf(X,Y, color,cmap='bone_r',alpha=0.6, norm=norm,extend='neither',levels=levels) 
        
        ax[coords[index][0]][coords[index][1]].streamplot(X,Y,u,v, density=1.2,color='black',linewidth=0.5,arrowsize=0.3)
        ax[coords[index][0]][coords[index][1]].set_xlim([-100,100])            
        ax[coords[index][0]][coords[index][1]].set_ylim([0,70])    
        
    else:
        pass

plt.colorbar(CS, ax=ax, shrink=0.6, format='%.2f', location='bottom', pad=0.005, label='$\\tilde{U}_{vw}$')

commonSubPlotElements()

plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_streamlines.pdf', format='pdf', dpi=600)




"""
Swirl figure (Fig. 3 of article)
"""

coords=[[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]]
cases=['Case_1_','Case_2_','Case_4_','Case_5_','Case_7_','Case_8_']

fig, ax = plt.subplots(3,2,frameon=False,constrained_layout=True)

fig.set_figwidth(5.0,forward=True)
fig.set_figheight(4,forward=True)


for index, case in enumerate(cases):
    
    if case in treatedData['Scalar']:
        print(case)
        
        swirl=treatedData['Scalar'][case]['sw']['avg']['grid']*((lengthScale/velocityScale)**2)
        
        u=treatedData['Vector'][case]['avg']['u_grid']/velocityScale
        v=treatedData['Vector'][case]['avg']['v_grid']/velocityScale     
        

        X=treatedData['Vector'][case]['avg']['X']
        Y=treatedData['Vector'][case]['avg']['Y']  

        X=treatedData['Scalar'][case]['sw']['avg']['X']
        Y=treatedData['Scalar'][case]['sw']['avg']['Y']
        
        CS=ax[coords[index][0]][coords[index][1]].contourf(X,Y,swirl,cmap='bone_r')
        ax[coords[index][0]][coords[index][1]].streamplot(X,Y,u,v, density=1,color='black',linewidth=0.5,arrowsize=0.3)

        
    else:
        pass

plt.colorbar(CS, ax=ax,shrink=0.6,format='%.e',location='bottom',pad=0.005, label='$\overline{\\tilde{\lambda}}$')

commonSubPlotElements()

plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_swirl.pdf', format='pdf', dpi=600)



"""
TKE figure (Fig. 4 of article) - ATTENTION - This only provides the basics. The rest was done in Inkscape.
"""

coords=[[0,0],[1,0],[2,0],[3,0],[4,0],[5,0]]
cases=['Case_1_','Case_2_','Case_4_','Case_5_','Case_7_','Case_8_']

fig, ax = plt.subplots(6,3,frameon=False,constrained_layout=True)

fig.set_figwidth(5.3,forward=True)
fig.set_figheight(4,forward=True)


for index, case in enumerate(cases):
    
    if case in treatedData['Scalar']:
        print(case)
        
        stdv=0.5*(treatedData['Scalar'][case]['stdv']['avg']['grid']**2)/(velocityScale**2)
        stdw=0.5*(treatedData['Scalar'][case]['stdw']['avg']['grid']**2)/(velocityScale**2)
        TKE=stdv+stdw

        X=treatedData['Scalar'][case]['stdv']['avg']['X']
        Y=treatedData['Scalar'][case]['stdv']['avg']['Y']
        
        levels = np.linspace(0.0, 0.06, 8)
       
        CS=ax[coords[index][0]][0].contourf(X,Y,stdv,cmap='bone_r',extend='neither',levels=levels)
        CS=ax[coords[index][0]][1].contourf(X,Y,stdw,cmap='bone_r',extend='neither',levels=levels)
        CS=ax[coords[index][0]][2].contourf(X,Y,TKE,cmap='bone_r',extend='neither',levels=levels)

    else:
        pass

   
plt.setp(ax, 
         xticks=[-98, -49,0,49,98], 
         xticklabels=['-1.40', '-0.70', '0.00','0.70','1.40'],
         yticks=[0, 35, 70],
         yticklabels=['0', '0.5', '1'])

for index, i in enumerate(coords):
    ax[i[0]][0].set_xlim([-98,98])   
    ax[i[0]][1].set_xlim([-98,98]) 
    ax[i[0]][2].set_xlim([-98,98]) 
    
    ax[i[0]][0].set_ylim([0,70]) 
    ax[i[0]][1].set_ylim([0,70])
    ax[i[0]][2].set_ylim([0,70])
          

xRemoveTicks=[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2],[3,0],[3,1],[3,2],[4,0],[4,1],[4,2],[5,0],[5,1],[5,2]]
yRemoveTicks=[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2],[3,0],[3,1],[3,2],[4,0],[4,1],[4,2],[5,0],[5,1],[5,2]]


for i in xRemoveTicks:
    ax[i[0]][i[1]].xaxis.set_major_locator(plt.NullLocator())

for i in yRemoveTicks:
    ax[i[0]][i[1]].yaxis.set_major_locator(plt.NullLocator())


plt.subplots_adjust(wspace=0.00,hspace=0.000)

plt.colorbar(CS, ax=ax,shrink=0.6,format='%.e',location='bottom',pad=0.005, label='$\\tilde{k}_{v}$, $\\tilde{k}_{w}$, $\\tilde{k}_{vw}$ [-]')


plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_stdv.pdf', format='pdf', dpi=600)  








"""
Instantaneous vector magnitude and LIF figure
"""


def commonSubPlotElementsInstVec():
    
    
    xRemoveTicks=[[0,0],[0,1],[1,0],[1,1]]
    yRemoveTicks=[[0,1],[1,1],[2,1]]
    
    topLabels=[[0,0],[0,1]]
    topText=['$V_r = 1.00$','$V_r = 0.55$']
    for index, i in enumerate(topLabels):
        ax[i[0]][i[1]].xaxis.set_label_position("top")
        ax[i[0]][i[1]].set_xlabel(topText[index])
    
    sideLabels=[[0,1],[1,1],[2,1]]
    sideText=['$\\tilde{\\Delta\\rho}_{0.00}$','$\\tilde{\\Delta\\rho}_{0.33}$','$\\tilde{\\Delta\\rho}_{0.66}$']
    for index, i in enumerate(sideLabels):
        
        ax[i[0]][i[1]].yaxis.set_label_position("right")
        ax[i[0]][i[1]].set_ylabel(sideText[index])
    
    for i in xRemoveTicks:
        ax[i[0]][i[1]].xaxis.set_major_locator(plt.NullLocator())
    
    for i in yRemoveTicks:
        ax[i[0]][i[1]].yaxis.set_major_locator(plt.NullLocator())
    
    
    letters=['(a)','(b)','(c)','(d)','(e)','(f)']
    
    for index, i in enumerate(coords):
    
        ax[i[0]][i[1]].add_patch(Rectangle((75,0), 25, 15, facecolor="white",zorder=2))
        ax[i[0]][i[1]].text(87.5,7,"%s" % letters[index], fontsize=8, color='black', weight='bold', zorder=3 , ha='center', va='center')
    
    plt.setp(ax[-1, :], xlabel='$\\tilde{y}$')
    plt.setp(ax[:, 0], ylabel='$\\tilde{z}$')




"""
Import instantaneous vector data and LIF data

-The vector components for the specified times were extracted from DaVis
-The LIF image corresponding to the specified times were also exported from DaVis
"""


instVecPath=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\CombinedVectorLIF"
instVecFiles=piv.get_PIV_files(instVecPath)
instLIFimages=piv.get_LIF_files(instVecPath)
labels=['Case_1_','Case_2_','Case_7_','Case_8_']
times=['t_0','t_1','t_2','t_3','t_4','t_5','t_6']
fields=['u','v']
keys=['B00001','B00002']

instVec=piv.load_PIV_dataframes(instVecFiles,labels,times,fields,keys)


tr_instVec={}
tr_instVec['Vector']=piv.treatMeanVecData(instVec,[0,25.46,0],[-110,110,0,80],5,avg=False,norm=False)






"""
Make the basic elements needed for Fig. 6 of the article (Instantaneous vector + LIF), the rest is done in Inkscape
"""


cases=['Case_1_','Case_2_','Case_7_','Case_8_']
times=['t_2','t_3','t_4','t_5']


fig, ax = plt.subplots(4,4,frameon=False,constrained_layout=False)

fig.set_figwidth(7.5,forward=True)
fig.set_figheight(4,forward=True)


for condition, case in enumerate(cases):
    

    for t, time in enumerate(times):
        print(case, time)
        
        u=tr_instVec['Vector'][case][time]['u_grid']
        v=tr_instVec['Vector'][case][time]['v_grid']      
        color=np.sqrt(np.multiply(u,u)+np.multiply(v,v))
        color=np.nan_to_num(color,False,0.0)
        print(np.amax(color))

        X=tr_instVec['Vector'][case][time]['X']
        Y=tr_instVec['Vector'][case][time]['Y']  
        
        for path in instLIFimages:
        
            if '%s\\%s' %(case,time) in path:
                im=plt.imread(path)
                
                if 'Case_1_' in path:
                    im = im[235:614, 19:1140, :]
                if 'Case_2_' in path:
                    im = im[235:614, 19:1114, :]
                if 'Case_7_' in path:
                    im = im[235:614, 24:1137, :]
                if 'Case_8_' in path:
                    im = im[235:614, 19:1114, :]                    

        ax[condition][t].imshow(im,extent=[-130,70,0,70],origin='upper')        

        import matplotlib.colors
        norm = matplotlib.colors.Normalize(vmin=0,vmax=0.06,clip=False)
        
        CS=ax[condition][t].quiver(X,Y,u,v,color,cmap='jet',linewidth=2,headwidth=5,norm=norm)
        #CS=ax[condition][t].streamplot(X,Y,u,v, density=1,color=color,cmap='seismic',linewidth=0.25,arrowsize=0.3)

        ax[condition][t].set_xlim([-100,70])            
        ax[condition][t].set_ylim([0,70])    
        #ax[condition][t].xaxis.set_major_locator(plt.NullLocator())
        #ax[condition][t].yaxis.set_major_locator(plt.NullLocator())

plt.tight_layout()
plt.subplots_adjust(wspace=0.01,hspace=0.001)

plt.colorbar(CS, ax=ax, shrink=0.3, format='%.2f', location='right', pad=0.01, label='$\mathbf{v}_{yw}$ (m/s)')

plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\InstantaneousLIF\NoText_colorbar.pdf', format='pdf', dpi=600)






"""
Hi_3 analysis and figure
"""

path=r"C:\Users\Jason\Dropbox\2021\LabConfluence\Data\PlanarPIV_at_4D"


files=piv.get_PIV_files(path)

"""Provide data labels and load data to pandas dataframes"""
fieldKeys=['B00001','B00002'] #u,w,TKE, Rvw,sw,vor
fieldNames=['u','v']
shot_labels=['Shot_1']

#These are the names of the individual folders in path, you might only need one
plane_labels=['Hi_3_']



#main dictionary holding untreated data (i.e. output from DaVis)
hi_3=piv.load_PIV_dataframes(files,plane_labels,shot_labels,fieldNames,fieldKeys)



"""
Interpolate vector data to a regular grid 
"""
hi_3_treated={}
hi_3_treated['Vector']=piv.treatMeanVecData(hi_3,[0,25.46,0],[-110,110,0,70],2.5,avg=True,norm=False)


"""
Figure of hi_3
"""

coords=[[0,0],[0,1]]
cases=['Hi_3_']

fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)

fig.set_figwidth(5.0,forward=True)
fig.set_figheight(4,forward=True)



for index, case in enumerate(cases):
    
    
    if case in hi_3_treated['Vector']:
        print(case)

        u=hi_3_treated['Vector'][case]['avg']['u_grid']/velocityScale
        v=hi_3_treated['Vector'][case]['avg']['v_grid']/velocityScale     
        color=np.sqrt(np.multiply(u,u)+np.multiply(v,v))

        X=hi_3_treated['Vector'][case]['avg']['X']
        Y=hi_3_treated['Vector'][case]['avg']['Y']  
        
        import matplotlib.colors
        norm = matplotlib.colors.Normalize(vmin=0,vmax=0.65,clip=False)
        
        levels = np.linspace(0.0, 0.65, 8)
        CS=ax.contourf(X,Y, color,cmap='bone_r',alpha=0.6, norm=norm,extend='neither',levels=levels) 
        
        ax.streamplot(X,Y,u,v, density=1.2,color='black',linewidth=0.5,arrowsize=0.3)
        ax.set_xlim([-100,100])            
        ax.set_ylim([0,70])    
        
    else:
        pass

plt.colorbar(CS, ax=ax, shrink=0.6, format='%.2f', location='bottom', pad=0.005, label='$\\tilde{U}_{vw}$')

commonSubPlotElements()

plt.savefig(r'C:\Users\Jason\Dropbox\Apps\Overleaf\JFM laboratory article\Figures\fig_streamlines.pdf', format='pdf', dpi=600)






















































"