import numpy as np
import matplotlib.pyplot as plt  
from matplotlib.ticker import AutoMinorLocator

from matplotlib import rc
#plt.style.use('seaborn-notebook')

SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['lines.dashed_pattern'] = [10, 5]
plt.rcParams['lines.dashdot_pattern'] = [10, 3, 3,  3]
plt.rcParams['lines.dotted_pattern'] = [3, 3]
plt.rcParams['lines.scale_dashes'] = False 

#print(plt.style.available) 


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)
 

name = '../../Faddeev-components-Lamb=' 
cut = ['0200','0400' ,'0800','1200','1600'] 
 
F=['c','n']
  
x=0


cor=['gray','orange','b','r','g', 'k']
typ=[':','-','-.','--','-',':']
for indF in F: 
   fig, (ax) = plt.subplots(figsize=(7,5.5))
   #fig.suptitle(r'', fontsize=17)
   i=0    

   for j in cut:
    
    xfilea = name+j 
    DATA=np.loadtxt('%s.dat'%xfilea)
    q=DATA[:,0]
    


    if (indF=='n'): 
     ax.plot(q,DATA[:,1],'%s'%typ[i],linewidth=3, color='%s'%cor[i], label=r'$\Lambda=$%s MeV'%int(j))
    else: 
     ax.plot(q,DATA[:,2],'%s'%typ[i],linewidth=3, color='%s'%cor[i], label=r'$\Lambda=$%s MeV'%int(j))
    
    i=i+1
   
   ax.legend(frameon=False)  
   ax.set_xlim([0.,400]) 
   ax.set_xlabel(r'$q$ (MeV)') 
   ax.set_ylabel(r'$F_%s(q) / F_c(x)$'%indF)
   ax.text(0.2, 0.8,'$x =$ %s '%x,  transform = ax.transAxes)
    # ax.xaxis.set_minor_locator(AutoMinorLocator())
   #ax.tick_params(direction="inout", length=6, width=1, color="k") 
   plt.savefig('Faddeev-components-F%s.pdf'%indF,bbox_inches='tight')
 
    
