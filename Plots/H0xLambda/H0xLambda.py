import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import AutoMinorLocator
plt.style.use('souza')
from matplotlib import rc
plt.rcParams.update({'font.size': 20})
plt.rcParams['patch.linewidth']=0.7 

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
 
fig, (ax) = plt.subplots(figsize=(7,5.5))

#plt.subplots_adjust(wspace = .0001) 

name = 'H0xLambda'
fig.suptitle(r'', fontsize=17)



xfilea ='H0xLambda' 
DATA=np.loadtxt('./%s.dat'%xfilea)
Lambda=DATA[:,0];B3=DATA[:,1]   
   
ax.set_xlabel(r'$ \Lambda$ (MeV)') 
ax.set_ylabel(r'H${}_0 $ ($\Lambda$)')
 


ax.plot(Lambda,B3,'--',linewidth=1.9, color='g', label=r'')
 
ax.set_ylim([-40.,40])
ax.set_xlim([50.,10000])  
#ax.set_xticks(np.arange(10, 30, 5))
 
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(direction="inout", length=6, width=1, color="k")
#ax.tick_params(direction="inout", length=6, width=1, color="k")

ax.semilogx()

#for ax in fig.get_axes():
#    ax.label_outer()


fig.subplots_adjust(top=.93)
#plt.tight_layout()
plt.savefig('./%s.pdf'%name,bbox_inches='tight')
