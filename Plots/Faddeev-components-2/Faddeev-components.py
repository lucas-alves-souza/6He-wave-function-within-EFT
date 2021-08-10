import numpy as np
import matplotlib.pyplot as plt  
from matplotlib.ticker import AutoMinorLocator

from scipy import interpolate
from matplotlib import rc
plt.style.use('souza3') 


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
 


def plotar(fname)  :
   fig, (ax) = plt.subplots(figsize=(7,5.5))

#  plt.subplots_adjust(wspace = .0001) 

   name = '../../Faddeev-components-Lamb='
   fig.suptitle(r'', fontsize=17)

   cut = ['0300' ,'0800' ,'3000','8000']
  
   Fc=[]
   Fn=[] 

   cor=['g','orange','r','b','g', 'k']
   typ=[':','-','-.','--','-',':']
   i=0   
   for j in cut:
    
    xfilea = name+j 
    DATA=np.loadtxt('%s.dat'%xfilea)
    q=DATA[:,0];Fn=DATA[:,1];Fc=DATA[:,2]   
    fintc = interpolate.interp1d(q, Fc)   
    fintn = interpolate.interp1d(q, Fn) 
    #Fc=[np.sqrt(x) for x in Fc]
    x=50. 
    if (fname == 'Fc'):
      ax.plot(q,q* Fc/x/fintc(x) ,'%s'%typ[i],linewidth=3, color='%s'%cor[i], label=r'$\Lambda=$%s MeV'%int(j))
      ax.set_ylim([-1.71,1.61]) 
      ax.set_ylabel(r'$q$ $%s(q) / x%s(x)$'%('F_c','F_c'))
      ff='F_c'
    if (fname == 'Fn'):
      ax.plot( (q),np.sqrt(q)* Fn/ np.sqrt(x)/fintn(x) ,'%s'%typ[i],linewidth=3, color='%s'%cor[i], label=r'$\Lambda=$%s MeV'%int(j))
      ax.set_ylim([-4.2,4])
      ax.set_ylabel(r'$\sqrt{q}$ $%s(q) / \sqrt{x}%s(x)$'%('F_n','F_n'))
      ff='F_n'
    i=i+1
    
   ax.legend(loc='lower left',  frameon=False,fontsize=14)
  
    
   plt.semilogx() 
   ax.set_xlim([.4,10000])
   ax.set_xlabel(r'$q$ (MeV)')
  
   ax.tick_params(direction="inout", length=6, width=1, color="k") 
   ax.text(0.13, 0.8,'$x =$ %s MeV'%x,  transform = ax.transAxes)
   #ax.text(0.25, 0.25,r'$F_n(q) $',  transform = ax.transAxes)
  
   fig.subplots_adjust(top=.93)
   #plt.tight_layout() 

   plt.savefig('./%s.pdf'%fname,bbox_inches='tight')

F=['Fc' , 'Fn']

for  i in F :
 plotar(i) 