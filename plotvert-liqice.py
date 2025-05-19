import numpy as np
import matplotlib
import matplotlib.pyplot as plt

outdir = './outdir/'
infile = outdir + 'parceltime-liqice.dat'

tcolumns = ['tc','zalt','temp','press', 'qvap','qliq','qice','RH', 'Ntot',
 'Nctot', 'ravg', 'rcavg', 'Ntoti', 'ravgi','aavg','cavg']
datas = np.genfromtxt(fname = infile, dtype = 'float', names=tcolumns)
qtot = datas['qvap'] + datas['qliq'] + datas['qice']

infile = outdir + 'initialdrydist.dat'
tcolumns = ['re','Nl'  , 'rlo']
drydists = np.genfromtxt(fname = infile, dtype = 'float', names=tcolumns)

infile = outdir + 'initialwetdist.dat'
tcolumns = ['Nl','r', 'msalt']
wetdists = np.genfromtxt(fname = infile, dtype = 'float', names=tcolumns)

fig1 = plt.figure(figsize = (12, 7))
ax1 = fig1.add_subplot(2, 3, 1)
ax1.plot(datas['tc'],datas['temp'],linestyle = '-', color = 'k')
ax1.set_xlabel('time [s]')
ax1.set_ylabel('temperature [C]')

ax2 = fig1.add_subplot(2, 3, 2)
ax2.plot(datas['tc'],datas['press'],linestyle = '-', color = 'k')
ax2.set_xlabel('time [s]')
ax2.set_ylabel('pressure [mb]')

ax3 = fig1.add_subplot(2, 3, 3)
ax3.plot(datas['tc'],datas['RH'],linestyle = '-', color = 'k')
ax3.set_xlabel('time [s]')
ax3.set_ylabel('RH [%]')

ax4 = fig1.add_subplot(2, 3, 4)
ax4.plot(datas['tc'],datas['qvap'],linestyle = '-', color = 'k',label='vapor')
ax4.plot(datas['tc'],datas['qliq'],linestyle = '--', color = 'k',label='liquid')
ax4.plot(datas['tc'],datas['qice'],linestyle = '--', color = 'magenta',label='ice')
ax4.set_xlabel('time [s]')
ax4.set_ylabel('mixing ratio [g/kg]')
ax4.legend()

ax5 = fig1.add_subplot(2, 3, 5)
ax5.plot(datas['tc'],datas['Nctot'],linestyle = '-', color = 'k')
ax5.set_xlabel('time [s]')
ax5.set_ylabel(r'concentration [# cm$^{-3}$]')

ax6 = fig1.add_subplot(2, 3, 6)
ax6.plot(datas['tc'],datas['rcavg'],linestyle = '-', color = 'k',label='liquid')
ax6.plot(datas['tc'],datas['ravgi'],linestyle = '--', color = 'magenta',label='ice')
ax6.plot(datas['tc'],datas['aavg'],linestyle = '-', color = 'red',linewidth=1,label='ice,a')
ax6.plot(datas['tc'],datas['cavg'],linestyle = '-', color = 'blue',linewidth=1,label='ice,c')
ax6.set_xlabel('time [s]')
ax6.set_ylabel(r'average radius [$\mu$m]')
ax6.legend()

fig2 = plt.figure(figsize = (14, 7))
ax1 = fig2.add_subplot(2, 3, 1)
ax1.plot(datas['temp'],datas['zalt'],linestyle = '-', color = 'k')
ax1.set_ylabel('height [m]')
ax1.set_xlabel('temperature [C]')

ax2 = fig2.add_subplot(2, 3, 2)
ax2.plot(datas['press'],datas['zalt'],linestyle = '-', color = 'k')
ax2.set_ylabel('height [m]')
ax2.set_xlabel('pressure [mb]')

ax3 = fig2.add_subplot(2, 3, 3)
ax3.plot(datas['RH'],datas['zalt'],linestyle = '-', color = 'k')
ax3.set_ylabel('height [m]')
ax3.set_xlabel('RH [%]')

ax4 = fig2.add_subplot(2, 3, 4)
ax4.plot(datas['qvap'],datas['zalt'],linestyle = '-', color = 'k',label='vapor')
ax4.plot(datas['qliq'],datas['zalt'],linestyle = '--', color = 'k',label='liquid')
ax4.plot(datas['qice'],datas['zalt'],linestyle = '--', color = 'firebrick',label='ice')
ax4.plot(qtot,datas['zalt'],linestyle = '--', color = 'purple',label='total')
ax4.set_ylabel('height [m]')
ax4.set_xlabel('mixing ratio [g/kg]')
ax4.legend()

ax5 = fig2.add_subplot(2, 3, 5)
ax5.plot(datas['Nctot'],datas['zalt'],linestyle = '-', color = 'k',label='cloud')
ax5.plot(datas['Ntot'],datas['zalt'],linestyle = '--', color = 'k',label='CCN + cloud')
ax5.plot(datas['Ntoti'],datas['zalt'],linestyle = '--', color = 'firebrick',label='ice')
ax5.set_xscale('log')
ax5.set_xlim(0.001,500.)
ax5.set_ylabel('height [m]')
ax5.set_xlabel(r'concentration [# cm$^{-3}$]')
ax5.legend()

ax6 = fig2.add_subplot(2, 3, 6)
ax6.plot(datas['rcavg'],datas['zalt'],linestyle = '-', color = 'k',label='cloud')
ax6.plot(datas['ravg'],datas['zalt'],linestyle = '--', color = 'k',label='total')
ax6.plot(datas['ravgi'],datas['zalt'],linestyle = '--', color = 'magenta',label='ice')
ax6.plot(datas['aavg'],datas['zalt'],linestyle = '-', color = 'red',linewidth=1,label='ice,a')
ax6.plot(datas['cavg'],datas['zalt'],linestyle = '-', color = 'blue',linewidth=1,label='ice,c')
ax6.set_ylabel('height [m]')
ax6.set_xlabel(r'average radius [$\mu$m]')
ax6.set_xscale('log')
ax6.legend()

plt.savefig('vertplot-liqice.pdf')

fig4 = plt.figure(figsize = (9, 7))
ax1 = fig4.add_subplot(2, 1, 1)
ax1.plot(drydists['rlo']*1.e6,drydists['Nl'],linestyle = '-', color = 'k')
ax1.plot(wetdists['r']*1.e6,wetdists['Nl'],linestyle = '--', color = 'firebrick')
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.set_xlabel(r'radius [\$mu$m]')
ax1.set_ylabel('concentration density')

ax2 = fig4.add_subplot(2, 1, 2)
ax2.plot(wetdists['r']*1.e6,wetdists['Nl'],linestyle = '-', color = 'k')
ax2.set_xscale('log')
ax2.set_xlabel(r'radius [\$mu$m]')
ax2.set_ylabel('concentration density')

plt.show()
