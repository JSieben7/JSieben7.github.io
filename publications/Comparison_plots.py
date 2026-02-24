import numpy as np, matplotlib.pyplot as plt, scipy.optimize as opt
#cosmolopy.luminosityfunction as cplf, 
plt.rc('font', family='serif')
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('legend', framealpha=0.8)

#plt.style.use('/Users/Jennifer/Anaconda3/pkgs/matplotlib-2.1.0-py36h11b4b9c_0/Lib/site-packages/matplotlib/mpl-data/stylelib/mystyle.mplstyle')
#plt.style.use('/Users/jsieben/Box Sync/Lazy folder/mystyle.mplstyle')


#   If on laptop, switch to defined functions (without cplf.) and call data sets differently


def schechterL(luminosity, phiStar, alpha, LStar): 
      """Schechter luminosity function.""" 
      LOverLStar = (luminosity/LStar) 
      return (phiStar/LStar) * LOverLStar**alpha * np.exp(- LOverLStar)

def schechter_func(logL, alpha, phistar, logLstar):
    #log version
    return np.log(10.) * phistar * 10.**((alpha+1.)*(logL-logLstar)) * np.exp(-10.**(logL-logLstar))

def schechterM(magnitude, phiStar, alpha, MStar): 
      """Schechter luminosity function by magnitudes.""" 
      MStarMinM = 0.4 * (MStar - magnitude) 
      return (0.4 * np.log(10) * phiStar * 
              10.0**(MStarMinM * (alpha + 1.)) * np.exp(-10.**MStarMinM))

def DPL(M, phistar, Mstar, alpha, beta):
    return phistar/(10**(0.4*(alpha+1)*(M-Mstar))+10**(0.4*(beta+1)*(M-Mstar)))

#Designed to compare data between different sets, run in the main LF project file on box
 
Plot_Rcompare=False
Plot_KB_O3compare=False
Plot_Hacompare=False

#R band comparison

    #AHA
#ANGIE'S BINS FOR AHA
if Plot_Rcompare==True:
    ("plotting R band")
    AHAR=np.loadtxt("./AHA/Rband_test/AHA_LF_R.txt")
    dM=1.0
    low,up = 0,12 #formerly 7 when not testing
    sig=AHAR[low:up,2]
    lindat=np.linspace(-24,-14,50)
    fit=opt.curve_fit(schechterM,AHAR[low:up,0],AHAR[low:up,1],sigma=AHAR[low:up,2],p0=(5*10**-5,-1.2,-16.60))
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(abs(fit[1])))
    print ("Schecteraha: ps, a, ms")
    print (ps, a, ms)
    np.savetxt("Tablesforpaper/AHA_R_fits.txt",(ps,a,ms,dM),header="SchechterM: phi*, alpha, M*, dM")
    #print ("error:")
    #print (perr)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    fit = opt.curve_fit(DPL, AHAR[low:up,0], AHAR[low:up,1], sigma=AHAR[low:up,2], p0=(10e-3, -16, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    #DPLfit = DPL(lindat, phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    #ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    
    from matplotlib.ticker import AutoMinorLocator
    plt.figure(2)
    ax = plt.subplot(1, 1, 1)
    plt.yscale('log')
    plt.plot(AHAR[low:up,0],AHAR[low:up,1]*dM,'o',color = '#0ACF00',label='AHA R data')
    #plt.plot(Vblue[0:fitto,3],np.log10(Vblue[0:fitto,1]),':',color = '#0069bf')
    err=AHAR[low:up,2]
    yerr_lower = np.abs(AHAR[low:up,1]-AHAR[low:up,2]-AHAR[low:up,1])
    yerr_upper = np.abs(AHAR[low:up,1]+AHAR[low:up,2]-AHAR[low:up,1])
    plt.errorbar(AHAR[low:up,0],AHAR[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#2D9B27')
    #plt.plot(lindat,DPLfit*dM,'-.',color='#2c5d4a',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $M^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    #plt.plot(lindat,DPLfit*dM,'-.',color='#2c5d4a',label="Double Power Law")
    #plt.plot(lindat,shct*dM,'--',color='#078600',label=r"Schechter Function fit $log\ \phi^*$=%1.2f, $\alpha$=%1.4f, "
     #                                                    r"$M^*$=%1.2f" %(np.log10(ps), a, ms))
          #red
    """"""
    Vred=np.loadtxt("./Red/Rband/kissrR_2020.txt")
    dM=1.0
    low,up=0,6
    fitto=up
    sig=Vred[low:up,2]
    fit=opt.curve_fit(schechterM,Vred[low:up,3],Vred[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-20))
    np.set_printoptions(suppress=False, precision=4)
    lindat=np.linspace(-24,-14,50)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    #print ("Schecterred: ps, a, ms")
    #print (ps, a, ms)
    #print ("error:")
    #print (perr)
    fit = opt.curve_fit(DPL, Vred[low:up,0], Vred[low:up,1], sigma=Vred[low:up,2], p0=(10e-3, -19, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    #DPLfit = DPL(lindat, phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    #ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    ax = plt.subplot(1, 1, 1)
    #plt.plot(Vred[:,3],Vred[:,1],'^',color = '#FD0006',label='KISS Red data')
    #plt.plot(Vred[0:fitto,3],np.log10(Vred[0:fitto,1]),':',color = '#ff2610')
    err=Vred[:,2]
    yerr_lower = np.abs(Vred[:,1]-Vred[:,2]-Vred[:,1])
    yerr_upper = np.abs(Vred[:,1]+Vred[:,2]-Vred[:,1])
    plt.errorbar(Vred[:,3],Vred[:,1],yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#BE2F33')
    #plt.plot(lindat,DPLfit*dM,'-.',color='#b0331a',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, "
     #                                                              r"$M^*$=%1.2f, "
      #                                              r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
       # np.log10(ps),ms,a,beta))
    #plt.plot(lindat,DPLfit*dM,'-.',color='#b0331a',label="Double Power Law")
    #plt.plot(lindat,shct*dM,'--',color='#A40004',label=r"Schechter Function fit $log\ \phi^*$=%1.2f, $\alpha$=%1.4f, "
     #                                                    r"$M^*$=%1.2f" %(np.log10(ps), a, ms))
    np.savetxt("Tablesforpaper/KISSR_R_fits.txt",(ps,a,ms,dM),header="SchechterM: phi*, alpha, M*, dM")
    print( "**************")
    
    
          #blue
    
    Vblue=np.loadtxt("./Blue/Rband/kissbr_2020.txt")
    dM=1.0
    low,up=0,6
    #print (Vblue[low:up,0])
    plt.figure(2)
    sig=Vblue[low:up,2]
    lindat=np.linspace(-23,-14,50)
    fit=opt.curve_fit(schechterM,Vblue[low:up,3],Vblue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-17))
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(abs(fit[1])))
    #print ("Schecterblue: ps, a, ms")
    #print (ps, a, ms)
    #print ("error:")
    #print (perr)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    
    fit = opt.curve_fit(DPL, Vblue[low:up,0], Vblue[low:up,1], sigma=Vblue[low:up,2], p0=(10e-3, -17, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    #DPLfit = DPL(lindat, phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    #ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    from matplotlib.ticker import AutoMinorLocator
    ax = plt.subplot(1, 1, 1)
    plt.plot(Vblue[low:up,3],Vblue[low:up,1],'v',color = '#3714B0',label='KISS Blue data')
    #plt.plot(Vblue[0:fitto,3],np.log10(Vblue[0:fitto,1]),':',color = '#0069bf')
    err=Vblue[low:up,2]
    yerr_lower = np.abs(Vblue[low:up,1]-Vblue[low:up,2]-Vblue[low:up,1])
    yerr_upper = np.abs(Vblue[low:up,1]+Vblue[low:up,2]-Vblue[low:up,1])
    plt.errorbar(Vblue[low:up,3],Vblue[low:up,1],yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#402C84')
    #plt.plot(lindat,DPLfit*dM,'-.',color='#2b2177',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $M^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    #plt.plot(lindat,shct*dM,'--',color='#1F0772',label=r"Schechter Function fit $log\ \phi^*$=%1.2f, $\alpha$=%1.4f, "
     #                                                            r"$M^*$=%1.2f" %(np.log10(ps), a, ms))
    #plt.plot(lindat,shct*dM,'--',color='#00ccef',label="Schechter Function")
    np.savetxt("Tablesforpaper/KISSB_R_fits.txt",(ps,a,ms,dM),header="SchechterM: phi*, alpha, M*, dM")
    
    
    
    plt.xlabel("Absolute R band Magnitude")
    plt.xlim(-12,-24)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    plt.ylim(10**-5.5,10**-1)
    #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    #plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
    #plt.legend()
    #plt.savefig('ComparisonPlots/Rband_040920.eps', format='eps', dpi=1000)
    #plt.savefig('ComparisonPlots/Rband_040920.png', format='png', dpi=300)
    #plt.savefig('AHA/Rband_060920.eps', format='eps', dpi=1000)
    #plt.savefig('AHA/Rband_060920.png', format='png', dpi=300)
    #plt.show()
    plt.close()
    


#Compare Ha band blue and O3 band blue
      #Ha
if Plot_KB_O3compare==True:
    print("starting o3plot")
    Hablue=np.loadtxt("./Blue/Ha/kissb_LF_2020.txt")
    dM=1.0
    low,up=1,9
    #print (Hablue[low:up,0])
    sig=Hablue[low:up,2]
    fit=opt.curve_fit(schechterM,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
    np.set_printoptions(suppress=False, precision=4)
    lindat=np.linspace(-108,-97,50)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    #print ("Schecterha: ps, a, ms")
    #print (ps, a, ms)
    #print ("error:")
    #print (perr)
    from matplotlib.ticker import AutoMinorLocator
    plt.figure(3)
    ax = plt.subplot(1, 1, 1)
    plt.yscale('log')
    #plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,'v',color = '#3714B0',label=r'KISSB, [O III] selected, '
     #                                                                              r'H$\alpha$ luminosity')
    #plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,':',color = '#0069bf')
    err=Hablue[low:up,2]
    yerr_lower = np.abs(Hablue[low:up,1]*dM-Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
    yerr_upper = np.abs(Hablue[low:up,1]*dM+Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
    plt.errorbar(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#402C84')
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#1F0772',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
     #                                                     r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
    fit=opt.curve_fit(DPL,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    np.savetxt("Tablesforpaper/KISSB_HA_DPL_fits.txt",(ps, ms, a, beta, dM),header="DPL: phi*, M*, alpha, beta, dM")
    perr = np.sqrt(np.diag(abs(fit[1])))
    #print ("DPLha: ps, a, ms, beta")
    #print (ps, a, ms, beta)
    #print ("error:")
    #print (perr)
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'-.',color='#6949D7',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    
          #O3
    dM=1.0
    low,up=1,8
    sig=O3blue[low:up,2]
    fit=opt.curve_fit(schechterM,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
    lindat=np.linspace(-108,-97.5,50)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    #print ("Schectero3: ps, a, ms")
    #print (ps, a, ms)
    #print ("error:")
    #print (perr)
    ax = plt.subplot(1, 1, 1)
    #plt.plot(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,'s',color = '#BF008A',label='KISSB, [O III] selected, '
     #                                                                                     '[O III] luminosity')
    #plt.plot(O3blue[low:up,3]/-2.5,np.log10(O3blue[low:up,1]),':',color = '#0069bf')
    err=O3blue[low:up,2]
    yerr_lower = np.abs(O3blue[low:up,1]*dM-O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
    yerr_upper = np.abs(O3blue[low:up,1]*dM+O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
    plt.errorbar(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#8F2471')
    #plt.title("KISS NB blue Comparison, d > 35Mpc")
    #plt.xlabel(r"Log L(H$\alpha$ or [O III]) [$erg\ s^{-1}$]")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    plt.ylim(10**-6,10**-1)
    #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    #plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#DF38B1',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, "
     #                                                              r"$\alpha$=%1.4f, $log\ L^*$=%1.2f" %(np.log10(ps), a,
      #                                                                                                  ms))
    fit=opt.curve_fit(DPL,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    np.savetxt("Tablesforpaper/KISSB_O3_DPL_fits.txt",(ps, ms, a, beta, dM),header="DPL: phi*, M*, alpha, beta, dM")
    perr = np.sqrt(np.diag(abs(fit[1])))
    #print ("DPLo3: ps, a, ms, beta")
    #print (ps, a, ms, beta)
    #print ("error:")
    #print (perr)
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'-.',color='#7C005A',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    
    plt.xlim(39,43)
    #plt.legend()
    #plt.savefig('ComparisonPlots/Ha_v_O3_050120.eps', format='eps', dpi=1000)
    #plt.savefig('ComparisonPlots/Ha_v_O3_050120.png', format='png', dpi=300)
    #plt.show()
    plt.close()



#input x values in magnitude

#Compare HA across AHA, KISSR, and KISSB
#DOUBLE POWER LAW
#log axis

from matplotlib.ticker import AutoMinorLocator
plt.figure(3)
ax = plt.subplot(1, 1, 1)
plt.yscale('log')
      #AHA
if Plot_Hacompare==True:
    #print ("TEST")
    print("plotting Ha LF")
    #AHA=np.loadtxt("./AHA/AHAHa_weighted.txt")
    AHA=np.loadtxt("./AHA/AHA_LF_20.txt")
    dM=1.25 #normalisation, 1.0 if bins are every 0.4 dex in Lum or 1 in mag (coming from angie's script)
    #dM=.625 #coming from my script
    #print (AHA)
    low,up=6,18 #Indexing is currently not inclusive at upper end!!
    sig=AHA[low:up,2]
    #schecter fit
    fit = opt.curve_fit(schechter_func, AHA[low:up,0], AHA[low:up,1], sigma=AHA[low:up,2], p0=(-1,.0003,41))
    np.set_printoptions(suppress=False, precision=4)
    lindat = np.linspace(36, 43, 50)
    shct = schechter_func(lindat, alpha=fit[0][0], phistar=fit[0][1], logLstar=fit[0][2])
    a, ps, ms = fit[0][0], fit[0][1], fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    print ("Schecteraha: ps, a, ms")
    print (ps, a, ms)
    print ("error:")
    print (perr)
    #plt.plot(lindat,shct*dM,'--',color='#078600',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
     #                                                r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
    #plt.plot(lindat,shct,'--',color='#078600',label="Schechter Function")
    #DPL
    fit = opt.curve_fit(DPL, AHA[low:up,0]*-2.5, AHA[low:up,1], sigma=AHA[low:up,2], p0=(10e-3, -2.5*40, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    lindat = np.linspace(36, 43, 50)
    shct = DPL(lindat*-2.5, phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    perr = np.sqrt(np.diag(abs(fit[1])))
    print ("DPLaha: ps, a, ms, beta")
    print (ps, a, ms, beta)
    print ("error:")
    print (perr)
    np.savetxt("Tablesforpaper/AHA_HA_DPL_fits.txt",(ps, ms, a, beta, dM),header="DPL: phi*, M*, alpha, beta, dM")
    #plt.plot(lindat,shct*dM,'-.',color='#42E73A',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms/-2.5,a,beta))
    plt.plot(lindat,shct*dM,'-.',color='#42E73A',label="Double Power Law")
    #print (ps, ms, a, beta)
    #shct = powerlaw(AHA[low:up,0]*-2.5, alpha=21, k=-4)
    #a, k = 1e21, 15
    #plt.plot(lindat,10**shct, '-.', color='#b20ae8', label='weighted Schechter Function fit')
    ax = plt.subplot(1, 1, 1)
    #plt.plot(AHA[:,0],AHA[:,1]*dM,'o',color = '#C0F6AB')#,label='fit between 0 and 20')
    plt.plot(AHA[low:up,0],AHA[low:up,1]*dM,'o',color = '#0ACF00',label="AHA data")
    #plt.plot(AHA[:,0],AHA[:,1],':',color = '#134026')
    err=AHA[:,2]
    plt.errorbar(AHA[low:up,0],AHA[low:up,1]*dM,yerr=AHA[low:up,2],ecolor='#2D9B27',fmt='none')
    #plt.errorbar(AHA[:,0],AHA[:,1]*dM,yerr=AHA[:,2],ecolor='#2D9B27',fmt='none')
    #plt.title("Full Ha band Comparison")
    plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    #plt.ylim(10**-5.5,.1)AHA[low:up,0]
    #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.ylabel(r"Log $\phi$ [$1/Mpc^3$]")
    
             #red
    """"""
    Hared=np.loadtxt("./Red/Ha/kissr_2020.txt")
    #print (Hared)
    dM=2
    low,up=6,23
    sig=Hared[low:up,2]
    #print (Hared[24,3])
    #print (Hared[24,3]/-2.5)
    #schecter
    fit=opt.curve_fit(schechterM,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
    lindat=np.linspace(-108,-90,50)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    print ("Schecterred: ps, a, ms")
    print (ps, a, ms)
    print ("error:")
    print (perr)
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#A40004',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
     #                                                r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#A40004',label="Schechter Function")
    #DPL
    fit=opt.curve_fit(DPL,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    lindat=np.linspace(-108,-90,50)
    shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    perr = np.sqrt(np.diag(abs(fit[1])))
    print ("DPLred: ps, a, ms, beta")
    print (ps, a, ms, beta)
    print ("error:")
    print (perr)
    np.savetxt("Tablesforpaper/KISSR_HA_DPL_fits.txt",(ps, ms, a, beta, dM),header="DPL: phi*, M*, alpha, beta, dM")
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'-.',color='#FE3F44',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    plt.plot(lindat/-2.5,shct*dM,'-.',color='#FE3F44',label="Double Power Law")
    ax = plt.subplot(1, 1, 1)
    #plt.plot(Hared[:,3]/-2.5,Hared[:,1]*dM,'^',color = '#FE7276')
    plt.plot(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,'^',color = '#FD0006',label='KISS Red data')
    #plt.plot(Hared[low:up,3]/-2.5,Hared[low:up,1],':',color = '#ff2610')
    err=Hared[low:up,2]
    yerr_lower = np.abs(Hared[low:up,1]*dM-Hared[low:up,2]*dM-Hared[low:up,1]*dM)
    yerr_upper = np.abs(Hared[low:up,1]*dM+Hared[low:up,2]*dM-Hared[low:up,1]*dM)
    plt.errorbar(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#BE2F33')
    yerr_lower = np.abs(Hared[:,1]*dM-Hared[:,2]*dM-Hared[:,1]*dM)
    yerr_upper = np.abs(Hared[:,1]*dM+Hared[:,2]*dM-Hared[:,1]*dM)
    #plt.errorbar(Hared[:,3]/-2.5,Hared[:,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#BE2F33')
    #plt.title("KISS Ha band Comparison, d > 35Mpc")
    plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
    #ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    #plt.ylim(-5.5,-1.5)
    #ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    
         #blue
    """"""
    Hablue=np.loadtxt("./Blue/Ha/kissb_LF_2020.txt")
    dM = 1.0
    #print (Hablue)
    low,up=1,9
    sig=Hablue[low:up,2]
    #schecter
    fit=opt.curve_fit(schechterM,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
    lindat=np.linspace(-108,-90,50)
    shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
    ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
    perr = np.sqrt(np.diag(fit[1]))
    print ("Schecterblue: ps, a, ms")
    print (ps, a, ms)
    print ("error:")
    print (perr)
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#1F0772',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
     #                                                r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
    #plt.plot(lindat/-2.5,shct*dM,'--',color='#1F0772',label="Schechter Function")
    #DPL
    fit=opt.curve_fit(DPL,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
    np.set_printoptions(suppress=False, precision=4)
    lindat=np.linspace(-108,-90,50)
    shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
    ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
    perr = np.sqrt(np.diag(abs(fit[1])))
    print ("DPLblue: ps, a, ms, beta")
    print (ps, a, ms, beta)
    print ("error:")
    print (perr)
    ms=ms/-2.5
    #plt.plot(lindat/-2.5,shct*dM,'-.',color='#6949D7',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
     #                                               r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
      #  np.log10(ps),ms,a,beta))
    plt.plot(lindat/-2.5,shct*dM,'-.',color='#6949D7',label="Double Power Law")
    #plt.plot(Hablue[:,3]/-2.5,Hablue[:,1]*dM,'v',color = '#866FD7')#,label='fit between 1 and 8')
    plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,'v',color = '#3714B0',label='KISS Blue data')
    #plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1],':',color = '#0069bf')
    err=Hablue[low:up,2]*dM
    yerr_lower = np.abs(Hablue[low:up,1]*dM-Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
    yerr_upper = np.abs(Hablue[low:up,1]*dM+Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
    plt.errorbar(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#402C84')
    yerr_lower = np.abs(Hablue[:,1]*dM-Hablue[:,2]*dM-Hablue[:,1]*dM)
    yerr_upper = np.abs(Hablue[:,1]*dM+Hablue[:,2]*dM-Hablue[:,1]*dM)
    #plt.errorbar(Hablue[:,3]/-2.5,Hablue[:,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#402C84')
    plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
    plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
    
    """"""
    plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
    plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
    plt.legend()
    plt.xlim(39,44)
    plt.ylim(10**-5.5,10**-1)
    #plt.savefig('ComparisonPlots/Haband_040920.eps', format='eps', dpi=1000)
    #plt.savefig('ComparisonPlots/Haband_040920.png', format='png', dpi=300)
    
    #plt.savefig('AHA/AHA_Ha_040920.eps', format='eps', dpi=1000)
    #plt.savefig('AHA/AHA_Ha_040920.png', format='png', dpi=300)
    #plt.savefig('Red/Ha/KR_Ha_040920.eps', format='eps', dpi=1000)
    #plt.savefig('Red/Ha/KR_Ha_040920.png', format='png', dpi=300)
    #plt.savefig('Blue/Ha/KB_Ha_040920.eps', format='eps', dpi=1000)
    #plt.savefig('Blue/Ha/KB_Ha_040920.png', format='png', dpi=300)
    
    plt.show()
    plt.close()




#2019 KR Ha vs 2020 KR Ha
"""
from matplotlib.ticker import AutoMinorLocator
plt.figure(11)
ax = plt.subplot(1, 1, 1)
plt.yscale('log')
Hared=np.loadtxt("./Red/Ha/kissr_2020.txt")
#print (Hared)
dM=2
low,up=6,23
sig=Hared[low:up,2]
#schecter
fit=opt.curve_fit(schechterM,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-90,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecterred: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#fc927c',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
                                                 r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#fc927c',label="Schechter Function")
#DPL
fit=opt.curve_fit(DPL,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
lindat=np.linspace(-108,-90,50)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPLred: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#b0331a',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#b0331a',label="Double Power Law")
ax = plt.subplot(1, 1, 1)
plt.plot(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,'^',color = '#ff2610',label='KISS Red 2020 data')
err=Hared[low:up,2]
yerr_lower = np.abs(Hared[low:up,1]*dM-Hared[low:up,2]*dM-Hared[low:up,1]*dM)
yerr_upper = np.abs(Hared[low:up,1]*dM+Hared[low:up,2]*dM-Hared[low:up,1]*dM)
plt.errorbar(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#f25749')

Hared=np.loadtxt("./Red/Ha/kissr_2019.txt")
#print (Hared)
dM=2
low,up=6,23
sig=Hared[low:up,2]
#schecter
fit=opt.curve_fit(schechterM,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-90,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecterred: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
                                                 r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#fc927c',label="Schechter Function")
#DPL
fit=opt.curve_fit(DPL,Hared[low:up,3],Hared[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
lindat=np.linspace(-108,-90,50)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPLred: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#b0331a',label="Double Power Law")
ax = plt.subplot(1, 1, 1)
plt.plot(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,'v',color = 'black',label='KISS Red 2019 data')
err=Hared[low:up,2]
yerr_lower = np.abs(Hared[low:up,1]*dM-Hared[low:up,2]*dM-Hared[low:up,1]*dM)
yerr_upper = np.abs(Hared[low:up,1]*dM+Hared[low:up,2]*dM-Hared[low:up,1]*dM)
plt.errorbar(Hared[low:up,3]/-2.5,Hared[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='grey')

plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
plt.legend()
plt.xlim(39,44)
plt.ylim(10**-5.5,10**-1)
#plt.savefig('ComparisonPlots/KRHA_2.11.20.eps', format='eps', dpi=1000)
#plt.savefig('ComparisonPlots/KRHA_2.11.20.png', format='png', dpi=300)
#plt.show()
#plt.close()
"""

#2019 KB Ha vs 2020 KB Ha
"""
from matplotlib.ticker import AutoMinorLocator
plt.figure(12)
ax = plt.subplot(1, 1, 1)
plt.yscale('log')
Hablue=np.loadtxt("./Blue/Ha/kissb_LF_2020.txt")
dM = 1.0
#print (Hablue)
low,up=1,9
sig=Hablue[low:up,2]
#schecter
fit=opt.curve_fit(schechterM,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-90,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecterblue: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#00ccef',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
                                                 r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#00ccef',label="Schechter Function")
#DPL
fit=opt.curve_fit(DPL,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
lindat=np.linspace(-108,-90,50)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPLblue: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#2b2177',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))
#plt.plot(lindat/-2.5,shct*dM,'-.',color='#2b2177',label="Double Power Law")
plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,'v',color = '#0069bf',label='KISS Blue 2020 data')
err=Hablue[low:up,2]*dM
yerr_lower = np.abs(Hablue[low:up,1]*dM-Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
yerr_upper = np.abs(Hablue[low:up,1]*dM+Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
plt.errorbar(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#00a8e0')


Hablue=np.loadtxt("./Blue/Ha/kissb_LF_2019.txt")
dM = 1.0
#print (Hablue)
low,up=1,9
sig=Hablue[low:up,2]
#schecter
fit=opt.curve_fit(schechterM,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-90,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecterblue: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, $\alpha$=%1.4f, "
                                                 r"$log\ L^*$=%1.2f" %(np.log10(ps), a, ms))
#DPL
fit=opt.curve_fit(DPL,Hablue[low:up,3],Hablue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
lindat=np.linspace(-108,-90,50)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPLblue: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))
plt.plot(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,'^',color = '#0069bf',label='KISS Blue 2019 data')
err=Hablue[low:up,2]*dM
yerr_lower = np.abs(Hablue[low:up,1]*dM-Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
yerr_upper = np.abs(Hablue[low:up,1]*dM+Hablue[low:up,2]*dM-Hablue[low:up,1]*dM)
plt.errorbar(Hablue[low:up,3]/-2.5,Hablue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='grey')
plt.xlabel(r"Log L(H$\alpha$) [$erg\ s^{-1}$]")
plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
plt.legend()
plt.xlim(39,44)
plt.ylim(10**-5.5,10**-1)
#plt.savefig('ComparisonPlots/KBHA_2.11.20.eps', format='eps', dpi=1000)
#plt.savefig('ComparisonPlots/KBHA_2.11.20.png', format='png', dpi=300)
#plt.show()
#plt.close()
"""

#2019 KB O3 vs 2020 KB O3
"""
from matplotlib.ticker import AutoMinorLocator
plt.figure(13)
ax = plt.subplot(1, 1, 1)
plt.yscale('log')
O3blue=np.loadtxt("./Blue/O3/kissb_LF_2020.txt")
dM=1.0
low,up=1,8
sig=O3blue[low:up,2]
fit=opt.curve_fit(schechterM,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-97.5,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecter: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ax = plt.subplot(1, 1, 1)
plt.plot(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,'s',color = '#2eb516',label='KISSB 2020, [O III] selected, '
                                                                                      '[O III] luminosity')
err=O3blue[low:up,2]
yerr_lower = np.abs(O3blue[low:up,1]*dM-O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
yerr_upper = np.abs(O3blue[low:up,1]*dM+O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
plt.errorbar(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='#23De94')
plt.xlabel(r"Log L(H$\alpha$ or [O III]) [$erg\ s^{-1}$]")
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
plt.ylim(10**-6,10**-1)
plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#68d46d',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, "
                                                               r"$\alpha$=%1.4f, $log\ L^*$=%1.2f" %(np.log10(ps), a,
                                                                                                    ms))
fit=opt.curve_fit(DPL,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPL: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='#46ab46',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))

O3blue=np.loadtxt("./Blue/O3/kissb_LF_2019.txt")
dM=1.0
low,up=1,8
sig=O3blue[low:up,2]
fit=opt.curve_fit(schechterM,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(5*10**-3,-1.2,-104))
lindat=np.linspace(-108,-97.5,50)
shct=schechterM(lindat,phiStar=fit[0][0],alpha=fit[0][1],MStar=fit[0][2])
ps,a,ms=fit[0][0], fit[0][1],fit[0][2]
perr = np.sqrt(np.diag(fit[1]))
print ("Schecter: ps, a, ms")
print (ps, a, ms)
print ("error:")
print (perr)
ax = plt.subplot(1, 1, 1)
plt.plot(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,'d',color = 'black',label='KISSB 2019, [O III] selected, '
                                                                                      '[O III] luminosity')
err=O3blue[low:up,2]
yerr_lower = np.abs(O3blue[low:up,1]*dM-O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
yerr_upper = np.abs(O3blue[low:up,1]*dM+O3blue[low:up,2]*dM-O3blue[low:up,1]*dM)
plt.errorbar(O3blue[low:up,3]/-2.5,O3blue[low:up,1]*dM,yerr=[yerr_lower,yerr_upper], fmt='none',ecolor='grey')
plt.xlabel(r"Log L(H$\alpha$ or [O III]) [$erg\ s^{-1}$]")
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
plt.ylim(10**-6,10**-1)
plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Schechter Function fit $log\ \phi^*$=%1.4f, "
                                                               r"$\alpha$=%1.4f, $log\ L^*$=%1.2f" %(np.log10(ps), a,
                                                                                                    ms))
fit=opt.curve_fit(DPL,O3blue[low:up,3],O3blue[low:up,1],sigma=sig,p0=(10e-3, -2.5*40, -1.7, -3.8))
np.set_printoptions(suppress=False, precision=4)
shct=DPL(lindat,phistar=fit[0][0], Mstar=fit[0][1], alpha=fit[0][2], beta=fit[0][3])
ps, ms, a, beta=fit[0][0], fit[0][1], fit[0][2], fit[0][3]
perr = np.sqrt(np.diag(abs(fit[1])))
print ("DPL: ps, a, ms, beta")
print (ps, a, ms, beta)
print ("error:")
print (perr)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct*dM,'-.',color='black',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))

plt.xlim(39,43)
plt.legend()
#plt.savefig('ComparisonPlots/KBO3_2.13.20.eps', format='eps', dpi=1000)
#plt.savefig('ComparisonPlots/KBO3_2.13.20.png', format='png', dpi=300)
#plt.show()
plt.close()
"""


#DPL Sanity check
"""
print("checking sanity")
from matplotlib.ticker import AutoMinorLocator
plt.figure(11)
ax = plt.subplot(1, 1, 1)
plt.yscale('log')
lindat=np.linspace(-108,-97.5,50)
ps, ms, a, beta=2**-3, -99.61, 0.2099, -1.91
shct=DPL(lindat,phistar=ps, Mstar=ms, alpha=a, beta=beta)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct,'-.',color='red',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))

ps, ms, a, beta=2**-3, -99.61, 0.5, -1.91
shct=DPL(lindat,phistar=ps, Mstar=ms, alpha=a, beta=beta)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct,'-.',color='blue',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))

ps, ms, a, beta=2**-3, -99.61, 0.2099, -3
shct=DPL(lindat,phistar=ps, Mstar=ms, alpha=a, beta=beta)
ms=ms/-2.5
plt.plot(lindat/-2.5,shct,'-.',color='green',label=r"Double Power Law fit $log\ \phi^*$=%1.4f, $log\ L^*$=%1.2f, "
                                                r"$\alpha$=%1.2f, $\beta$=%1.2f" %(
    np.log10(ps),ms,a,beta))


plt.xlabel(r"Log L(H$\alpha$ or [O III]) [$erg\ s^{-1}$]")
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
plt.ylim(10**-6,10**-1)
#ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.ylabel(r"Log $\phi$ [$Mpc^{-3}$]")

plt.xlim(39,43)
plt.legend()
#plt.savefig('ComparisonPlots/sanitycheck.eps', format='eps', dpi=1000)
#plt.savefig('ComparisonPlots/sanitycheck.png', format='png', dpi=300)
#plt.show()
plt.close()

"""
