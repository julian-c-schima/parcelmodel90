program parcelmodel

  implicit none
  
  ! stuff for vode
  integer neq,neqtot,iequ,ic,nbliq,nbice
  integer itol,itask,istate,iopt,lrw,liw,mf
  real tc,tout,t
  parameter (neq=5,nbliq=100,nbice=5000,neqtot=nbliq+nbice+neq)
  integer ipar(20*neqtot)
  real y(neqtot),dy(neqtot),iwork(40),rwork(22+16*neqtot)  ! y = radius,  dy = dr/dt growth rate
  ! except for y(neqtot-neq:neqtot)
  real rpar(20*neqtot)
  real atol(neqtot)
  real rtol
  external jac,FunDiff,vode

  ! index/tolerance variables
  ! Note: don't be silly and accidentally define itmax as a real in some function.
  integer k,ivap,ipress,itemp,iliq,iice,nt
  integer iwvel,isolute,i,j,iicemodel,itmax,ihollow
  real para, parb

  ! liquid variables
  real radiusin,Nlo,Nui,ravgl,sigmaccn
  real rlo(nbliq), Nmrl(nbliq),rhodrop(nbliq),Nl(nbliq),msalt(nbliq)
  real rdryaer(nbliq),rhosolute,molmass,vhoff
  real nlinit(nbliq),nhetfrz(nbliq),rdrp(nbliq),rdrp0(nbliq)
  real Nltot,NUl,lwc,ml,mrl,dmrl,es,rdry,msalt2,tempk,Sl
  real Nctot,ravg,rcavg,Ntot

  ! ice variables
  integer nbicelim,hetmodel
  real alpha_a(nbice),alpha_c(nbice),Nmri(nbice),Ni(nbice)
  real ao(nbice),co(nbice),cinit(nbice),ainit(nbice),a(nbice),c(nbice),aavg,cavg
  real rhoxtal(nbice),rhodep(nbice),phi(nbice),rio(nbice),xhol(nbice)
  real IWC,MRi,nfr(nbice),iwcnuc,nhetin,rnhetin,nacmr,thet
  real sui,suiratio,naer,naermr,MRnucdep,MRnucfreez,areain
  real nhetlo,nhethi,nactmr,ei,qeqi,qeql,m,Ntoti,ravgi
  real Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin,alphaSphere
  real avp, athp, cvp, cthp, cap, phio, Ab, Absol, Ap, IGR, Vo, V
  real Rk, Rt, Fmax, Fa, Fc, drho, negsui, xholinit, xholavg
  ! extra rosette-related variables
  real scrit_a_temp, scrit_c_temp, scritrat, adj_fraction, adj_factor, hpyr, pyrang, xfac, bsize
  integer habit(nbice), ihabitchoice
  
  ! other variables / constants / functions to get variables
  real Tinit,TinitK,Pinit,ztopinit,RHinit
  real pi,Rd,Rv,cp,grav,dt,To,w,rhoa,rhoi
  real press,wmax,zalt,vapmr,nosc,tau_osc
  real qlsat,RH,FLHS
  real esi_new,es_new,funslocal,capacitance,hollowfrac
  real errabs, rel1
  
  ! thermodynamic variables
  real Dv, Kt, v_v, drhovidT, Ls

  logical endloop,firsthomog

  ! Common blocks
  COMMON/neqtotvars/itemp,ipress,ivap,iliq,iice,nbicelim ! Used to attain correct indices in y
  COMMON/envvars/w,rhodrop, Nl,rdrp,rdrp0,rdryaer,Nmrl,Nmri 
  COMMON/ccnvars/msalt
  COMMON/aerosolzzbren/Sl, tempk, es, rdry, msalt2
  common/ccnchar/vhoff,molmass,rhosolute  
  COMMON/inputvars/iicemodel
  ! Arrays for information pertaining to the a and c dimensions of ice
  COMMON/icevars/ao,co,a,c,rio,alpha_a,alpha_c,habit,xhol,rhoxtal,rhodep
  ! Input variables into root finding method for determining alpha
  COMMON/edgezzbren/sui,Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin
  
  ! iwval: pattern of vertical motion
  parameter(iwvel=1)      ! 1 = cosine w profile
                          ! 2 = constant w
                          
  ! isolute: aerosol comprising CCN                        
  parameter(isolute=2)    ! 1 = NaCl (sodium chloride)
                          ! 2 = (NH4)2SO4 (ammonium sulfate)
                          
  ! hetmodel: method for computing heterogenous freezing
  parameter(hetmodel=1)   ! 1 = instantaneous heterogeneous freezing
                          ! 2 = Demott (2010) PNAS
                          ! 3 = Knopf aw immersion freezing (2013)
                          ! 4 = Classical Theory (KC, 2004, 2009)
                          ! 5 = No Heterogeneous Freezing
                          
  ! ihabitchoice: habit class for ice crystals				  
  ihabitchoice=6          ! 1 = single crystals
  					      ! 2 = rosettes (random number of arms)
  						  ! 3 = combination of single crystals and rosettes (random num of arms)
  						  ! 4 - 12 = n-arm rosettes
                          
  ! iicemodel: method for computing ice growth
  iicemodel=2			  ! 0 = Basic spherical model (no effective diffusivity)
                          ! 1 = Spheres (can choose constant deposition coeff alphaSphere)
						  ! 2 = Faceted growth with spheroid approximation to adjust phi
						  ! 3 = Faceted growth with adjusting a and c each full timestep
						  ! 4 = Faceted growth with recalculating dimensions within VODE
  ! ihollow: allow ice crystals to hollow 
  ihollow=2				  ! 0 = No hollowing
  						  ! 1 = User-selected constant hollowing fraction
  						  ! 2 = Parameterization for hollowing prediction (function hollowfrac)
  ! Size at which rosettes begin to branch (set scrits to be the same for smaller crystals)
  bsize=1.e-5


  ! output file
  open(20,file='outdir/parceltime-liqice.dat')
  
  ! physical constants and initial conditions
  pi = 3.14159             ! key lime is the best flavor
  Rd = 287.                ! dry air gas constant
  Rv = 461.5    		   ! water vapor gas constant
  cp = 1004.5			   ! specific heat capacity of air at constant pressure
  grav = 9.81			   ! acceleration due to gravity
  To = 273.15              ! Reference temperature (0C)
  rhoi = 917.              ! density of solid ice
  
  ! initial conditions / model parameters
  dt = 1.0                 ! parcel model time step
  Ma = 10.                 ! growth mode for prism (a) face
  Mc = 10.                 ! growth mode for basal (c) face
  alphaSphere = 1          ! constant deposition coefficient if using spherical model
  xholinit = 0.7	       ! fraction of basal face width that is hollow (if ihollow = 1)
  Tinit = -30.             ! initial temperature [C]
  TinitK = Tinit + To
  Pinit = 40000.0          ! initial pressure [Pa]
  ztopinit = 2000.0        ! Total cloud depth [m]
  rhinit = 0.95             ! Initial relative humidity with respect to liquid
  wmax = 0.5               ! maximum vertical motion
  endloop = .false.        ! determine when to stop model based on vertical motion input

  !     initial conditions for microphysics
  radiusin = 1.e-6         ! radius of heterogeneous nuclei
  Nlo = 200.*100.**3       ! total # of dry CCN per volume[=]#/m^3 
  NUi = 4.                 ! determines ice distribution shape
  ravgi = 10.e-6           ! initial average radius of ice particles (m)
  ravgl = 0.02e-6          ! geometric radius of dry CCN distribution (m)
  sigmaccn = 2.2           ! geometric standard deviation of CCN distribution

  ! initialization for ice nucleation routines
  
  nbicelim = 0             ! start of with no ice bins
  areain = 1.e-8           ! cm^2
  nactmr = 0.0

  ! heterogeneous nucleation for different models
  ! typical high and low values are given
  if(hetmodel.eq.1)then   ! instantaneous concentration per liter
     nhetlo = 1.e-4         
     nhethi = 400.         
     nhetin = 1.
  endif
  if(hetmodel.eq.2)then
     nhetlo = 0.01          ! Demott (2015), concentration of aerosol per cc 
     nhethi = 400.
     nhetin = 100.
  endif
  if(hetmodel.eq.4)then   ! both Knopf and classical immersion theory
     nhetlo = 1.e-4          ! concentration of aerosol per cc
     nhethi = 200.
     nhetin = 50.
  endif
  naer = nhetin/1000.
  rnhetin = nhetin/rhoa
  thet = tempk + 1.1


  ! parameters for aerosol
  if(isolute.eq.1)then
     vhoff = 2.0            !van't hoff factor
     Molmass = 0.058443     ! molecular mass of NaCl
     rhosolute = 2165.      ! solute density (kg/m^3)
  else
     vhoff = 2.4            !van't hoff factor ! Modify Later
     Molmass = 132.13/1000. ! molecular mass of ammonium sulfate [(NH4)2SO4]
     rhosolute = 1769.
  endif

  ! oscillation scales
  tau_osc = ztopinit/(2.*wmax)               ! Parcel oscillation time-scale
  nosc = 1.                                  ! Number of 1/2 parcel oscilations (1/2 is up, second is down)
  nt = int(nosc*tau_osc/dt*pi)               ! Total number of time-steps
  if(iwvel.eq.2)then
     nt = int(nosc*ztopinit/wmax/dt)         ! if using a constant vertical motion
  endif
  print*,"Model initial conditions"
  print*,"T = ", Tinit, " C"
  print*,"pinit = ",  pinit/100., " hPa" 
  print*,"rhinit = ", + rhinit * 100., " percent"
  print*,"cloud depth = ", + ztopinit, " m"
  print*,"max w = ", wmax, " m/s"
  print*,"Elapsed time in model is ", int(nt * dt), "seconds"

  ! "pointers" to environmental variables
  itemp = nbliq + nbice + 1
  ipress = nbliq + nbice + 2
  ivap =  nbliq + nbice + 3
  iliq =  nbliq + nbice + 4
  iice =  nbliq + nbice + 5
  
  !     Set up VODE tolerances
  !     Note: Reduce tolerances if oscillations or jumps appear       
  rtol=1.e-7
  DO iequ=1,neqtot,1
     atol(iequ)=1.e-18
  ENDDO

  !     Set constants for solution of df(t)/dt in VODE
  itol   = 2
  itask  = 1
  istate = 1
  iopt   = 0
  lrw    = 22+10*neqtot
  liw    = 30
  mf     = 23
  
  ! initialize everything in y array (input to vode) as zeros
  
  ! environmental derivatives
  dy(itemp) = 0.0
  dy(ipress) = 0.0
  dy(ivap) = 0.0
  dy(iliq) = 0.0
  dy(iice) = 0.0
  
  DO i = 1, nbliq   ! liquid distribution
     y(i) = 0.0
     dy(i) = 0.0
     Nmrl(i) = 0.0
     Nl(i) = 0.0
     rlo(i) = 0.0
     msalt(i) = 0.0
     nhetfrz(i) = 0.0
  END DO

  DO i = 1, nbice   ! ice distribution
     j = i + nbliq  ! index for y array starts at nbliq + 1
     y(j) = 0.0
     dy(j) = 0.0
     alpha_a(i) = 0.0
     alpha_c(i) = 0.0
     Nmri(i) = 0.0
     Ni(i) = 0.0
     ao(i) = 0.0
     co(i) = 0.0
     rio(i) = 0.0
	 xhol(i) = 0.0
     ainit(i) = 0.0
     cinit(i) = 0.0
     a(i) = 0.0
     c(i) = 0.0
     rhoxtal(i) = 0.0
     rhodep(i) = 0.0
     phi(i) = 0.0
  END DO

  ! initial conditions
  press = Pinit
  zalt = 0.0
  vapmr = (es_new(TinitK)*(rhinit)*Rd)/(Rv*press)
  y(ivap) = vapmr
  y(iliq) = 0.0
  y(iice) = 0.0
  y(ipress) = Pinit
  y(itemp) = TinitK
  rhoa = press/(Rd*TinitK)
  firsthomog = .false.
  Sl = rhinit
  tempk = TinitK
  es = es_new(TinitK)
  
  ! Form solution drops from dry aerosol distribution.
  ! Solution drops assumed to be in equilibrium with the vapor field initially.
  CALL LiqNucleation(rlo, y, Nl, Nmrl, rhoa, ravgl, &
       rhodrop, ml, LWC, MRl, dMRl, Nlo, & 
       NUl, rhinit,nlinit,rdrp,rdrp0,rdryaer,sigmaccn)

  ! time integration loop
  tc = 0.0
  k = 0
  Do while (tc.lt.float(nt)*dt.and..not.endloop)

	w= wmax*sin(tc/tau_osc)
	if(iwvel.eq.2) w = wmax
	
	! get some needed values
	call CalcVars(y(itemp), y(ipress), y(ivap), &
	Sl, es, sui, rhoa,ei,suiratio,Dv,Kt,v_v,Ls,drhovidT)
	
	! nucleate ice
	if(nbicelim.lt.nbice.and.y(itemp).lt.To)then
	call homhetfreezing(ao, co, ainit, cinit, rio, y, Ni, Nmri, rhoa, &
		 rhoxtal, m, IWC, MRi,nbicelim,msalt, &
		 ivap,nl,dt,itemp,nfr,nlinit,nmrl, &
		 iwcnuc,nhetin,radiusin,ipress,rnhetin,nactmr,thet, &
		 sui,sl,naer,hetmodel,MRnucdep,MRnucfreez, &
		 areain,Nlo,nhetfrz,firsthomog,iice,iliq,rhodrop, &
		 ihabitchoice, habit)
	endif
     
	! Critical supersaturations for edge/edge growth assumption (M=10) from Julian's qualifying exam
	CALL findcritsats(y(itemp) - 273.15, scrit_a, scrit_c)
	
	! Calculate alpha
	do i = 1,nbicelim
		j = nbliq + i
		if(y(j) .ne. 0.0)then
		
		   ! Calculate alpha values for each ice crystal bin.
		   if(iicemodel.eq.0)then

		   	alpha_a(i) = 1.
		   	alpha_c(i) = 1.
		    rhodep(i) = rhoi

		   elseif(iicemodel.eq.1)then

		   	alpha_a(i) = alphaSphere
		   	alpha_c(i) = alphaSphere
		    rhodep(i) = rhoi

		   elseif(iicemodel.eq.2.or.iicemodel.eq.3.or.iicemodel.eq.4)then
		   
		    ! Adjustment of scrit for rosettes in platelike regimes to avoid impossible geometry

		   	if(habit(i).ge.4)then
		   	
		   		if(habit(i).le.6)then
                	pyrang = 45. * (pi/180.)
        		else 
                	pyrang = asin(1. - pi/(2.*REAL(habit(i))))
       		    endif
        		hpyr = ao(i)/tan(pi/2. - pyrang)
		   	
		   		scrit_a_temp = scrit_a
		   		scrit_c_temp = scrit_c
		   	
		   		adj_factor = 10.
		   		scritrat = scrit_c/scrit_a
		   		adj_fraction = (2.*co(i) - hpyr * 1.2) / hpyr
		   		scritrat = MIN(scritrat, 1 + adj_fraction ** 3. * adj_factor) 
		   				   	
		   		if(ao(i).lt.bsize)then
		   		    ! Forced constraint to avoid blowing up the model
		   			scritrat = MIN(scritrat, 1.)
		   		endif
		   		
		   		scrit_a = scrit_c / scritrat

		   	endif
		   	
		   	! Predict hollowing fraction from laboratory measurements at -40C
		   	! (effects of hollowing on growth not implemented yet)
		   	if(ihollow.eq.0)then
		   		xhol(i) = 0.
		   	elseif(ihollow.eq.1)then
		   		xhol(i) = xholinit
		   	elseif(ihollow.eq.2)then
		   		xhol(i) = hollowfrac(sui, scrit_c)
		   	endif
		   
		   	call calcfluxes(Fa, Fc, ao(i), co(i), .true., Dv, Kt, Ls, drhovidT, ei, &
		   	y(itemp), v_v, habit(i), xhol(i), Ab, Absol, Ap)
		   					
			! Calculate deposition density (from ratio of flux rates)
			rhodep(i) = rhoi * (Fa * Ap + Fc * Ab)/(Fa * Ap + Fc * Absol)
		   	
		   	! Revert scrit to previous value for next crystal
		   	if(habit(i).ge.4)then
		   		scrit_a = scrit_a_temp
		   		scrit_c = scrit_c_temp
		   	endif
	
		   	! Update alpha values for this ice bin 
		   	! (alpha_ain and alpha_cin updated through common block)
		   	alpha_a(i) = alpha_ain
		   	alpha_c(i) = alpha_cin

		   endif
		   
		   if(iicemodel.eq.0.or.iicemodel.eq.1)then
		    ! spheres; a = c.
			a(i) = y(j)
			c(i) = y(j)
		   elseif(iicemodel.eq.2.or.iicemodel.eq.4)then
			! a and c will be updated after the vode call instead
			a(i) = ao(i)
            c(i) = co(i)
		   elseif(iicemodel.eq.3)then
			! Update a and c by simple timestepping
			! a is an array to store the updated a and c dimensions in
			! until the timestep is completed
			a(i) = ao(i) + (Fa/rhoi) * dt
			c(i) = co(i) + (Fc/rhoi) * dt
			if(a(i).le.1.e-7)then
				a(i) = 1.e-7
			endif
			if(c(i).le.1.e-7)then
				c(i) = 1.e-7
			endif		   	            
			if(y(j).le.0.)then
				a(i) = 0.
				c(i) = 0.
			endif
		   endif
	
		   endif
	enddo

     istate = 1
     t = tc
     tout = t + dt 
     
     CALL vode(FunDiff, neqtot, y, t, tout, itol, rtol, atol, &
          itask, istate, iopt, rwork, lrw, iwork, liw, jac, &
          mf, rpar, ipar)
     
     ! find some characteristic information

     rhoa = y(ipress)/(Rd*y(itemp))
     ravg = 0.0
     rcavg = 0.0
     Nctot = 0.0
     Ntot = 0.0
     do i = 1,nbliq

        Nl(i) = Nmrl(i)*rhoa
        if (y(i) .ne. 0.0) then
           Ntot = Ntot + Nmrl(i)*rhoa
           ravg = ravg + y(i)*1.e6*Nmrl(i)*rhoa
           if (y(i) .gt. 1.5e-6)then
              Nctot = Nctot + Nmrl(i)*rhoa
              rcavg = rcavg + y(i)*1.e6*Nmrl(i)*rhoa
           endif
        endif
        
     enddo
     ravg = ravg/(max(1.0,Ntot))
     rcavg = rcavg/(max(1.0,Nctot))

     Ntoti = 0.0
     ravgi = 0.0
     aavg = 0.0
     cavg = 0.0
     do i = 1,nbicelim
        j = nbliq + i
        if(y(j) .ne. 0.0)then
        	
           if (iicemodel.eq.0.or.iicemodel.eq.1) then
           
           ! y = a = c
           rio(i) = y(j)
           ao(i) = y(j)
		   co(i) = y(j)
           
           elseif (iicemodel.eq.2) then
              ! Calculate volume change
              Vo = (4./3.)*pi*rio(i)**3.
              V =  (4./3.)*pi*y(j)**3.             
                          
              !print*, 'depd, effd', rhodep(i), rhoxtal(i)
              
              ! Update effective density
              rhoxtal(i) = rhoxtal(i) * (Vo/V) + rhodep(i) * (1. - Vo/V)
           
              ! Update phi
              phio = co(i)/ao(i)
              IGR = alpha_c(i)/alpha_a(i)
              !IGR = 2.1
              
              phi(i) = phio * (V/Vo) ** ((IGR/phio - 1.)/(IGR/phio + 2.))
              !phi(i) = phio * (V/Vo) ** ((IGR - 1.)/(IGR + 2.))
              !print*,a(i),c(i),phi(i),y(j)
              
              if(habit(i).ge.4)then

		   		if(habit(i).le.6)then
                	pyrang = 45. * (pi/180.)
        		else 
                	pyrang = asin(1. - pi/(2.*REAL(habit(i))))
       		    endif
       		    ! I did the math wrong!!
       		    
       		    ! xfac is a constant value corresponding to the volume loss from the pyramidal region
        		xfac = (2./3.)/tan(pi/2. - pyrang)
        		        		
        		if((2 * phi(i) - xfac).le.0)then

        			print*,"Error: rosette Aprism < 0! (stopping)"
        		
        			STOP
        		
        		endif
        		
        		a(i) = (V/(pi * REAL(habit(i)) * (2 * phi(i) - xfac)))**(1./3.)
        		c(i) = a(i) * phi(i)
        		
        		hpyr = a(i) / tan(pi/2. - pyrang)

        	  elseif(habit(i).eq.1)then
              	a(i) = (V/(2. * pi*phi(i)))**(1./3.)
              	c(i) = a(i) * phi(i)

              endif
              
              ! Update phi, initial radius of an equivalent volume sphere, and a and c arrays
              rio(i) = y(j)
              ao(i) = a(i)
              co(i) = c(i)

              if (y(j).le.0.) then
                 print*,"Warning: ice crystal volume is negative!"
                 y(j) = 0.
                 ao(i) = 0.
                 co(i) = 0.
              endif
		   	
           elseif(iicemodel.eq.3.or.iicemodel.eq.4)then

              ! Update beginning-of-timestep dimensions
              rio(i) = y(j)
              ao(i) = a(i)
              co(i) = c(i) 
		   	
           endif
           
           !phi(i) = phio*((V/Vo)**((IGR(i) - 1.)/(IGR(i) + 2.)))
           
           !ao(i) = ((3./(4.*pi))*(V/phi(i)))**(1./3.)   
           !co(i) = phi(i)*ao(i) 
        
           Ntoti = Ntoti + Nmri(i)*rhoa
           ravgi = ravgi + y(j)*1.e6*Nmri(i)*rhoa
           aavg =  aavg + ao(i)*1.e6*Nmri(i)*rhoa
           cavg =  cavg + co(i)*1.e6*Nmri(i)*rhoa
           xholavg =  xholavg + xhol(i)*Nmri(i)*rhoa

        endif
     enddo
     ravgi = ravgi/(max(1.e-10,Ntoti))
     aavg = aavg /(max(1.e-10,Ntoti))
     cavg = cavg /(max(1.e-10,Ntoti))
     xholavg = xholavg / (max(1.e-10,Ntoti))
     if(ihollow.eq.0.or.xholavg.ge.1.)then
     	xholavg = 0.
     endif

     zalt = zalt + w*dt
     qlsat = es_new(y(itemp))*Rd/(Rv*y(ipress))
     RH = y(ivap)/qlsat * 100.
     write(20,'(45E16.8)')tc,zalt,y(itemp)-To,y(ipress)/100. &
          ,y(ivap)*1000.,y(iliq)*1000.,y(iice)*1000.,RH, &
          Ntot/100.**3,Nctot/100.**3,ravg,rcavg,Ntoti/100.**3,ravgi, &
          aavg, cavg, xholavg
     
     k = k + int(dt)
     tc = tc + dt
  enddo
     
  stop
end program parcelmodel

SUBROUTINE homhetfreezing(ao, co, ainit, cinit, rio, y, Ni, Nmri, rhoa,  &
     rhoxtal, m, IWC, MRi,nbicelim,msalt,ivap,nl,dt,itemp, &
     nfr,nlinit,nmrl,iwcnuc,nhetin,radiusin,ipress,rnhetin,nactmr, &
     thet,sui,sl,naerin,hetmodel,MRnucdep,MRnucfreez, &
     areain,Nccn,nhetfrz,firsthomog,iice,iliq,rhodrop,ihabitchoice,habit) 
  
  IMPLICIT NONE

  INTEGER i,ji,nbice,ne, initdistflag, neqtot, nbliq, il, &
       nbicelim, ili,itemp,ivap,ipress,hetmodel,iice,iliq,ihabitchoice
  PARAMETER(nbice = 5000, ne = nbice + 1, nbliq = 100, &
       neqtot = nbice + nbliq + 5)
  REAL ao(nbice), co(nbice), rio(nbice), y(neqtot), phi, Ni(nbice), &
  	   ainit(nbice), cinit(nbice), &
       Nmri(nbice), rhoa, rhoxtal(nbice), m, pi, IWC, MRi,  &
       msalt(nbliq),aw(nbliq),rhosalt,vdrop,sl, &
       nl(nbliq),dt,vdropcm,nfr(nbice),nlinit(nbliq), &
       nmrl(nbliq),iwcnuc,lf,FLHM,cp,nhetin,radiusin,flhs,flhv,ls,lv, &
       evap,es_new,Nfrz,rfrz,Mfrz,ncnt,nfrzinst,nmrlsave,rnhetin, &
       nactmr,thet,sui,ad,bd,cd,dd,pow,naer,rhostp,naerin,nact,nin, &
       Msol,Mwat,vhf,MRnucdep,MRnucfreez,iwcnucfreez,iwcnucdep, &
       areain,chet,mhet,Nfrzhet,Jhet,xhet,volin,radin,fracin,Nccn, &
       nhetfrz(nbliq),sfctens,rhol,Acurve,rhodrop(nbliq)
  real aw_i,daw,x,j,ffr,sumni,sumnl,P,Rd,tempk,presspa,MRnuc
  real aw2,vhoff,molmass,rhosolute,rhoi
  real tc,tempc,mmolec,dgact,dgactmks,sigiw,lem,lhm0,rgg,tmelt, &
       t0,rgerm1,rhoihet,rgerm,xx,xm,tfmx1,tfmx2,tfmx3, &
       tfmx4,fmx,deltag,Jhetclss,kboltz,hplnk,cis,kb,rstar, &
       cangle,deltan,nltotinit,numinaer,numindrop,MRliq,MRice
  real a0m,a1m,a2m,a3m,a4m
  real randomnum
  integer habit(nbice)
  logical firsthomog
  parameter (a0m=79.7,a1m=-0.12,a2m=-8.0481e-2,a3m=-3.2376e-3)
  common/ccnchar/vhoff,molmass,rhosolute

  ! for classical heterogeneous theory
  kboltz = 1.381e-23
  kb = 1.381e-16
  hplnk = 6.626e-34
  Rstar = 8.314
  rhoihet = 0.95
  cangle = 58.6
  cis = 1.e15 !* 100.**2
  T0 = 273.15

  ! some initial constants and values
  cp =  1005.0 !J/KgK
  Rd = 287.0
  pi = 3.141529
  tempk = y(itemp)
  tc = tempk-273.15
  presspa = y(ipress)
  P = presspa

  rhosalt = rhosolute
  vhf = vhoff
  Msol = molmass
  Mwat = 0.018015             !kg
  sfcTens = 7.5e-2          !N/m
  rhol = 1000.0
  rhoi = 917.

  Acurve = (2.0*Mwat*sfcTens)/(Rstar*tempk*rhol) !curvature term

  ! Do Deposition Freezing First
  iwcnuc = 0.0
  
  if(nbicelim.eq.0.and.hetmodel.eq.1)then  ! nucleate only first time for het
     print*,'Instantaneous heterogenous nucleation selected'
     print*,'Heterogenous INP concentration ', nhetin, ' per L'
     print*,'Initial ice radius (assumed a = c) ', radiusin * 1.e6, ' microns' 
     nbicelim = nbicelim + 1
     ji = nbliq + nbicelim
     ni(nbicelim) = nhetin*1000.  ! concentration comes in per liter
     y(ji) = radiusin
     ao(nbicelim) = y(ji)
     co(nbicelim) =  y(ji)
     ainit(nbicelim) = y(ji)
     cinit(nbicelim) = y(ji)
     rio(nbicelim) = y(ji)
     rhoxtal(nbicelim) = 917.0
     nmri(nbicelim) = ni(nbicelim)/rhoa
     !nactmr = nactmr + nmr(nbicelim)*rhoa
     iwcnuc = iwcnuc + ni(nbicelim)*4./3.*pi*y(ji)**3.* &
          rhoxtal(nbicelim)/dt
     nactmr = nactmr + nhetin/rhoa
     
     ! Assign habit
     if(ihabitchoice.eq.1)then
     	habit(nbicelim) = 1
     elseif(ihabitchoice.eq.2)then
        ! rosette with a random number of arms (4-12)
     	CALL RANDOM_NUMBER(randomnum)
        habit(nbicelim) = 4 + INT(randomnum * 9.)
     elseif(ihabitchoice.eq.3)then
        CALL RANDOM_NUMBER(randomnum)
        ! rosette fraction specified by greater than term
        if(randomnum.gt.0.5)then
        	habit(nbicelim) = 1
        else
        ! rosette with a random number of arms (4-12)
     	CALL RANDOM_NUMBER(randomnum)
        habit(nbicelim) = 4 + INT(randomnum * 9.)
        endif
     elseif(ihabitchoice.ge.4)then
     	habit(nbicelim) = ihabitchoice
     endif

  endif

  rhostp = (1013.25*100.)/(289.*273.15)
  !naer = rnhetin*rhoa*1000. / 100.**3 !400.  ! assume 400 per cc
  naer = naerin             ! stp value of naer
  ad = 0.0000594
  bd = 3.33
  cd = 0.0264
  dd = 0.0033

  if((nbicelim.lt.nbice).and.(sui*100..gt.1.0).and. &
       (thet-tempk).ge.0.1)then
     
     thet = tempk
     pow = cd*(273.16-tempk) + dd
     nin = min((ad*(273.16-tempk)**bd * naer**pow)/rhostp*rhoa, &
          1000.0)           ! units of per liter
     nact = max(nin-nactmr*rhoa,0.0)
     if(hetmodel.eq.2.and.nact.gt.0.0)then
        nactmr = nactmr + nact/rhoa
        nbicelim = nbicelim + 1
        ji = nbliq + nbicelim
        ni(nbicelim) = nact*1000.  ! convert to #/m^3
        y(ji) = radiusin
        ao(nbicelim) = y(ji)
        co(nbicelim) =  y(ji)
    	ainit(nbicelim) = y(ji)
        cinit(nbicelim) = y(ji)
        rio(nbicelim) = y(ji)
        rhoxtal(nbicelim) = rhoi
        nmri(nbicelim) = ni(nbicelim)/rhoa
        iwcnuc = iwcnuc + ni(nbicelim)*4./3.*pi*y(ji)**3.* &
             rhoxtal(nbicelim)/dt
        
     endif

  endif

  iwcnucdep = iwcnuc   ! store for latent heating calculation
  iwcnuc = 0.0

  Nfrzhet = 0.0
  !homogeneous freezing
  rfrz = 0.0
  Nfrz = 0.0
  Mfrz = 0.0
  ncnt = 0.0
  do i = 1,nbliq

     vdrop = 4./3.*pi*y(i)**3
     vdropcm = vdrop*100.**3
     aw(i) = 1.0 - vhf*(msalt(i)/rhosalt)/(vdrop &
          + vhf*(msalt(i)/rhosalt))
     aw(i) = 1.0 - vhf*(msalt(i)/Molmass)/ & ! from P&K eqn 4-62
          ((vdrop-msalt(i)/rhosalt)*1000./Mwat + &
          vhf*msalt(i)/Molmass)

     aw_i=exp((210368.+131.438*tempk- 3.32373e6/tempk- & ! JYH Fixed Coefficent
          41729.1*log(tempk))/(8.31441*tempk))

     daw=aw(i)-aw_i
     if(daw.ge.0.26)then    ! JYH Set Limit
        x=-906.7+8502.*daw-26924.*daw**2+29180.*daw**3
        J=10**x
     else
        J = 0.0
     endif

     Ffr=1.-exp(-J*vdropcm*dt)
     nmrlsave = nmrl(i)
     nmrl(i) = max(0.0,nmrl(i) - ffr*nmrl(i))
     nfrzinst = max(0.0,nmrlsave - nmrl(i))
     nl(i) = nmrl(i)*rhoa
     
     if(.not.firsthomog)then
        if(nfrzinst.gt.0.0) firsthomog = .true.
     endif
     Nfrz = Nfrz + nfrzinst
     rfrz = rfrz + nfrzinst * y(i)*(rhodrop(i)/rhoi)**0.333333
     Mfrz = Mfrz + nfrzinst * &
          4./3.*pi*(y(i))**3. * rhodrop(i)
     if(ffr.gt.0.0) ncnt =ncnt + 1.

     ! Knopf heterogeneous immersion freezing  
     if(hetmodel.eq.3)then
        chet = -10.54758    !Kaolinite
        mhet = 54.58834
        numinaer =  100. * 100.**3  ! # m^-3 aerosol range of 10^-4 to 1 cm^-3
        numindrop = 1.
        fracin = numinaer/Nccn
        xhet = chet+mhet*daw
        Jhet = 10**xhet     ! units of cm^2/s 
        volin = 0.1*msalt(i)/rhosalt ! 10% of CCN volume is IN
        radin = 100.e-9
        volin = 4./3.*pi*radin**3
        numindrop = msalt(i)/rhosalt/volin       
        areain = numindrop*4.*pi*radin**2 * 100**2
        Ffr=1.-exp(-Jhet*areain*dt)   
        nmrlsave = nmrl(i)
        nmrl(i) = max(0.0,nmrl(i) - ffr*nmrl(i)*fracin) !nlinit(i)*fracin)
        nfrzinst = max(0.0,nmrlsave - nmrl(i))
        deltan = nlinit(i) - nmrlsave
        deltan = nhetfrz(i)
        if(deltan/nlinit(i).gt.fracin)then
           nfrzinst = 0.0
           nmrl(i) = nmrlsave
        else if(deltan/nlinit(i).lt.fracin.and. &
             (deltan+nfrzinst)/nlinit(i).gt.fracin)then
           nfrzinst = max(0.0,fracin*nlinit(i) - deltan)
           nmrl(i) = nmrlsave-nfrzinst
        endif
        nl(i) = nmrl(i)*rhoa
            
        nhetfrz(i) = nhetfrz(i) + nfrzinst
        Nfrzhet = Nfrzhet + nfrzinst ! ffr*nlinit(i)
        Nfrz = Nfrz + nfrzinst
        rfrz = rfrz + nfrzinst * y(i)*(rhodrop(i)/rhoi)**0.33
        Mfrz = Mfrz + nfrzinst * &
             4./3.*pi*(y(i))**3. * rhodrop(i)
        if(ffr.gt.0.0) ncnt =ncnt + 1.
     endif

     !Classical nuc theory
     if(hetmodel.eq.4)then
        numinaer =  naer * 100.**3  ! # m^-3 aerosol range of 10^-4 to 1 cm^-3
        fracin = numinaer/Nccn
        if(i.eq.100) print*,'numinaer is ',numinaer/1000.,fracin
        !stop
        volin = 0.1*msalt(i)/rhosalt  ! 10% of CCN volume is IN
        radin = (3./(4.*pi)*volin)**(1./3.) * 100.
        areain = 4.*pi*radin**2 
        
        tempc = max(tc,-45.0)
        if(tempc.gt.-30.0)then
           mmolec = 18.015/1000.
           dgact = 5.55*exp(-8.423e-3*tempc + 6.384e-4*tempc**2 + &
                7.891e-6*tempc**3)/(6.02e23*2.39e-11) ! kcal/mol
           dgactmks = dgact 
        else if(tempc.le.-30.0)then
           mmolec = 18.015/1000. !/6.022e23
           dgact = 0.694e-12*(1. + 0.027*(tempc+30.)* &
                exp(0.01*(tempc+30.d0))) ! erg per molec
           dgactmks = dgact !* 1e-7/mmolec
        endif
        sigiw = (28.+0.25*(tempc)) !*1.e-7*100**2
        lem = 79.7 + 0.708*tempc - 2.5e-3*tempc**2 ! cal/gram
        lhm0 = lem*4184.
        lem = lem * 4.184 *1.e7 ! mks
        rgg =  461.5*tempk/lhm0
        tmelt = t0*(aw(i)*exp(Acurve/y(il)))**rgg
        rgerm1 = 2.*sigiw/(lem*rhoihet*alog(tmelt/tempk))
        rgerm = max(rgerm1,1.e-40)
        xx = radin/rgerm
        xm = 0.5            
        phi = (1.-2.*xm*xx + xx**2)**(0.5)
        tfmx1 = ((1.-xm*xx)/phi)**3
        tfmx2 = 3.*((xx-xm)/phi)
        tfmx3 = ((xx-xm)/phi)**3
        tfmx4 = 3.*xm*xx**2*(((xx-xm)/phi) -1.)
        fmx = 0.5*(1.+tfmx1+xx**3*(2.-tfmx2+tfmx3)+tfmx4)
        deltaG = 4./3.*pi*sigiw*rgerm**2*fmx
        Jhetclss = 0.0
        if(rgerm1.gt.0.0)then
           Jhetclss = kboltz/hplnk*tempk*areain*cis* &
                exp( -dgactmks/(kb*tempk)-  &
                deltaG/(kb*tempk))
        endif
        if(Jhetclss.lt.1.e-90)then
           Jhetclss = 0.0
        endif
            
        Ffr = 1.-exp(-Jhetclss*dt)
        
        nmrlsave = nmrl(i)
        nmrl(i) = max(0.0,nmrl(i) - ffr*nmrl(i)*fracin)
        nfrzinst = max(0.0,nmrlsave - nmrl(i))
        deltan = nlinit(i) - nmrlsave
        deltan = nhetfrz(i)
        if(deltan/nlinit(i).gt.fracin)then
           nfrzinst = 0.0
           nmrl(i) = nmrlsave
        else if(deltan/nlinit(i).lt.fracin.and. &
             (deltan+nfrzinst)/nlinit(i).gt.fracin)then
           nfrzinst = max(0.0,fracin*nlinit(i) - deltan)
           nmrl(i) = nmrlsave-nfrzinst
        endif

        nl(i) = nmrl(i)*rhoa
        if(i.eq.1) print*,'Jhetclss = ',Jhetclss,tempc,Ffr,nfrzinst, &
             nmrl(i),nmrlsave,nlinit(i),ffr*nlinit(i)*fracin
        
        nhetfrz(i) = nhetfrz(i) + nfrzinst
        Nfrzhet = Nfrzhet + nfrzinst ! ffr*nlinit(i)
        Nfrz = Nfrz + nfrzinst
        rfrz = rfrz + nfrzinst * y(i)*(rhodrop(i)/rhoi)**0.33
        Mfrz = Mfrz + nfrzinst *  &
             4./3.*pi*(y(i))**3. * rhodrop(i)
        if(ffr.gt.0.0) ncnt = ncnt + 1.
     endif

  enddo
  
  if(hetmodel.eq.3.or.hetmodel.eq.4) nactmr = nactmr + Nfrzhet/1000.

  if(Nfrz.gt.0.0.and.nbicelim.lt.nbice)then
     nbicelim = nbicelim + 1
     if(nbicelim.ge.nbice)then
        print*,"ERROR: Out of ice bins!"
        STOP
  	 endif
     ji = nbicelim + nbliq
     nfr(nbicelim) = Nfrz
     ni(nbicelim) = nfr(nbicelim)*rhoa
     y(ji) = (Mfrz/Nfrz * 3./(4.*pi*rhoi))**0.33333333
     ao(nbicelim) = y(ji)
     co(nbicelim) = y(ji)
     ainit(nbicelim) = y(ji)
     cinit(nbicelim) = y(ji)
     rio(nbicelim) = y(ji)
     nmri(nbicelim) = ni(nbicelim)/rhoa
     rhoxtal(nbicelim) = 917.0
     iwcnuc = iwcnuc + ni(nbicelim)*4./3.*pi*y(ji)**3.* &
          rhoxtal(nbicelim)/dt
          
     ! Assign habit
     if(ihabitchoice.eq.1)then
     	habit(nbicelim) = 1
     elseif(ihabitchoice.eq.2)then
        ! rosette with a random number of arms (4-12)
     	CALL RANDOM_NUMBER(randomnum)
        habit(nbicelim) = 4 + INT(randomnum * 9.)
     elseif(ihabitchoice.eq.3)then
        CALL RANDOM_NUMBER(randomnum)
        ! rosette fraction specified by greater than term
        if(randomnum.gt.0.5)then
        	habit(nbicelim) = 1
        else
        ! rosette with a random number of arms (4-12)
     	CALL RANDOM_NUMBER(randomnum)
        habit(nbicelim) = 4 + INT(randomnum * 9.)
        endif
     elseif(ihabitchoice.ge.4)then
     	habit(nbicelim) = ihabitchoice
     endif
  endif

  iwcnucfreez = iwcnuc
      
  ! change temperature due to freezing of liquid water

  Ls = FLHS(tempk)
  Lv = FLHV(tempk)
  !lf = FLHM(tempk)/1.e4
  lf = ls-lv
  !evap = y(isat)*es_new(tempk)
  !rhoa = P/(Rd*tempk)
  MRnuc = iwcnuc/rhoa
  MRnucdep = iwcnucdep/rhoa
  MRnucfreez = iwcnucfreez/rhoa
  !print*,'temp before freezing ',tempk-273.15,iwcnuc*1000.
  tempk = tempk + lf*iwcnucfreez/cp*dt
  y(itemp) = tempk
  y(iice) = y(iice) + min(Mfrz,y(iliq)) + MRnucdep
  y(iliq) = y(iliq) - min(Mfrz,y(iliq))
  y(ivap) = y(ivap) - MRnucdep   ! remove deposition nucleation from vapor mixing ratio

  ! recompute liquid and ice mixing ratios
  MRliq = 0.0
  do i = 1,nbliq
     MRliq = MRliq + 4./3.*pi*y(i)**3 * rhodrop(i) * Nmrl(i)
  enddo
  y(iliq) = MRliq

  MRice = 0.0
  do i = 1,nbicelim
     ji = i  + nbliq
     MRice = MRice + 4./3.*pi*y(ji)**3 * rhoxtal(i) * Nmri(i)
  enddo
  y(iice) = MRice

  return
end SUBROUTINE homhetfreezing

SUBROUTINE FunDiff(neqtot, t, y, dy, rpar, ipar)
       
  IMPLICIT NONE
  
  integer neqtot,neq,itemp,ipress,ivap,iliq,iice,nbliq,nbice
  integer nbicelim,i,j
  PARAMETER(nbliq = 100,nbice=5000,neq=5)
  real y(neqtot),dy(neqtot),t,w
  real rhodrop(nbliq),nl(nbliq),rdrp(nbliq),rdrp0(nbliq)
  real Nmrl(nbliq),msalt(nbliq),vhoff,molmass,rhosolute
  real rdryaer(nbliq),Nmri(nbice),GTP,Sui,qisat,es,ei,es_new,esi_new
  real Rd,Rv
  real ao(nbice), co(nbice), rio(nbice), alpha_a(nbice), alpha_c(nbice), xhol(nbice)
  real a(nbice), c(nbice), rhoxtal(nbice), rhodep(nbice)
  integer ipar(*)
  integer habit(nbice)
  real rpar(*)

  COMMON/neqtotvars/itemp,ipress,ivap,iliq,iice,nbicelim
  COMMON/envvars/w,rhodrop, Nl,rdrp,rdrp0,rdryaer,Nmrl,Nmri
  COMMON/ccnvars/msalt
  common/ccnchar/vhoff,molmass,rhosolute
  COMMON/icevars/ao,co,a,c,rio,alpha_a,alpha_c,habit,xhol,rhoxtal,rhodep

  CALL liqGrowth(neqtot, y, dy, rhodrop, itemp, &
       ipress, ivap, Nl,msalt,rdrp,rdrp0, &
       rdryaer,GTP)

  call iceGrowth(neqtot, nbice, nbliq, y, dy, itemp, &
       ipress, ivap, Nmri, ao, co, a, c, rio, alpha_a, alpha_c, habit, xhol, rhoxtal, rhodep)
  
  CALL EnvEqns(neqtot,y,dy,w,itemp,ipress,ivap,iliq,iice, &
       nbliq,Nmrl,rhodrop,nbicelim,nbice,Nmri,rhoxtal,rhodep)
       
  return
end SUBROUTINE FunDiff

SUBROUTINE iceGrowth(neqtot, nbice, nbliq, y, dy, itemp, &
     ipress, ivap, Nmri, ao, co, a, c, rio, alpha_a, alpha_c, habit, xhol, rhoxtal, rhodep)

  implicit none
  integer neqtot,nbice,nbliq,i,j,itemp,ipress,ivap,iicemodel
  real y(neqtot),dy(neqtot),Nmri(nbice)
  real Rd,Rv,es,ei,qisat,Sui,P,Po,tempk,Ls,G,rhoi,Kt,Dv,Dvkin
  real es_new,esi_new,flhs
  
  !for flux calculations
  real v_v, pi, drhovidT, Ab, Absol, Ap
  real avp, athp, atot, cvp, cthp, ctot, hc, ha, cap
  real phio, phi, Rk, Rt, Fmax, drho, Fa, Fc, Favg, suiratio, rhoa, sl
  real Ma, Mc, scrit_a, scrit_c, alpha_ain, alpha_cin
  real ao(nbice), co(nbice), rio(nbice), alpha_a(nbice), alpha_c(nbice), xhol(nbice)
  real a(nbice), c(nbice), rhodep(nbice), rhoxtal(nbice)
  real Vo, V, IGR
  real rprev
  real pyrang, xfac
  integer habit(nbice)
  
  COMMON/inputvars/iicemodel
  COMMON/edgezzbren/sui,Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin
  
  pi = 3.14159265
  rhoi = 917.0
  
  ! Calculate values of many derived atmospheric variables
  CALL CalcVars(y(itemp), y(ipress), y(ivap), &
  Sl, es, sui, rhoa,ei,suiratio,Dv,Kt,v_v,Ls,drhovidT)
  
  do i = 1,nbice
  
     j = i + nbliq
     dy(j) = 0.0

     if (y(j).gt.0.0) then
     
                 
        ! When not recalculating alpha, set input variables (through common block)
        ! to previously calculated alpha value (outside of vode)
        alpha_ain = alpha_a(i)
        alpha_cin = alpha_c(i)
        
        if (iicemodel.eq.0) then
           
           G = ((((Ls/(Rv*y(itemp)))-1.)*(Ls/(Kt*y(itemp)))) + ((Rv*y(itemp)) &
             	/(Dv*ei)))**(-1.) ! maximum diffusivity
           
           dy(j) = (G*Sui)/(y(j)*rhoi)  
           
     	elseif (iicemodel.eq.1) then
         
         ! Kinetically modified diffusivity
         ! Note that we set alpha_a = alpha_c = const. in this model
           Dvkin = Dv / ((4 * Dv)/(alpha_a(i) * v_v * y(j)) + 1)
			
         ! The equation with a and c dimensions (currently not used here)
         !Dvkin = (2./3.) * Dv / ((4 * Dv * cap)/(alpha_a(i) * Ls * ao(i) * co(i)) + 1) + &
         !(1./3.) * Dv / ((4 * Dv * cap)/(alpha_c(i) * Ls * ao(i) * ao(i)) + 1)
         
            G = ((((Ls/(Rv*y(itemp)))-1.)*(Ls/(Kt*y(itemp)))) + ((Rv*y(itemp)) &
              /(Dvkin*ei)))**(-1.) ! effective diffusivity
         
            dy(j) = (G*Sui)/(y(j)*rhoi)
         
         
           elseif (iicemodel.eq.2.or.iicemodel.eq.3) then
            
            ! Calculate fluxes to prism and basal faces (appropriate variables updated)
            ! through common block edgezzbren
            call calcfluxes(Fa, Fc, ao(i), co(i), .false., Dv, Kt, Ls, drhovidT, ei, & 
            y(itemp), v_v, habit(i), xhol(i), Ab, Absol, Ap)
	        
	        if(habit(i).ge.4)then
				Ab = Ab * REAL(habit(i))
				Ap = Ap * REAL(habit(i))
			endif
        	
        	if (iicemodel.eq.3) then
        	    dy(j) = (Ab * Fc + Ap * Fa) / (4 * Pi * rhoi * y(j)**2.)
        	else
        		dy(j) = (Ab * Fc + Ap * Fa) / (4 * Pi * rhodep(i) * y(j)**2.)
        	endif
            
        elseif (iicemodel.eq.4) then
        
        	! Update dimensions according to output from the last VODE call
        	! (which apparently seems to work even if y decreases between iterations)
          	Vo = (4./3.)*pi*rio(i)**3.
           	V =  (4./3.)*pi*y(j)**3.
           
           	! Update phi
           	phio = co(i)/ao(i)
           	IGR = alpha_c(i)/alpha_a(i)

            phi = phio * (V/Vo) ** ((IGR/phio - 1.)/(IGR/phio + 2.))
                
            if(habit(i).ge.4)then

		   		if(habit(i).le.6)then
                	pyrang = 45. * (pi/180.)
        		else 
                	pyrang = asin(1. - pi/(2.*REAL(habit(i))))
       		    endif
       		    ! xfac is a constant value corresponding to the volume loss from the pyramidal region
        		xfac = (2./3.)/tan(pi/2. - pyrang)
        		        		
        		if((2 * phi - xfac).le.0)then

        			print*,"Error: rosette Aprism < 0! (stopping)"
        		
        			STOP
        		
        		endif
        		
        		a(i) = (V/(pi * REAL(habit(i)) * (2 * phi - xfac)))**(1./3.)

        	  elseif(habit(i).eq.1)then
              	a(i) = (V/(2. * pi*phi))**(1./3.)


            endif

           	c(i) = a(i) * phi
           	
           	! Calculate basal and prism fluxes ()
           	call calcfluxes(Fa, Fc, a(i), c(i), .false., Dv, Kt, Ls, drhovidT, ei, &
           	y(itemp), v_v, habit(i), Ab, Absol, Ap)
			
			! Account for surface area of all of the arms
			if(habit(i).ge.4)then
				Ab = Ab * REAL(habit(i))
				Ap = Ap * REAL(habit(i))
			endif
        	
        	dy(j) = (Ab * Fc + Ap * Fa) / (4 * Pi * rhodep(i) * y(j)**2.)
        	    
        	! Update value for previous radius        
        	!rio(i) = y(j)
        
      endif
   endif
  enddo

  return
end SUBROUTINE iceGrowth

SUBROUTINE liqGrowth(neqtot, y, dy, &
     rhodrop, itemp, ipress, ivap, Nl,msalt, &
     rdrp,rdrp0,rdryaer,G)

  IMPLICIT NONE
  
  INTEGER i, neqtot, nbice, itemp, ipress, isat, ivap, il, nbliq, &
       nbicelim
  PARAMETER(nbliq = 100)
  REAL Ls, Rv, ei, es, Dv, Kt, Lv, Sl, G, Nl(nbliq), &
       To, dt, eo, pi, gam, y(neqtot),t,rhodrop(nbliq), &
       es_new, esi_new, m, celsius, drho,IGR_DEN, dy(neqtot), &
       A, B, Mw, vhf, Ms, R, sfcTens, rhol, msalt(nbliq), &
       rdrp(nbliq),rdrp0(nbliq),fsoln,Slavg,dtime,aw,vdrop,nliq, &
       nsol,rhosalt,vsol,vliq,rhosoldrop,tempk,rdryaer(nbliq),epsr, &
       vhoff,molmass,rhosolute,qlsat,Rd,suiratio,v_v,rhoa,drhovidT,sui
  real FLHV
  common/ccnchar/vhoff,molmass,rhosolute

  dtime = 1.0
  To = 273.15               ! [=]K
  Rd = 287.0
  Rv = 461.5                ! Individual gas constant for water vapor[=]J/kgK
  Lv = 2.5E6                ! Latent heat of vaporization[=]J/kg
  pi = 3.14159
  Mw = 0.018015             !kg
  R = 8.314                 ! J/molK
  sfcTens = 7.5e-2          !N/m
  rhol = 1000.0
  tempk = y(itemp) 
  Lv = FLHV(tempk)

  Ms = molmass
  rhosalt = rhosolute
  vhf = vhoff
  
  CALL CalcVars(y(itemp), y(ipress), y(ivap), &
  Sl, es, sui, rhoa,ei,suiratio,Dv,Kt,v_v,Ls,drhovidT)
        
  G = ((((Lv/(Rv*tempk))-1.)*(Lv/(Kt*tempk))) + ((Rv*tempk) &
             /(Dv*es)))**(-1.)
      
  DO i = 1, nbliq

     il = i 

     epsr = rdryaer(i)*0.1
     If (y(i) .lt. rdryaer(i)+epsr .or. Nl(i) .eq. 0.0) THEN
        
        dy(i) = 0.0
        if (Nl(i) .eq. 0.0) y(i) = 0.0

     ELSE

        A = (2*Mw*sfcTens)/(R*tempk*rhol)

        B = (3*vhf*msalt(i)*Mw)/(4*pi*Ms*rhol)
                     
        call getdropchar(y(i),msalt(i),vsol,vliq,aw,rhosoldrop)

        rhodrop(i) = rhosoldrop

        fsoln = 1.-aw*exp(A/y(i))

        dy(i) = (G*(Sl - 1. + fsoln))/(y(i)*rhodrop(i))

        if(vliq.lt.0.0)then
           print*,'inside liqgrowth vsol > vdrop'
           print*,il,i,msalt(i),vsol,vliq,aw,y(i),rdrp0(i)
           dy(i) = 0.0               
           stop
        endif

     END IF

  END DO

END SUBROUTINE liqGrowth

subroutine EnvEqns(neqtot,y,dy,w,itemp,ipress,ivap,iliq,iice, &
     nbliq,Nmrl,rhodrop,nbicelim,nbice,Nmri,rhoxtal,rhodep)

  implicit none

  integer neqtot,itemp,ipress,ivap,iliq,iice,nbliq,i,j,nbicelim,nbice,iicemodel
  real w
  real y(neqtot),dy(neqtot)
  real g,Rd,Mw,Ma,eps,Rm,Cp,Lv,Ls,Rv,qlsat,rhoi
  real Nmrl(nbliq),rhodrop(nbliq),pi,Nmri(nbice),rhoxtal(nbice),rhodep(nbice)
  real FLHV,FLHS,es_new

  g = 9.81
  pi = 3.1415926
  Rd = 287.0
  Rv = 461.5
  rhoi = 917.0 
  Mw = 0.01801528	! Molar mass of water [kg/mole]
  Ma = 0.0289647	! Molar mass of dry air [kg/mole]
  eps = Mw/Ma
  Rm = (1.0 + (y(ivap)/eps))/(1.0 + y(ivap)) * Rd 		! Gas constant for moist air.
  Cp = 1004.5
  Lv = FLHV(y(itemp))
  Ls = FLHS(y(itemp))

  qlsat = es_new(y(itemp))*Rd/(Rv*y(ipress))
  
  dy(iliq) = 0.0

  do i = 1,nbliq
     dy(iliq) = dy(iliq) + 4.*pi*y(i)**2 * rhodrop(i) * Nmrl(i)*dy(i)
  enddo

  dy(iice) = 0.0
  do i = 1,nbicelim
     j = i + nbliq
     if (iicemodel.eq.3)then
     	dy(iice) = dy(iice) + 4.*pi*y(j)**2 * rhoi * Nmri(i)*dy(j)
     else
     	dy(iice) = dy(iice) + 4.*pi*y(j)**2 * rhodep(i) * Nmri(i)*dy(j)
     endif

  enddo
  dy(ipress) = -((g * w)/(Rm * y(itemp))) * y(ipress)
  dy(itemp) = -(g * w)/Cp + (Lv/((1.0 + y(ivap)) * Cp) * dy(iliq)) + (Ls/((1.0 + dy(ivap)) * Cp) * dy(iice))
  dy(ivap) = -(dy(iliq) + dy(iice))

  return
end subroutine enveqns

!*********************************************************************
!     This subroutine will initiate liquid drops based on the initial 
!     solute mass and liquid drop distribution.
!**********************************************************************

SUBROUTINE LiqNucleation(rlo, y, Nl, Nmrl, rhoa, ravgl, &
     rhodrop, ml, LWC, MRl, dMRl, Nlo, NUl,rhinit, &
     nlinit,rdrp,rdrp0,rdryaer,sigmaccn)
  
  IMPLICIT NONE

  INTEGER i, j, nbliq, ne, neqtot, nbice, jl,itmax
  PARAMETER(nbice = 5000, nbliq = 100, ne = nbliq + 1, &
       neqtot = nbliq + nbice + 5)
  REAL rlo(nbliq), y(neqtot), Nl(nbliq),Nlo, NUl, rnl, ravgl, &
       Nmrl(nbliq), rhoa, rhodrop(nbliq), ml, pi, LWC, MRl, dMRl, & 
       re(ne), Sl, tempk, msalt(nbliq), es, sigma, rlow, &
       rhigh, errabs, rel, rdry, msalt2, sig,nlinit(nbliq), &
       rdrp(nbliq),rdrp0(nbliq),vsol,vliq,aw,rhosoldrop, &
       rdryaer(nbliq),vhoff,molmass,rhosolute,sigmaccn,rhinit
  COMMON/aerosolzzbren/Sl, tempk, es, rdry, msalt2
  COMMON/ccnvars/msalt
  common/ccnchar/vhoff,molmass,rhosolute

  EXTERNAL Kohler

  Sl = rhinit
  pi = 3.141529
  sig = sigmaccn  ! geometric standard deviation 1.5 to 2.3
  sigma = log(sig)
  NUl = exp((sigma**2)/2.)
  rnl = ravgl   !ravgl/NUl
  open(51, file = "outdir/initialdrydist.dat")
  OPEN(52, file = "outdir/initialwetdist.dat")

  CALL RadEdgeLiq(Nlo, NUl, rnl, nbliq, ne, re, Nl, rlo, sig, ravgl) !Get dry aerosol distribution

  DO i = 1, ne-1
     WRITE(51,'(45E16.8)') re(i),Nl(i), rlo(i)
  END DO

  CLOSE(51)

  LWC = 0.0
  DO j = 1, nbliq
     rdry = rlo(j)
     rhodrop(j) = 1000.0
     rdryaer(j) = rdry
     itmax = 100
     rlow = rdry ! haze drops 1% larger than CCN due to swelling
     rhigh = 5.E-5
     errabs = 0.0
     rel = 1.E-11
     CALL ZZBREN (Kohler,errabs,rel,rlow,rhigh,itmax)
     y(j) = rhigh
     rdrp(j) = y(j)
     rdrp0(j) = rdrp(j)
     msalt(j) = msalt2
     call getdropchar(y(j),msalt(j),vsol,vliq,aw,rhosoldrop)
     if(vliq.lt.0.0)then
        print*,'problem in liquid nucleation...'
        print*,'vsol > vliq ',vsol,vliq
        stop
     endif
      
     Nmrl(j) = Nl(j)/rhoa
     nlinit(j) = nmrl(j)         
     rhodrop(j) = rhosoldrop
     ml = (4./3.)*pi*(y(j)**3)*rhodrop(j)
     LWC = LWC + Nl(j)*ml
  END DO
  DO j = 1,nbliq
     WRITE(52,'(45E16.8)') Nl(j), y(j), msalt(j)
  END DO

  MRl = LWC/rhoa
  
END SUBROUTINE LiqNucleation

!*********************************************************************
! Subroutine to calculate some derived variables from temperature, pressure and vapor pressure
! in standard SI units
! Outputs: liquid saturation ratio (Sl), liquid saturation vapor pressure (es), 
! ice supersaturation (sui), density of air (rhoa), ice saturation vapor pressure (ei),
! diffusivity of water vapor (Dv), thermal conductivity of air (Kt), mean thermal velocity of
! a vapor molecule (v_v), latent heat of sublimation (Ls), and the first derivative
! of saturation ice vapor density with respect to temperature (drhovidT)
!*********************************************************************
SUBROUTINE CalcVars(tempk, press, vpres, &
					Sl, es, sui, rhoa,ei,suiratio,Dv,Kt,v_v,Ls,drhovidT)

	 implicit none
	 
	 real tempk, press, vpres ! Input
	 real rhoa,ei,suiratio,Dv,Kt,v_v,Ls,drhovidT ! Output vars through common block derivedvars
	 real Sl, es, sui ! Output vars that are already in a common block elsewhere
	 real qeqi, qeql ! Intermediate vars
	 real esi_new, es_new, FLHS ! Functions
	 real Pi, Rv, Rd ! Constants
  	  
	 Pi = 3.14159265
	 Rd = 287.     ! dry air gas constant
  	 Rv = 461.5    ! water vapor gas constant

     rhoa = press/(Rd*tempk)
     ei = esi_new(tempk)
     es = es_new(tempk)
     qeqi = (ei*Rd)/(Rv*press)
     qeql = (es*Rd)/(Rv*press)
     sl = vpres/qeql
     sui = vpres/qeqi - 1.
     suiratio = sui/(es/ei-1.)
     	 
	 Dv = 0.211*(tempk/273.15)**1.94 * (1013.25/(press/100.)) * 1./100.**2. ! vapor diffusivity 
			
	 Kt = (5.69 + 0.017*(tempk-273.15)) * 418.684/1.e5 ! thermal diffusivity
				 
	 v_v = sqrt(8.*Rv*tempk/pi) ! mean speed of a vapor molecule
		 
	 Ls = FLHS(tempk)
		 
     drhovidT = (ei * ((Ls/(Rv*tempk)) - 1.))/(Rv * (tempk ** 2.))  ! change in eq vap density with temp  ! change in eq vap density with temp
     
END SUBROUTINE CalcVars

!*********************************************************************
! Get the fluxes to the basal and prism faces of an ice crystal; also gives the option to
! calculate the alpha values; otherwise make sure to set input values
! alpha_ain and alpha_cin to their last calculated value before calling this subroutine.
! These fluxes follow faceted growth theory with ledges nucleated at the crystal edge.
!*********************************************************************
SUBROUTINE calcfluxes(Fa, Fc, aoin, coin, icalcalpha, Dv, Kt, Ls, drhovidT, esi, tempk, &
						v_v, ihabit, xholin, Ab, Absol, Ap)

 implicit none
 
 logical icalcalpha
 integer itmax, ihabit
 real sui,Ma,Mc,scrit_a,scrit_c,ha,hc,hcinner,atot,ctot,alpha_ain,alpha_cin
 real errabs, rel1, para, parb
 real aoin, coin, xholin, ahol, phio, pi, Rk, Rt
 real Dv, Kt, drhovidT, esi, Rv, v_v, tempk, Ls
 real capacitance, cap, caphol
 real Fa, Fc, Ab, Absol, Ap, avp, cvp, athp, cthp, fmax
 
 ! variables specific to a rosette
 real pyrang, hpyr  
 
 external funslocal
 
 COMMON/edgezzbren/sui,Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin
 
 	pi = 3.14159265
 	Rv = 461.5

	! aspect ratio
	phio = coin/aoin
	
	! dimension of hollow center of basal face (will be zero for unhollowed crystals)
	ahol = xholin * aoin
	
	! facet areas and capacitance
	if(ihabit.eq.1)then
		Ab = 2. * Pi * (aoin**2. - ahol**2.)
		Absol = 2. * Pi * (aoin**2.)
		Ap = 4. * Pi * aoin * coin
		! capacitance of a cylinder
		cap = capacitance(aoin, coin, 1, 2)
		
		! capacitance of a cylinder with radius equal to that of the hollowed portion
		! of the basal face
		if(xholin.gt.0.01)then
			caphol = capacitance(aoin * xholin, coin, 1, 2) 
			hcinner = Ab / (4 * pi * caphol * coin)
		else
			! The capacitance function doesn't like it when a = 0.
			caphol = 0.
			hcinner = 0.
		endif
		
		! h-functions
		hc = Ab / (4. * pi * cap * coin) - hcinner
		hc = Ab / (4. * pi * cap * coin)

		ha = Ap / (4. * pi * cap * aoin)

	elseif(ihabit.ge.4)then
	    ! Facet areas (for one arm) 
        Ab = 2. * Pi * (aoin**2. - ahol**2.)
        Absol = 2. * Pi * (aoin**2.)
		Ap = 2. * Pi * aoin * (2.*coin - hpyr)
	
		! capacitance of a rosette
		cap = capacitance(aoin, coin, 3, ihabit)
		
		! capacitance of a rosette with radius equal to that of the hollowed portion
		! of the basal face
		if(xholin.gt.0.01)then
			caphol = capacitance(aoin * xholin, coin, 3, ihabit) 
			hcinner = (Ab * REAL(ihabit)) / (4 * pi * caphol * coin)
		else
			! The capacitance function doesn't like it when a = 0.
			caphol = 0.
			hcinner = 0.
		endif
		
		! Slope of pyramidal region of a rosette (use n-arms solution for 7+ arms)
		if(ihabit.le.6)then
                pyrang = 45. * (pi/180.)
        else 
                pyrang = asin(1. - pi/(2.*REAL(ihabit)))
        endif
        
        hpyr = aoin/tan(pi/2. - pyrang)
		
		! h-functions
		hc = (Absol * REAL(ihabit)) / (4. * pi * cap * coin) - hcinner
		ha = (Ap * REAL(ihabit)) / (4. * pi * cap * aoin)
	else
		print*,"In calcfluxes habit", ihabit, " is not a valid habit choice."
		STOP
	endif
	
	! thermodynamic quantities
	avp = (aoin * v_v) / (4. * Dv)
	athp = (aoin * v_v * drhovidT * Ls) / (4. * Kt)
	
	atot = avp + athp
	
	cvp = avp * phio
	cthp = athp * phio
	ctot = cvp + cthp
	
	! Inputs to Brent's routine and function call
	errabs = 0.0
	rel1 = 1.e-7
	itmax = 50
	para = sui*1.e-4
	parb = sui
	
	! Option to not call Brent's routine if run within VODE
	if(icalcalpha)then
		call zzbren(funslocal,errabs,rel1,para,parb,itmax)
	endif
	
	! Calculate resistances
	Rk = alpha_cin * cvp * hc + alpha_ain * avp * ha
	Rt = alpha_cin * cthp * hc + alpha_ain * athp * ha
	
	! Calculate fluxes
	Fmax = (0.25 * v_v * (sui * esi)/(Rv * tempk)) / (1. + Rk + Rt)
	Fa = Fmax * alpha_ain
	Fc = Fmax * alpha_cin
	! Only one basal face for rosette
	if(ihabit.ge.4)then
		Fc = Fc / 2.
	endif
	if(sui.le.0.)then
	   if(phio.gt.1)then
		  Fa = Fmax/phio
	   else
		  Fc = Fmax*phio
	   endif
	endif
 
 END SUBROUTINE calcfluxes

! Subroutine to get the critical supersaturations as a function of temperature
subroutine findcritsats(tempc, scrit_a, scrit_c)

	  implicit none
	  
	  real tempc, sc, sa, scrit_c, scrit_a, scrit_ratio

      if(tempc.lt.-30.0)then
         ! original magee (high scrit) 
         sc = 1.8115 + 0.15585*tempc + 0.011569*tempc**2
         sa = 1.8115 + 0.15585*tempc + 0.011569*tempc**2
         
         ! with bailey and hallett
         sc = 3.7955 + 0.10614*tempc + 0.0075309*tempc**2
         sa = sc
         ! modified magee
         !sc = 3.3649 + 0.13266*tempc + 0.008959*tempc**2
         !sa = sc
      elseif((tempc.ge.-30.0) .and. (tempc.lt.-22.0))then
         sc = 753.63 + 105.97 * tempc + 5.5532 * tempc**2 + 0.12809 * tempc**3 + 0.001103 * tempc**4
      elseif((tempc.ge.-22.0) .and. (tempc.lt.-1.0))then
         sc = 1.1217 + 0.038098 * tempc - 0.083749 * tempc**2 - 0.015734 * tempc**3 - 0.0010108 &
         * tempc**4 - 2.9148e-05 * tempc**5 - 3.1823e-07 * tempc**6
      endif
    
      if((tempc.ge.-30.0) .and. (tempc.lt.-22.0))then
         sa =  -0.71057 - 0.14775*tempc + 0.0042304*tempc**2
         if(((sc/sa).gt.0.92) .and. (tempc.lt.-25.0))then
             sc = 0.92 * sa
         endif
            
            
      elseif((tempc.ge.-22.0) .and. (tempc.lt.-15.0))then
         sa = -5.2367 - 1.3184*tempc - 0.11066*tempc**2 - 0.0032303*tempc**3
      elseif((tempc.ge.-15.0) .and. (tempc.lt.-10.0))then
         !sa = 0.755 + 0.03325*tempc + 0.00125*tempc**2   ! fit to Woods
         sa = 0.34572 - 0.0093029*tempc + 0.00030832*tempc**2 ! fit to Nelson & Night
      elseif((tempc.ge.-10.0) .and. (tempc.lt.-1.0))then
         !          sa = 0.37445 + 0.029636*tempc + 0.011575*tempc**2 &      ! Wood
         !               + 0.00069444*tempc**3
         sa = 0.34572 - 0.0093029*tempc + 0.00030832*tempc**2 ! Nelson & Knight
      
      elseif((tempc.lt.-30.0).and.(tempc.gt.-100.0))then
    
          scrit_ratio = 0.00334*tempc + 0.9825

          sa = -0.71057 - 0.14775*tempc + 0.0042304*tempc**2
          sc = sa * scrit_ratio
      else
      	 print*, "Temperature problem in findcritsats, stopping!", tempc
      	 sa = -32767
      	 sc = -32767
      	 STOP
      endif
      
      scrit_a = sa/100.
      scrit_c = sc/100.

  return
end subroutine findcritsats

!*********************************************************************
!     This subroutine will calculate and return the mass of the salt,
!     msalt, with a given relative humidity and radius of the liquid
!     drop.
!**********************************************************************

REAL FUNCTION Kohler(rdrop)
      
  IMPLICIT NONE

  REAL es, A, ravgl, B, Mw, sfcTens, R, tempk, rhol,vhf, &
       Ms,es_new, esol, pi, Sl, rdry, rhosalt, rdrop, Ssol, diff, &
       msalt2,vhoff,molmass,rhosolute,vsol,vliq,aw,rhosoldrop

  COMMON/aerosolzzbren/Sl, tempk, es, rdry, msalt2
  common/ccnchar/vhoff,molmass,rhosolute

  Mw = 0.018015             !kg
  R = 8.314                 ! J/molK
  sfcTens = 7.5e-2          !N/m
  rhol = 1000.0
  pi = 3.14159

  Ms = molmass
  vhf = vhoff
  rhosalt = rhosolute
      
  msalt2 = (4./3.)*pi*(rdry**3)*rhosalt 

  A = (2.0*Mw*sfcTens)/(R*tempk*rhol) !curvature term
  
  B = (3.0*vhf*msalt2*Mw)/(4*pi*Ms*rhol) !solution term

  call getdropchar(rdrop,msalt2,vsol,vliq,aw,rhosoldrop)         

  IF(rdrop .ne. 0.0) THEN 
     esol = Sl*es
     Ssol = exp(A/rdrop - B/rdrop**3)
     Ssol = aw*exp(A/rdrop)
  END IF

  diff = Ssol - Sl
  Kohler = diff
  
  RETURN
END FUNCTION Kohler

subroutine getdropchar(r,msol,vsol,vliq,aw,rhosoldrop)

  implicit none
  real r,msol,rhosol,vsol,vliq,aw,Ms,vhf
  real Mw,pi,rhol,vdrop,nsol,nliq,rhosoldrop
  real vhoff,molmass,rhosolute
  common/ccnchar/vhoff,molmass,rhosolute

  pi = 3.14159
  Mw = 0.018015
  rhol = 1000.
  
  vhf = vhoff
  Ms = molmass
  rhosol = rhosolute
  
  vdrop = 4./3.*pi*r**3
  vsol = msol/rhosol
  vliq = vdrop - vsol
  nsol =  msol/Ms/vdrop
  nliq = rhol*vliq/Mw/vdrop
  aw = nliq/(nliq + vhf*nsol)
  rhosoldrop = (rhosol*vsol + rhol*vliq)/vdrop
  
  return
end subroutine getdropchar


!**********************************************************************
!     This subroutine is used to locate multiple bin edges within a 
!     distribution. From this, the total # of drops, and their
!     sizes, for each bin are found and passed back into LiqNucleation.
!**********************************************************************
      
SUBROUTINE RadEdgeLiq(Nio, NU, rn, nb, ne, re, Ni, rlo, sig,rg)

  IMPLICIT NONE
      
  INTEGER i, nb, ne 
  REAL Nio, NU, rn, re(ne), rlow, rhigh,deltar,Ni(nb),rlo(nb)
  real sig,rg
  REAL factnorm, ntotal,val, pi
  
  pi = 3.14159
  rlow =  2.0e-8 
  rhigh = 1.0e-6                                                      
  re(1) = rlow
  re(ne) = rhigh
  deltar = (rhigh - rlow)/float(nb)
  
  DO i = 1, nb  ! Loop used to establish bin edges (201) for 200 bins
         
        ! re(i) = re(i-1) + deltar
     re(i) = val(rlow,rhigh,nb+1,i,'log')
         
  END DO
      
  DO i = 1, nb
         
         !Nigam = Nio*(gammp(NU,re(i+1)/rn) - gammp(NU,re(i)/rn)) 
!    # of particles in each bin
     rlo(i) = (re(i+1) + re(i))/2. ! size of particles in each bin
     Ni(i) = Nio/( (2.*pi)**0.5 * alog(sig)*rlo(i) ) * &
          exp( -(alog(rlo(i)/rg))**2/(2.*(alog(sig))**2) ) * &
          (re(i+1) - re(i))
     ntotal = ntotal + Ni(i)

  END DO
      
        !Normalize concentration so constant for all average radii
  factnorm = Nio/ntotal

  DO i = 1,nb
     Ni(i)=Ni(i)*factnorm
  END DO

  ntotal = 0.0
  DO i = 1,nb
     ntotal = ntotal + Ni(i)
  END DO
      
END SUBROUTINE RadEdgeLiq

function val(x0,x1,nx,ix,csca)
  character*(*) csca
  if(nx.le.1) print*,' error in val '
  fact=(real(ix)-1.)/(real(nx)-1.)
  if(csca(1:3).eq.'lin') then
     val = x0+(x1-x0)*fact
  elseif(csca(1:3).eq.'log') then
     val= x0*(x1/x0)**fact
  endif
  return
end function val


! using function from Murphy and Koop (2005)
REAL FUNCTION esi_new(tmp)
  REAL ci(0:3)
  data ci/9.550426, 5723.265, 3.53068, 0.00728332/

  alogpice = ci(0) - ci(1)/tmp + ci(2)*alog(tmp) - &
       ci(3)*tmp
     
  esi_new = exp(alogpice)

  RETURN
END FUNCTION esi_new

!using function from Murphy and Koop (2005)
real function es_new(tmp)

  DATA C0,C1,c2,c3,c4,c5,c6,c7,c8/54.842763, 6763.22, &
       4.210, 0.000367, 0.0415, 53.878, 1331.22, &
       9.44523, 0.014025/
  DATA CP,RCP/1004.0, 0.285856574/

  alpliq = c0 - c1/tmp - c2*alog(tmp) + c3*tmp + &
       tanh(c4*(tmp-218.8))*( c5 - c6/tmp - c7*alog(tmp) + &
       c8*tmp )
  es_new = exp(alpliq)

  RETURN
END function es_new



!**********************************************************************
!     Function for finding the surface supersaturations and
!	  deposition coefficients for a growing ice crystal with
!	  a root finding method.
!**********************************************************************
function funslocal(slocalguess)

	implicit none

	real sui,Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin, resistance
	real funslocal, slocalguess

    COMMON/edgezzbren/sui,Ma,Mc,scrit_a,scrit_c,ha,hc,atot,ctot,alpha_ain,alpha_cin

	alpha_ain = (slocalguess/scrit_a)**Ma * tanh((scrit_a/slocalguess)**Ma)
	alpha_cin = (slocalguess/scrit_c)**Mc * tanh((scrit_c/slocalguess)**Mc)

	resistance = alpha_ain * ctot + alpha_cin * atot
	
	funslocal = sui * (1/(1 + resistance)) - slocalguess
  
RETURN
END FUNCTION funslocal


!**********************************************************************
!     Parameterization for predicting the unhollowed basal face
!     ring width from environmental conditions and crystal dimensions.
!     Return: aring/a (1 - hollowing fraction)
!**********************************************************************
function hollowfrac(sicein, scritcin)

	implicit none

	real sicein, sice, scritcin, scritc, scrithol, A, B
	real hollowfrac
	
	! Scaling parameters for equation
	A = 0.315
	B = 1.537
	
	! Multiply supersats by 100 (allows parameters A and B to be closer to unity)
	sice = sicein * 100.
	scritc = scritcin * 100.
	scrithol = B * scritc

	! Can only have hollowing for supersat above scrithol
	if(sice.ge.scrithol)then
		hollowfrac = 1. - (A * scritc)/(sice - scrithol)
	else
		hollowfrac = 0.
	endif
	
	! Predicted ring width is larger than a; crystal is not hollow
	if(hollowfrac.lt.0.)then
		hollowfrac = 0.
	endif
  
RETURN
END FUNCTION hollowfrac


!**********************************************************************
!     This function calculates the enthalpy of sublimation for       
!     water as per Bohren and Albrecht (1999, pg 197, eq 5.64). 
!     Arguement is in kelvin and the returned result is in mks 
!     Coded: JYH Feb 2005    
!**********************************************************************
FUNCTION FLHS (tmp)
        
  cpv = 1850.                ! specific heat of vapor const p
  ci = 2106.                 ! specific heat of ice
  t0 = 273.15
  xlvs0 = 2.834e6            ! enthalpy of sublimation at 0 C 
  delcv = cpv - ci
  delt = t0 - tmp
  xlvs = (xlvs0 - delcv*delt) ! B&A equation 5.64
  FLHS = xlvs
  
  RETURN
END FUNCTION FLHS

!******************************************************************************
!******************************************************************************
!*
function FLHV (t)
!+
!+ purpose:  This function calculates the latent heat of vaporization for     +
!+           water as per P&K eqn. 4-85a. The result is then converted to     +
!+           ergs/g then to MKS.  T is deg K.                                 +
!+                                                                            +
  rga  = 0.167 + (3.67e-4 * t)
  rlhv  = 597.3 * (273.15 / t)**rga
  FLHV = rlhv / 2.38844e-08 /1.e4
                                                                           
  return
end function FLHV

!******************************************************************************
!******************************************************************************
!*
function capacitance (adim, cdim, shape, nfacets)
!+
!+ purpose:  This function calculates the electrostatic capacitance for       +
!+           faceted crystals of various shapes with basal semidimension a    +
!+           and prism semidimension c.                                       +
!+           Allowed shapes:                                                  +
!+		     1 - cylinder                                                     +
!+           2 - hexagonal prism                                              +
!+           3 - rosette                                                      +
!+           4 - plate-polycrystal (if I ever figure out how to implement it) +
!+       
 
  integer shape, nfacets  
  real adim, cdim, capacitance

  if(shape.eq.1)then
  	capacitance = 0.637*(1. + 0.868*((cdim/adim)**0.76))*adim
  elseif(shape.eq.2)then
  	capacitance = 0.58 * (1 + 0.95 * (cdim/adim)**0.75)*adim
  elseif(shape.eq.3)then
    if(nfacets.lt.5)then
        capacitance = 0.35 * max((cdim/adim), 1.) ** -0.27 * (sqrt((4. * cdim)**2. + (2.*adim)**2.))
    elseif(nfacets.lt.7)then
        capacitance = 0.40 * max((cdim/adim), 1.) ** -0.25 * (sqrt((4. * cdim)**2. + (2.*adim)**2.)) 
    else
        capacitance = 0.40 * ((REAL(nfacets)/6.) ** 0.257) * (cdim/adim) ** -0.25 * &
        (sqrt((4. * cdim)**2. + (2.*adim)**2.))
    endif
  elseif(shape.eq.4)then
  print*,"Hah, you thought I figured out the capacitance of polycrystals..."
  	capacitance = -1.
  else
  print*,'not a valid shape'
  	capacitance = -1.
  endif
                                                                           
  return
end function capacitance

!******************************************************************************
!*
!******************************************************************************
!*
function FLHM(tk)
!     This function calculates the latent heat of fusion as per P&K eqn. 3. 26.The result is then converted to ergs/g.t in Celsius

  t = tk - 273.15
  a0=79.7
  a1=-0.1200
  a2=-8.0481*1.e-2
  a3=-3.2376*1.e-3
  a4=-4.2553*1.e-5
  xlhm= a0+a1*t**1+a2*t**2+a3*t**3+a4*t**4 ! IT cal/g
  xlhm=xlhm*4.18            !J/g  1 cal(IT)=4.18 J
  xlhm=xlhm*1.e7            !ergs/g   ergs=1.e-7 J
      
  FLHM=xlhm
  return
end function FLHM

!**********************************************************************
!******************************************************************************
!*
      subroutine ZZBREN (f, errabs, errrel, a, b, maxfn)
!
!  modification to IMSL subroutine ZBREN to find root of an equation.
!  f = externally defined function of radius. F passes back
!      the saturation ratio for a soln drop of the radius passed in
!
!  errabs  = absolute error tolerance (real)
!  errrel  = relative error tolerance (real)
!  a       = lowest possible radius   (real)
!  b       = largest possible radius  (real)
!  itmax   = max number of iterations (integer)

        external F
        save

        t = errrel
        ic = 2
        s = a
        fa = f(s)
        s = b
        fb = f(s)
!                                  test for same sign

        if (fa*fb .gt. 0.0) then
        !write (10,*) ' ERROR IN ARGUMENTS TO ZZBREN'
        !write (10,*) ' fa, fb =',fa,fb
        !   stop
        end if
!

10      c = a
        fc = fa
        d = b - c
        e = d
20      if (ABS(fc) .lt. ABS(fb)) then
           a = b
           b = c
           c = a
           fa = fb
           fb = fc
           fc = fa
        end if

        tol = t*AMAX1(ABS(b),0.1)
        rm = (c-b)*0.5
!                                  test for first convergence criteria
        if (ABS(fb) .le. errabs) go to 9000
!                                  test for second convergence criteria
        if (ABS(c-b) .le. tol) go to 9000
!                                  check evaluation counter
        if (ic .ge. maxfn) then
!        write (10,*) 'failure to converge in zzbren'
!        write (10,*) ' f1, f2 =',c,b
           go to 9000
        end if
!                                  is bisection forced
        if (ABS(e) .lt. tol) then
           e = rm
           d = e
        else
           if (ABS(fa) .le. ABS(fb)) then
              e = rm
              d = e
           else
              s = fb/fa
              if (a .ne. c) then
!                                  inverse quadratic interpolation
                 q = fa/fc
                 r = fb/fc
                 rone = r - 1.0
                 p = s*((c-b)*q*(q-r)-(b-a)*rone)
                 q = (q-1.0)*rone*(s-1.0)
              else
!                                  linear interpolation
                 p = (c-b)*s
                 q = 1.0 - s
              end if
              if (p .gt. 0.0) q = -q
              if (p .lt. 0.0) p = -p
              s = e
              e = d
!                                  if abs(p/q).ge.75*abs(c-b) then
!                                     force bisection
              if (p+p .ge. 3.0*rm*q) then
!                                  bisection
                 e = rm
                 d = e
              else
!                                  if abs(p/q).ge..5*abs(s) then force
!                                     bisection. s = the value of p/q
!                                     on the step before the last one
                 if (p+p .ge. ABS(s*q)) then
!                                  bisection
                    e = rm
                    d = e
                 else
                    d = p/q
                 end if
              end if
           end if
        end if
!                                  increment b
        a = b
        fa = fb
        temp = d
        if (ABS(temp) .le. 0.5*tol) temp = SIGN(0.5*tol,rm)
        b = b + temp
        s = b
        fb = f(s)
        ic = ic + 1
        if (fb*fc .le. 0.0) go to 20
        go to 10
!
9000    a = c
!
        return
      end subroutine ZZBREN
!-----------------------------------------------------
!-----------------------------------
