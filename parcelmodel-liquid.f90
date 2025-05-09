program parcelmodel

  implicit none
  
  ! stuff for vode

  integer neq,neqtot,iequ,ic,nbliq
  integer itol,itask,istate,iopt,lrw,liw,mf
  real tc,tout,t
  parameter (neq=5,nbliq=100,neqtot=nbliq+neq)
  integer ipar(20*neqtot)
  real y(neqtot),dy(neqtot),iwork(40),rwork(22+16*neqtot)  ! y = radius,  dy = dr/dt growth rate
  real rpar(20*neqtot)
  real atol(neqtot)
  real rtol
  external jac,FunDiff,vode

  ! other variables
  integer k,ivap,ipress,itemp,iliq,iice,nt
  integer iwvel,isolute,i

       !rhodrop, ml, LWC, MRl, dMRl, initdistflagliq, Nlo, & 
       !NUl, isat,msaltflag,nlinit,rdrp,rdrp0,rdryaer,sigmaccn
  real radiusin,Nlo,Nui,ravgi,ravgl,sigmaccn
  real rlo(nbliq), Nmrl(nbliq),rhodrop(nbliq),Nl(nbliq),msalt(nbliq)
  real rdryaer(nbliq),rhosolute,molmass,vhoff,MRnucdep,MRnucfreez
  real nlinit(nbliq),nhetfrz(nbliq),rdrp(nbliq),rdrp0(nbliq)
  real Nltot,NUl,lwc,ml,mrl,dmrl,es,rdry,msalt2,temp,Sl
  real Nctot,ravg,rcavg,Ntot
  
  real Tinit,TinitK,Pinit,ztopinit,RHinit
  real pi,Rd,Rv,cp,grav,dt,To,w,rhoa
  real tempk,press,wmax,zalt,vapmr,nosc,tau_osc
  real qlsat,RH
  real esi_new,es_new

  logical endloop

  COMMON/neqtotvars/itemp,ipress,ivap,iliq,iice
  COMMON/envvars/w,rhodrop, Nl,rdrp,rdrp0,rdryaer,Nmrl
  COMMON/ccnvars/msalt
  COMMON/aerosol/Sl, temp, es, rdry, msalt2
  common/ccnchar/vhoff,molmass,rhosolute
  
  
  parameter(iwvel=1)      ! 1 = cosine w profile
                          ! 2 = constant w
  parameter(isolute=2)    ! 1 = NaCl (sodium chloride)
                          ! 2 = (NH4)2SO4 (ammonium sulfate)


  ! output file

  open(20,file='outdir/parceltime.dat')
  
  ! physical constants and initial conditions
  pi = 3.14159   
  Rd = 287.                 ! dry air gas constant
  Rv = 461.5    ! water vapor gas constant
  cp = 1004.5
  grav = 9.81
  dt = 1.0                 ! time step
  To = 273.15               ! [=]K
 
  Tinit = -30.              ! initial temperature [C]
  TinitK = Tinit + To
  Pinit = 40000.0           ! initial pressure [Pa]
  ztopinit = 2000.0         ! Total cloud depth [m]
  rhinit = 0.9!*esi_new(TinitK)/es_new(TinitK)
  wmax = 0.5   ! maximum vertical motion
  endloop = .false.

  !     initial conditions for microphysics
  radiusin = 1.e-6         ! radius of heterogeneous nuclei
  Nlo = 200.*100.**3 !20.*100.**3          ! total # of dry CCN per volume[=]#/m^3 
  NUi = 4.                  ! determines ice distribution shape
  ravgi = 10.e-6            ! initial average radius of ice particles (m)
  ravgl = 0.02e-6           ! geometric radius of dry CCN distribution (m)
  sigmaccn = 2.2            ! geometric standard deviation of CCN distribution


  ! parameters for aerosol
  if(isolute.eq.1)then
     vhoff = 2.0            !van't hoff factor
     Molmass = 0.058443     ! kg NaCl
     rhosolute = 2165.      !kg/m^3
  else
     vhoff = 2.4            !van't hoff factor ! Modify Later
     Molmass = 132.13/1000. ! molecular mass ammonium sulfate (NH4)2SO4
     rhosolute = 1769.      !kg/m^3
  endif

  ! oscillation scales
  tau_osc = ztopinit/(2.*wmax)               ! Parcel oscillation time-scale
  nosc = 1.                              ! Number of 1/2 parcel oscilations
  nt = int(nosc*tau_osc/dt*pi)              ! Total number of time-steps
  if(iwvel.eq.2)then
     nt = int(nosc*ztopinit/wmax/dt)         ! if using a constant vertical motion
  endif

  itemp = nbliq + 1
  ipress =  nbliq + 2
  ivap =  nbliq + 3
  iliq =  nbliq + 4
  iice =  nbliq + 5
  
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

  ! initial conditions

  tempk = TinitK
  press = Pinit
  zalt = 0.0
  vapmr = (es_new(tempk)*(rhinit)*Rd)/(Rv*press)
  y(ivap) = vapmr
  y(iliq) = 0.0
  y(iice) = 0.0
  y(ipress) = Pinit
  y(itemp) = TinitK
  rhoa = press/(Rd*tempk)

  ! form solution drops from dry aerosol distribution
  ! assumes NaCl for now. Solution drops assumed to be
  ! in equilibrium with the vapor field initially.

  Sl = rhinit
  temp = tempk
  es = es_new(tempk)
  CALL LiqNucleation(rlo, y, Nl, Nmrl, rhoa, ravgl, &
       rhodrop, ml, LWC, MRl, dMRl, Nlo, & 
       NUl, rhinit,nlinit,rdrp,rdrp0,rdryaer,sigmaccn)

  ! time integration loop
  tc = 0.0
  k = 0
  Do while (tc.lt.float(nt)*dt.and..not.endloop)

     w= wmax*sin(tc/tau_osc)
     if(iwvel.eq.2) w = wmax

     print*,y(100),y(99)
     t = tc
     tout = t + dt
     istate = 1
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

     zalt = zalt + w*dt
     qlsat = es_new(y(itemp))*Rd/(Rv*y(ipress))
     RH = y(ivap)/qlsat * 100.
     print*,w,RH,ravg,rcavg
     write(20,'(45E16.8)')tc,zalt,y(itemp)-To,y(ipress)/100. &
          ,y(ivap)*1000.,y(iliq)*1000.,y(iice)*1000.,RH, &
          Ntot/100.**3,Nctot/100.**3,ravg,rcavg
     
     k = k + int(dt)
     tc = tc + dt
  enddo
     
  stop
end program parcelmodel

SUBROUTINE FunDiff(neqtot, t, y, dy, rpar, ipar)
       
  IMPLICIT NONE
  
  integer neqtot,itemp,ipress,ivap,iliq,iice,nbliq
  PARAMETER(nbliq = 100)
  real y(neqtot),dy(neqtot),t,w
  real rhodrop(nbliq),nl(nbliq),rdrp(nbliq),rdrp0(nbliq)
  real Nmrl(nbliq),msalt(nbliq),vhoff,molmass,rhosolute
  real rdryaer(nbliq)
  integer ipar(*)
  real rpar(*)

  COMMON/neqtotvars/itemp,ipress,ivap,iliq,iice
  COMMON/envvars/w,rhodrop, Nl,rdrp,rdrp0,rdryaer,Nmrl
  COMMON/ccnvars/msalt
  common/ccnchar/vhoff,molmass,rhosolute

  print*,'in', y(100),y(99)
  CALL liqGrowth(neqtot, y, dy, rhodrop, itemp, &
       ipress, ivap, Nl,msalt,rdrp,rdrp0, &
       rdryaer)
  
  CALL EnvEqns(neqtot,y,dy,w,itemp,ipress,ivap,iliq,iice, &
       nbliq,Nmrl,rhodrop)
       
  return
end SUBROUTINE FunDiff

SUBROUTINE liqGrowth(neqtot, y, dy, &
     rhodrop, itemp, ipress, ivap, Nl,msalt, &
     rdrp,rdrp0,rdryaer)

  IMPLICIT NONE
  
  INTEGER i, neqtot, nbice, itemp, ipress, isat, ivap, il, nbliq, &
       nbicelim
  PARAMETER(nbliq = 100)
  REAL Ls, Rv, ei, es, Dv, Kt, Lv, Sl, G, Nl(nbliq), &
       temp, To, dt, eo, pi, gam, P, Po, y(neqtot),t,rhodrop(nbliq), &
       es_new, esi_new, m, celsius, drho,IGR_DEN, dy(neqtot), &
       A, B, Mw, vhf, Ms, R, sfcTens, rhol, msalt(nbliq), &
       rdrp(nbliq),rdrp0(nbliq),fsoln,Slavg,dtime,aw,vdrop,nliq, &
       nsol,rhosalt,vsol,vliq,rhosoldrop,tempk,rdryaer(nbliq),epsr, &
       vhoff,molmass,rhosolute,qlsat,Rd
  real FLHV
  common/ccnchar/vhoff,molmass,rhosolute

  dtime = 1.0
  To = 273.15               ! [=]K
  !     rho = 920.                ! [=]kg/m^3
  Rd = 287.0
  Rv = 461.5                ! Individual gas constant for water vapor[=]J/kgK
  Lv = 2.5E6                ! Latent heat of vaporization[=]J/kg
  pi = 3.14159
  Mw = 0.018015             !kg
  R = 8.314                 ! J/molK
  sfcTens = 7.5e-2          !N/m
  rhol = 1000.0   
  temp = y(itemp) 
  Lv = FLHV(temp)

  Ms = molmass
  rhosalt = rhosolute
  vhf = vhoff
  
  P = y(ipress)
  Po = P/100.
      
  DO i = 1, nbliq

     il = i 
!         print*,'testing rdrp0,',i,rdrp0(i)*1.e6,y(il)*1.e6,
!     +        rdryaer(i)*1.e6
!         IF(y(il) .lt. 9.9E-10) y(il) = 0.0 

         !IF (y(il) .lt. 2.E-8 .or. Nl(i) .eq. 0.0) THEN
     epsr = rdryaer(i)*0.1
     If (y(i) .lt. rdryaer(i)+epsr .or. Nl(i) .eq. 0.0) THEN
        
        dy(i) = 0.0
            !y(il) = rdrp0(i)

     ELSE

        es = es_new(temp)         
        
        ei = esi_new(temp)

            !rhodrop(i) = 1000.0
        qlsat = es*Rd/(Rv*y(ipress))
        Sl = y(ivap)/qlsat !*(ei/es)
        
        Dv = 0.211*(temp/273.15)**1.94 * (1013.25/Po)* 1/100.**2 ! vapor diffusivity 
        
        Kt = (5.69 + 0.017*(temp-273.15)) * 418.684/1.e5 ! thermal diffusivity
        
        G = ((((Lv/(Rv*temp))-1.)*(Lv/(Kt*temp))) + ((Rv*temp) &
             /(Dv*es)))**(-1.)
        
        A = (2*Mw*sfcTens)/(R*temp*rhol)

        B = (3*vhf*msalt(i)*Mw)/(4*pi*Ms*rhol)
                     
        call getdropchar(y(i),msalt(i),vsol,vliq,aw,rhosoldrop)
!            print*,'testing again ',i,y(il),msalt(i),aw,aw*exp(A/y(il))
        rhodrop(i) = rhosoldrop
        !dy(il) = (G*(Sl - 1.- ((A/y(il)) - (B/y(il)**3))))/ &
        !     (y(il)*rhodrop(i))

        fsoln = 1.-aw*exp(A/y(i))
        !fsoln = 0.0
        dy(i) = (G*(Sl - 1. + fsoln))/(y(i)*rhodrop(i))
        !print*,i,y(i),dy(i),aw,Sl-1.,fsoln
        !dy(i) = 0.0
        if(vliq.lt.0.0)then
           print*,'inside liqgrowth vsol > vdrop'
           print*,il,i,msalt(i),vsol,vliq,aw,y(i),rdrp0(i)
           dy(i) = 0.0               
           stop
        endif

!            if(i.eq.24)then
!               print*,'inside liqgrowth vsol > vdrop'
!               print*,il,i,msalt(i-1),msalt(i),msalt(i+1),
!     +              vsol,vliq,aw,y(il),rdrp0(i)
!            endif
            !print*,'testing ',i,rhosoldrop
        IF( i .eq. 50)THEN
           !            WRITE(*,*) 'A', Sl - 1.,(A/y(il)),(B/y(il)**3)
        END IF
     END IF
         !WRITE(*,*) 'before', dy(il), y(il)

  END DO
      !stop
END SUBROUTINE liqGrowth

subroutine EnvEqns(neqtot,y,dy,w,itemp,ipress,ivap,iliq,iice, &
     nbliq,Nmrl,rhodrop)

  implicit none

  integer neqtot,itemp,ipress,ivap,iliq,iice,nbliq,i
  real w
  real y(neqtot),dy(neqtot)
  real g,Rd,Mw,Ma,eps,Rm,Cp,Lv,Ls,Rv,qlsat
  real Nmrl(nbliq),rhodrop(nbliq),pi
  real FLHV,FLHS,es_new

  g = 9.81
  pi = 3.1415926
  Rd = 287.0
  Rv = 461.5
  Mw = 0.01801528	! Molar mass of water [kg/mole]
  Ma = 0.0289647	! Molar mass of dry air [kg/mole]
  eps = Mw/Ma
  Rm = (1.0 + (y(ivap)/eps))/(1.0 + y(ivap)) * Rd 		! Gas constant for moist air.
  Cp = 1004.5
  Lv = FLHV(y(itemp))
  Ls = FLHS(y(itemp))

  qlsat = es_new(y(itemp))*Rd/(Rv*y(ipress))
  
  dy(iliq) = 0.0
  !if (qlsat .lt. y(ivap))then
  !   dy(iliq) = (y(ivap) - qlsat)*0.5/2.0
  !endif
  do i = 1,nbliq
     dy(iliq) = dy(iliq) + 4.*pi*y(i)**2 * rhodrop(i) * Nmrl(i)*dy(i)
     !print*,'nmrl = ',nmrl(i),dy(i),y(i)
  enddo
  dy(iice) = 0.0
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
  PARAMETER(nbice = 1000, nbliq = 100, ne = nbliq + 1, &
       neqtot = nbliq + 5)
  REAL rlo(nbliq), y(neqtot), Nl(nbliq),Nlo, NUl, rnl, ravgl, &
       Nmrl(nbliq), rhoa, rhodrop(nbliq), ml, pi, LWC, MRl, dMRl, & 
       re(ne), Sl, temp, msalt(nbliq), es, sigma, rlow, &
       rhigh, errabs, rel, rdry, msalt2, sig,nlinit(nbliq), &
       rdrp(nbliq),rdrp0(nbliq),vsol,vliq,aw,rhosoldrop, &
       rdryaer(nbliq),vhoff,molmass,rhosolute,sigmaccn,rhinit
  COMMON/aerosol/Sl, temp, es, rdry, msalt2
  COMMON/ccnvars/msalt
  common/ccnchar/vhoff,molmass,rhosolute

  EXTERNAL Kohler

  Sl = rhinit
  pi = 3.141529
  sig = sigmaccn  ! 2.2 !1.8  ! geometric standard deviation 1.5 to 2.3
  !sig = 1.5
  sigma = log(sig)
  NUl = exp((sigma**2)/2.)
  rnl = ravgl   !ravgl/NUl
  open(51, file = "outdir/initialdrydist.dat")
  OPEN(52, file = "outdir/initialwetdist.dat")

  CALL RadEdgeLiq(Nlo, NUl, rnl, nbliq, ne, re, Nl, rlo, sig, ravgl) !Get dry aerosol distribution

 !     write(0,*)'\n\n',re,'\n\n',rlo
  DO i = 1, ne-1
     WRITE(51,'(45E16.8)') re(i),Nl(i), rlo(i)
     !         write(0,*)4./3.*pi*re(i)**3,4./3.*pi*rlo(i)**3
  END DO

  CLOSE(51)

!      write(0,*)'\n\n'

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
            !print*,'msalt is ',j,msalt(j),rdrp(j)*1.e6,rel
      
     Nmrl(j) = Nl(j)/rhoa
     nlinit(j) = nmrl(j)         
     rhodrop(j) = rhosoldrop
     ml = (4./3.)*pi*(y(j)**3)*rhodrop(j)
     LWC = LWC + Nl(j)*ml
  END DO
  DO j = 1,nbliq
     WRITE(52,'(45E16.8)') Nl(j), y(j), msalt(j)
     !         write(0,*)(4./3.*pi*y(jl)**3-msalt(j)/2165.)*1000./0.018,
!     1        msalt(j)/0.0584
  END DO

  MRl = LWC/rhoa
  
  !stop
  
END SUBROUTINE LiqNucleation

!*********************************************************************
!     This subroutine will calculate and return the mass of the salt,
!     msalt, with a given relative humidity and radius of the liquid
!     drop.
!**********************************************************************

REAL FUNCTION Kohler(rdrop)
      
  IMPLICIT NONE

  REAL es, A, ravgl, B, Mw, sfcTens, R, temp, rhol,vhf, &
       Ms,es_new, esol, pi, Sl, rdry, rhosalt, rdrop, Ssol, diff, &
       msalt2,vhoff,molmass,rhosolute,vsol,vliq,aw,rhosoldrop

  COMMON/aerosol/Sl, temp, es, rdry, msalt2
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

  A = (2.0*Mw*sfcTens)/(R*temp*rhol) !curvature term
  
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
      
SUBROUTINE RadEdgeLiq(Nio, NU, rn, nb, ne, re, Ni, rio,sig,rg)

  IMPLICIT NONE
      
  INTEGER i, nb, ne 
  REAL Nio, NU, rn, re(ne), rlow, rhigh,deltar,Ni(nb),rio(nb)
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
     rio(i) = (re(i+1) + re(i))/2. ! size of particles in each bin
     Ni(i) = Nio/( (2.*pi)**0.5 * alog(sig)*rio(i) ) * &
          exp( -(alog(rio(i)/rg))**2/(2.*(alog(sig))**2) ) * &
          (re(i+1) - re(i))
     ntotal = ntotal + Ni(i)

  END DO
      
                                !    print*,'testing gamm = ',gammtest
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
