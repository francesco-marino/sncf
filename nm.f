      subroutine nm
      implicit double precision (a-h,o-z)
      PARAMETER (nnp=10000)
      common/bunits/pi,hbc,amc2,hbdm
      common/bctrnm/i_plt_sym,i_plt_effmass
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
      common/b_nm_const/alpha
      COMMON/DER/FCT(NNP),DF(NNP),H       
      COMMON/BMAX/NMAXT
      external der_ea
    4 format(1x,5(e12.5,1x))

c-----Constants
      alpha = (3.d0/5.d0)*(24.d0*(datan(1.d0))**2)**(2.d0/3.d0)
      gamma = 1.d0/3.d0
      t1x = t1*(1+0.5d0*x1)
      t2x = t2*(1+0.5d0*x2)
      tx = t1x + t2x
      sumt12 = 0.25d0*(1.d0/hbdm)*(t1x+t2x)
      sumt12d = sumt12
      w1=2.d0*w2
      vso = -(9.d0/32.d0)*(1.d0/hbdm)*w2**2
c-----Calculation of rho_NM,E/A,radius and effective mass
      a = 0.05d0
      b = 2.95d0
c     ifail = 0
c     call c05adf(a,b,eps1,eps2,der_ea,rhonm,ifail)
c     if (ifail.eq.0) go to 5
c     write(*,*) '>>> ERROR IN C05ADF: ifail = ',ifail
c     stop
c   5 continue
      accur=1e-8
      call bsec(der_ea,a,b,accur,rhonm)
c     call bsec(f,xmin,xmax,accur,xzero)
      write(2,*)
      write(2,*) 'Nuclear matter density rho_nm ............ ',rhonm
      enm = densen_iso(rhonm) / rhonm
      write(2,*) 'Nuclear matter energy E/A ................ ',enm
      rnm = (3.d0/(4.d0*pi*rhonm))**(1.d0/3.d0)
      write(2,*) 'Radius r_nm .............................. ',rnm
      rapp = 1.d0/rapm(rhonm)
      write(2,*) 'm*/m ..................................... ',rapp
      eff_mass = rapp
      thetas = 3.d0*t1+t2*(5.d0+4.d0*x2)
      write(2,*) 'Theta_s .................................. ',thetas
      thetav = t1*(x1+2.d0)+t2*(x2+2.d0)
      write(2,*) 'Theta_v .................................. ',thetav
c-----Plot effective masses
      if(i_plt_effmass.eq.0)go to 50
      open(unit=97,status='unknown',file='mstar2.dat')
      open(unit=98,status='unknown',file='mstar.dat')
      beta_asym = 0.2d0
      rhomin=0.005d0
      rhomax=0.450d0
      drho=0.005d0
      nptrho=(rhomax-rhomin)/drho+1
      do 51 irho=1,nptrho
      rho=rhomin+float(irho-1)*drho
      ratio_n = 1.d0 + 1.d0/(16.d0*hbdm)*(thetas -
     &(2*thetav-thetas)*beta_asym)*rho
      ratio_n = 1.d0 / ratio_n
      ratio_p = 1.d0 + 1.d0/(16.d0*hbdm)*(thetas +
     &(2*thetav-thetas)*beta_asym)*rho
      ratio_p = 1.d0 / ratio_p
      write(98,4) rho,ratio_n,ratio_p
   51 continue
      rho=0.18d0
      do 52 i_asym=1,10
      beta_asym=float(i_asym-1)*0.1
      ratio_n = 1.d0 + 1.d0/(16.d0*hbdm)*(thetas -
     &(2*thetav-thetas)*beta_asym)*rho
      ratio_n = 1.d0 / ratio_n
      ratio_p = 1.d0 + 1.d0/(16.d0*hbdm)*(thetas +
     &(2*thetav-thetas)*beta_asym)*rho
      ratio_p = 1.d0 / ratio_p
      write(97,4) beta_asym,ratio_n,ratio_p
   52 continue
      close(98)
      close(97)
c-----Landau parameters
c-----Coefficients: d,d1,d2,d3,J,K,L,K_sym,K_tau,K_C
   50 continue
      write(2,*)
      d = (1.d0/64.d0)*(9.d0*t1-t2*(5.d0+4.d0*x2))+((beta-gamma)/16.d0)*
     1    (3.d0*t1+t2*(5.d0+4.d0*x2))
      d1 = (3.d0/32.d0)*(t1*(1.d0-x1)-t2*(1.d0+x2))+(1.d0/8.d0)*
     1     (3.d0*t2*(1.d0+x2)+t1*(1.d0-x1))*(beta-gamma)
      d2 = (1.d0/8.d0)*(3.d0*t1*(1.d0+x1/2.d0)-t2*(1.d0+x2/2.d0))-
     1     (gamma/2.d0)*(t1*(1.d0+x1/2.d0)+t2*(1.d0+x2/2.d0))
      d3 = (1.d0/4.d0)*(t1*(1.d0+x1/2.d0)+t2*(1+x2/2.d0))*beta
      cj = (5.d0/9.d0)*hbdm*alpha*rhonm**(2.d0/3.d0)-(t0/4.d0)*(x0+
     1     0.5d0)*rhonm-(t3/24.d0)*(x3+0.5)*rhonm**(1.d0+alfe)+
     2     (5.d0/72.d0)*(t2*(4.d0+5.d0*x2)-3.d0*t1*x1)*alpha*
     3     rhonm**(5.d0/3.d0) 
      write(2,*) 'Volume symmetry energy J ................. ',cj
      ck = -2.d0*hbdm*alpha*rhonm**(2.d0/3.d0)*(1.d0-(5.d0/16.d0)*
     1     (1.d0/hbdm)*(3.d0*t1+t2*(5.d0+4.d0*x2))*rhonm)+
     2     (9.d0/16.d0)*(1.d0+alfe)*alfe*t3*rhonm**(1.d0+alfe) 
      write(2,*) 'Nuclear matter incompressibility K_infty . ',ck
      cl = (10.d0/9.d0)*hbdm*alpha*rhonm**(2.d0/3.d0)-(3.d0*t0/4.d0)*
     1     (x0+0.5)*rhonm-(1.d0+alfe)*(t3/8.d0)*(x3+0.5)*
     2     rhonm**(1.d0+alfe)+(25.d0/72.d0)*(t2*(4.d0+5.d0*x2)-
     3     3.d0*t1*x1)*alpha*rhonm**(5.d0/3.d0)
      write(2,*) 'Coefficient L ............................ ',cl
C_TEST_22_08_13
      write(94,*) cl
      cksym = (5.d0/3.d0)*alpha*rhonm**(2.d0/3.d0)*((-2.d0/3.d0)*
     1        hbdm+(5.d0/12.d0)*rhonm*(t2*(4.d0+5.d0*x2)-3.d0*t1*x1))
     2        -(3.d0/8.d0)*t3*(0.5d0+x3)*alfe*(alfe+1.d0)*
     3        rhonm**(alfe+1.d0)
      write(2,*) 'Coefficient K_sym ........................ ',cksym
      d3h = (30.d0/72.d0)*alpha*((3.d0*t1/2.d0)+2.d0*t2*(x2+5.d0/4.d0))*
     1      rhonm**(-1.d0/3.d0) - alpha*hbdm*(1+(rhonm/hbdm/8.d0)*
     2      ((3.d0*t1/2.d0)+2.d0*t2*(x2+5.d0/4.d0)))*(10.d0/27.d0)*
     3      rhonm**(-4.d0/3.d0) + t3/16.d0*(alfe+2)*(alfe+1)*alfe*
     4      rhonm**(alfe-1)
      r2d3h = rhonm**2*d3h
      write(2,*) 'rho_nm**2 * (d^3 h/d rho^3) .............. ',r2d3h
      cktau = cksym + 3.d0*cl - 27.d0*cl*r2d3h/ck
      write(2,*) 'Coefficient K_tau ........................ ',cktau
      ckc = 0.764d0 * (1.d0 - 27.d0*r2d3h/ck)
      write(2,*) 'Coefficient K_Coul ....................... ',ckc
      akf=(3.d0*pi**2*rhonm/2.d0)**(1./3.)
      write(2,*) 'k_F ...................................... ',akf
      an0=2.d0*akf*rapp*amc2/pi**2/hbc**2
      write(2,*) 'N0 ....................................... ',an0
      write(2,*) '1/N0 ..................................... ',1.d0/an0
      f0=0.75d0*t0+t3/16.d0*rhonm**alfe*(alfe+1.d0)*(alfe+2.d0)
     &+1.d0/8.d0*akf**2*(3.d0*t1+t2*(5.d0+4.d0*x2))
      f0=an0*f0
      write(2,*) 'f0 ....................................... ',f0
      f0p=0.25d0*t0*(1.d0+2.d0*x0)
     &+t3/24.d0*rhonm**alfe*(1.d0+2.d0*x3)
     &+1.d0/8.d0*akf**2*(t1*(1.d0+2.d0*x1)-t2*(1.d0+2.d0*x2))
      f0p=-an0*f0p
      write(2,*) 'f0p ...................................... ',f0p
      g0p=0.25d0*t0+t3/24.d0*rhonm**alfe+1.d0/8.d0*akf**2*(t1-t2)
      g0p=-an0*g0p
      write(2,*) 'g0p ...................................... ',g0p
c-----Parameter F_D of PLB 363 (1995) 5 for 208Pb
      f1d=tx/694.48d0
      f2d=5.d0/3.d0*(-5.599d0*rapp+6.899d0)*0.168777d0
      fd=sqrt(f1d*cj*f2d)
      write(2,*)
      write(2,*) 'Parameter F_D of PLB 363 (1995) 5 ........ ',fd
      ed=11.352d0+0.771d0*fd
      write(2,*) 'Predicted E(-1) for the IVGDR in 208Pb ... ',ed
c-----Symmetry energy   
      sym0=sym(rhonm)
      write(2,*)
      write(2,*) 'J calculated from the explicit S(rho) .... ',sym0
      sym1=sym(0.1d0)
      write(2,*) 'Value of the symmetry energy at 0.1 fm^-3. ',sym1
      a = 208.d0
      a3 = a**(-1.d0/3.d0) 
      gA16 = cj/(1.d0+a3*(5.d0/cj)*((cl/3.d0)-(cksym/36.d0)))
      write(2,*) 'g_A(0.16) as defined in paper with L.T. .. ',gA16
c     if(i_plt_sym.ne.1)go to 1
      open(unit=99,status='unknown',file='sym.dat')
      rhomin=0.005d0
      rhomax=0.350d0
      drho=0.001d0
      nptrho=(rhomax-rhomin)/drho+1
      if(nptrho.gt.10000) stop '>>> TOO MANY POINTS FOR THE S(rho)'
      do 2 irho=1,nptrho
      rho=rhomin+float(irho-1)*drho
      sss=sym(rho)
      eee=ea(rho)
      enm=eee+sss
cc    HERE   writing energies    ->   sss: symmetry energy
cc    eee: symmetric matter;   enm: neutron matter
      write(99,4) rho,sss,eee,enm,sss/cj
      fct(irho)=sss
      if(rho.lt.rhonm) isign=-1
      if(rho.gt.rhonm) isign=1
      if(irho.eq.1)go to 3
      if((isign*isign_old).lt.0)then
       iinf=irho-1
       isup=irho
      end if 
    3 isign_old=isign
    2 continue
      h=drho
      nmaxt=nptrho
      call deriv
c     df(irho) contains the derivative of S(rho)
      sym1=df(iinf)+(df(isup)-df(iinf))*
     1(rhonm-float(iinf-1)*drho-rhomin)/drho
      sym1=3.d0*rhonm*sym1
      write(2,*) 'L calculated from the explicit S(rho) .... ',sym1
      close(99)
    1 continue
     
      return
      end

      function densen_iso(rho)
      implicit double precision (a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
      common/b_nm_const/alpha
ccc   h(rho) defined in (A.5b)
      densen_iso = alpha*hbdm*rapm(rho)*rho**(5.d0/3.d0) +
     1             (3.d0/8.d0)*t0*rho**2 + (t3/16.d0)*rho**(2.d0+alfe)
      return
      end

      function ea(rho)
      implicit double precision (a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
      common/b_nm_const/alpha
      thetas = 3.d0*t1+t2*(5.d0+4.d0*x2)
      ea = alpha*hbdm*rho**(2.d0/3.d0)
     1   + (3.d0/8.d0)*t0*rho 
     2   + alpha*1.d0/16.d0*
     3     (3.d0*t1+t2*(5.d0+4.d0*x2))*rho**(5.d0/3.d0)
     4   + (t3/16.d0)*rho**(1+alfe)
      return
      end

      function der_ea(rho)
      implicit double precision (a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
      common/b_nm_const/alpha
ccc   The derivative of h(rho)/rho, with h(rho) defined in (A.5b)
ccc   and calculated in densen_iso
      der_rapm = (1.d0/hbdm)*((3.d0/2.d0)*t1+2.d0*t2*
     1           ((5.d0/4.d0)+x2))*(1.d0/8.d0)
c     der_ea   = (2.d0/3.d0)*alpha*hbdm*rapm(rho)*rho**(-1.d0/3.d0)
c    1                 + (3.d0/8.d0)*t0 +
c    2                 (t3/16.d0)*(1.d0+alfe)*rho**alfe
c    3                 + alpha*hbdm*der_rapm*rho**(2.d0/3.d0)
      thetas = 3.d0*t1+t2*(5.d0+4.d0*x2)
      der_ea = (2.d0/3.d0)*alpha*hbdm*rho**(-1.d0/3.d0)
     1                 + (3.d0/8.d0)*t0 
     2                 + alpha*5.d0/48.d0*
     3                 (3.d0*t1+t2*(5.d0+4.d0*x2))*rho**(2.d0/3.d0)
     4                 + (t3/16.d0)*(1.d0+alfe)*rho**alfe
      return
      end
 
      function sym(rho)
      implicit double precision (a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
      common/b_nm_const/alpha
      s1 = alpha*hbdm*5.d0/9.d0*rho**(2.d0/3.d0)
      s2 = -t0/8.d0*(2.d0*x0+1.d0)*rho
      s3 = -t3/48.d0*(2.d0*x3+1.d0)*rho**(alfe+1)
      s4 = alpha/16.d0*( (t1*(x1+2.d0)+t2*(x2+2.d0))*10.d0 +
     1     0.5d0*(t2*(2.d0*x2+1.d0)-t1*(2.d0*x1+1.d0))*40.d0 ) * 
     1     rho**(5.d0/3.d0)/9.d0
      sym = s1 + s2 + s3 + s4
      return
      end

      function rapm(rho)
      implicit double precision (a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
ccc   rapm, that is, m/m*
      rapm = 1.d0 + (1.d0/hbdm)*((3.d0/2.d0)*t1+2.d0*t2*
     1((5.d0/4.d0)+x2))*(rho/8.d0)
      return
      end

