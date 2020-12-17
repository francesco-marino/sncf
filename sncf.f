c     Program Skyrme Consistent Nuclear Field                  |
c     G. Colo`, 2003                                           |
c     ----------------------------------------------------------
      implicit double precision(a-h,o-z)
      common/b_to_do/ichf,ihf,irpa
      common/btime/time1,sumtime
      common/b_gr0/igr0_j
 
      open(unit=1,status='old',file='sncf.in')
      open(unit=2,status='unknown',file='sncf.out')
      call units
      call reader
      if(ichf.eq.1)then
       call cbcs_new 
       stop
      end if
      if(ihf.eq.1)then
       call bcsgen
      end if       
      if(irpa.ne.0.and.igr0_j.eq.0)then
       call qrpa
      end if       
      if(igr0_j.ne.0)then
       stop
      end if       
      close(2)
      close(1)
      stop
      end 



      subroutine units
      implicit double precision(a-h,o-z)
      common/bunits/pi,hbc,amc2,hbdm
      pi=4.d0*datan(1.d0)
      hbc=197.327053d0            ! hbar
      amc2=938.91869d0            ! mass
c     obtained as (m_p + m_n)/2.
      hbdm=hbc**2/2.d0/amc2
      return
      end



      subroutine reader
      implicit double precision(a-h,o-z)
      character purpose*80,skyrme*5,shells*10,purpose2*80
      character spurious_proc*80,potential*80
      character(3) :: dummy
      integer bcsp,bcsn

      common/b_to_do/ichf,ihf,irpa
      common/bmax/nmaxt
      common/b_box/ibox
      common/b_chf_1/irpa_chf,n_constr_step
      common/b_chf_2/constr_step
      common/b_mesh_r/step
      common/b_forceparam/skyrme
      common/b_imore_n/i_more_n
      common/b_occup/nos,npo,ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
     &ioc04,ideg4,ioc05,ideg5,ioc06,ideg6 
      common/b_potential/icoul_d,icoul_e,iso,isj,icm1
      common/b_iter/nhf,niter
      common/b_extra/npplus,nnplus
      common/b_bcs/bcsp,bcsn
      common/b_vpair/vzero,xx,densc,gamma
      common/bcons/icons,iovert
      common/bctrnm/i_plt_sym,i_plt_effmass
      common/bdeb/ideb0,ideb
      COMMON/CTR2/IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON/PAR/GPP,GPH,DGPP,DGPH
      COMMON/PAR1/NGPP,NGPH
      common/hf/nmax,nocc,nunocc,norb
      common/b_isotope/aa,zz
      common/bafa0/afa0
      common/bnosc/nosc
      common/becutp/ecutp_p,ecutp_n
      common/brc/e_centr_min,e_centr_max
      common/b_gr0/igr0_j 
      common/b_tensor/stc_t,stc_u
      common/b_tensor_i/itens

    1 format(75('*'))
      inm=0
      ihf=0
      ichf=0
      read(1,err=1000,end=1000,fmt='(a80)') purpose
      if(index(purpose,'HF').ne.0) ihf=1  
      if(index(purpose,'CHF').ne.0) ichf=1  
      read(1,err=1000,end=1000,fmt='(a5)') skyrme
      call forceparam
cc    infinite matter
      if(index(purpose,'NM').ne.0) inm=1
      if(inm.eq.1) then
       i_plt_sym=0
       i_plt_effmass=0
       read(1,err=1000,end=3,fmt='(a80)') purpose2
       if(index(purpose2,'PLOT SYMMETRY ENERGY').ne.0) i_plt_sym=1 
       if(index(purpose2,'PLOT EFFECTIVE MASS').ne.0) i_plt_effmass=1 
    3  call nm
       stop
      end if

cc    HERE   reading input file
      read(1,*,err=1000,end=1000) nmaxt,step    ! max radius, step
c     read(1,*,err=1000,end=1000) nos,npo
c     read A and Z
      read(1,*,err=1000,end=1000) aa,zz
      if (command_argument_count().ge.2) then   ! read A and Z from command line
        call get_command_argument(1, dummy); read(dummy,*) aa    !convert ro real
        call get_command_argument(2, dummy); read(dummy,*) zz        
      end if
cccc      print *, "AA=",aa,"  ZZ=",zz

      read(1,*,err=1000,end=1000) nhf,niter      !  ???, max. n. iter
      ioc01=0
      ideg1=0
      ioc02=0
      ideg2=0
      ioc03=0
      ideg3=0
      ioc04=0
      ideg4=0
      ioc05=0
      ideg5=0
      ioc06=0
      ideg6=0

      read(1,*,err=1000,end=1000) shells
      if(index(shells,'OPEN').ne.0)then
       read(1,*,err=1000,end=1000) ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
     & ioc04,ideg4,ioc05,ideg5,ioc06,ideg6
c      write(2,*) ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
c    & ioc04,ideg4,ioc05,ideg5,ioc06,ideg6
      end if

cc    set to 1 -> use BCS      
      read(1,*,err=1000,end=1000) bcsp,bcsn
cc    npplus, nnplus: add s.p. levels to HF calculation
      read(1,*,err=1000,end=1000) npplus,nnplus,ecutp_p,ecutp_n
cc    V(pair) = V0 [ 1 + x0 (rho/rho_0)^gamma ] 
      read(1,*,err=1000,end=1000) vzero,xx,densc,gamma 
      
      icoul_d=1
      icoul_e=1
      iso=1
      isj=0
      icm1=1
c     HERE reading options
      read(1,err=1000,end=1000,fmt='(a80)') potential
      if(index(potential,'NO COULOMB').ne.0) then 
      icoul_d=0
      icoul_e=0
      end if
      if(index(potential,'NO COULOMB EXCHANGE').ne.0) then
      icoul_d=1
      icoul_e=0
      end if
      if(index(potential,'NO SPIN-ORBIT').ne.0) iso=0
      if(index(potential,'INCLUDE J2 TERMS').ne.0) isj=1
      if(index(potential,'NO CENTER OF MASS CORRECTION').ne.0) icm1=0
      itens=0
      if(index(potential,'INCLUDE TENSOR').ne.0) itens=1
      if(itens.eq.1)read(1,*,err=1000,end=1000) t_alpha,t_beta
      write(2,*) 'isj,itens=',isj,itens
      stc_t=12.d0*(2.d0*t_beta-t_alpha)/5.d0
      stc_u=12.d0*t_alpha/5.d0
      if(ichf.ne.1)go to 2
      read(1,*,err=1000,end=1000) ideb0
      read(1,*,err=1000,end=1000) constr_step,n_constr_step,irpa_chf
    2 read(1,*,err=1000,end=1001) purpose2
      irpa=0
      igr0_j=0
      if(index(purpose2,'TDA').ne.0) irpa=1  
      if(index(purpose2,'RPA').ne.0) irpa=2  
      if(index(purpose2,'UNP').ne.0) irpa=3  
      if(index(purpose2,'GR0').ne.0) igr0_j=1
      if(index(purpose2,'HO').ne.0) then
       read(1,*,err=1000) afa0,nosc
      end if
      read(1,*,err=1000,end=1000) ispin,iorb,isflip,ipar,irad,iovert
      ibox = 1
      read(1,*,err=1000,end=1000) i_more_n,ecfmin,ecfmax,espmin,espmax
      write(2,*) 'i_more_n,ecfmin,ecfmax,espmin,espmax',
     1i_more_n,ecfmin,ecfmax,espmin,espmax
      nocc = nos
      gpp = 1.d0
      dgpp = 0.d0
      ngpp = 1
      gph = 1.d0
      dgph = 0.d0
      ngph = 1
      read(1,err=1000,end=1001,fmt='(a80)') spurious_proc
      if(index(spurious_proc,'DIP').ne.0) then
       read(1,*,err=1000,end=1000) GPP,DGPP,NGPP,GPH,DGPH,NGPH
      end if
      if(index(spurious_proc,'RANGE CENTROID').ne.0) then
       read(1,*,err=1000,end=1000) e_centr_min,e_centr_max
      end if
      go to 1001
1000  stop '>>> PROBLEM IN READER '
1001  write(2,1)
      write(2,fmt='(a80)') purpose
      if(irpa.ne.0)write(2,fmt='(a80)') purpose2
      write(2,1)
      write(2,*)
      return
      end subroutine reader



      subroutine forceparam
      implicit double precision(a-h,o-z)
      character skyrme*5
      common/b_forceparam/skyrme
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
cc    HERE coefff. a la Dobaceski 
      common/b_Doba/c0rho,c0rho_dens,c1rho,c1rho_dens,c0deltarho,
     + c1deltarho,c0tau,c1tau,c0s,c0s_dens,c1s_dens,c0nablaj,c1nablaj,
     + c0deltas,c1deltas,c0t,c1t


    2 format(10(1x,e12.5))
      write(2,*) 'Interaction used: '
      write(2,fmt='(a5)') skyrme
      write(2,*)
      if(index(skyrme,'oSIo').ne.0)then 
       t0=-1057.3d0
       t1=235.9d0
       t2=-100.d0
       t3=14463.5d0
       x0=0.56d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=120.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
      if(index(skyrme,'oSII').ne.0)then 
       t0=-1169.9d0
       t1=585.6d0
       t2=-27.1d0
       t3=9331.1d0
       x0=0.34d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=105.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
      if(index(skyrme,'SIII').ne.0)then 
       t0=-1128.75d0
       t1=395.3d0
       t2=-95.d0
       t3=14000.d0
       x0=0.45d0
       x1=0.d0
       x2=0.d0
       x3=1.d0
       w0_0=120.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
      if(index(skyrme,'oS3*').ne.0)then 
       t0=-1121.d0
       t1=400.d0
       t2=-533.33d0
       t3=14000.d0
       x0=0.43225d0
       x1=0.35d0
       x2=-0.9875d0
       x3=0.d0
       w0_0=120.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
      if(index(skyrme,'SV-').ne.0)then 
       t0=-1248.29d0
       t1=970.56d0
       t2=107.22d0
       t3=0.d0
       x0=-0.17d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=150.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
      if(index(skyrme,'SVII').ne.0)then 
       t0=-1096.8d0
       t1=246.2d0
       t2=-148.d0
       t3=17726.d0
       x0=0.62d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=112.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=1.d0
       go to 1001
      end if
C     Ska: NPA 258 (1976) 301
      if(index(skyrme,'Ska').ne.0)then
       t0=-1602.78d0
       t1=570.88d0
       t2=-67.7d0
       t3=8000.d0
       x0=-0.02d0
       x1=0.d0
       x2=0.d0
       x3=-0.286d0
       w0_0=125.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     SkM: NPA 336 (1980) 155
      if(index(skyrme,'SkM0').ne.0)then
       t0=-2645.d0
       t1=385.d0
       t2=-120.d0
       t3=15595.d0
       x0=0.09d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=130.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.166666666d0
       go to 1001
      end if
C     SkMP: PRC 40 (1989) 2834
      if(index(skyrme,'SkMP').ne.0)then
       t0=-2372.24d0
       t1=503.623d0
       t2=57.2783d0
       t3=12585.3d0
       x0=-0.157563d0
       x1=-0.402886d0
       x2=-2.95571d0
       x3=-0.267933d0
       w0_0=160.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.166666666d0
       go to 1001
      end if
C     SGI: PLB 106 (1981) 379
      if(index(skyrme,'SG1').ne.0)then 
       t0=-1603.d0
       t1=515.9d0
       t2=84.5d0
       t3=8000.d0
       x0=-0.02d0
       x1=-0.5d0
       x2=-1.731d0
       x3=0.1381d0
       w0_0=115.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     SGII: PLB 106 (1981) 379
      if(index(skyrme,'SGII').ne.0)then 
       t0=-2645.d0
       t1=340.d0
       t2=-41.9d0
       t3=15595.d0
       x0=0.09d0
       x1=-0.0588d0
       x2=1.425
       x3=0.06044
       w0_0=105.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     SKM*: NPA 386 (1982) 79
      if(index(skyrme,'SKM*').ne.0.or.index(skyrme,'SkM*').ne.0)then 
       t0=-2645.d0
       t1=410.d0
       t2=-135.d0
       t3=15595.d0
       x0=0.09d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=130.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     RATP: Astron. Astrophys. 116 (1982) 183
      if(index(skyrme,'RATP').ne.0)then 
       t0=-2160.d0
       t1=513.d0
       t2=121.d0
       t3=11600.d0
       x0=0.418d0
       x1=-0.360d0
       x2=-2.290d0
       x3=0.586d0
       w0_0=120.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.20d0
       go to 1001
      end if
C     SkT1: NPA 420 (1984) 297
      if(index(skyrme,'Tond1').ne.0)then 
       t0=-1794.d0
       t1=298.d0
       t2=-298.d0
       t3=12812.d0
       x0=0.154d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.089d0
       w0_0=110.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3333333d0
       go to 1001
      end if
C     SkT2: NPA 420 (1984) 297
      if(index(skyrme,'Tond2').ne.0)then 
       t0=-1791.6d0
       t1=300.d0
       t2=-300.d0
       t3=12792.d0
       x0=0.154d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.089d0
       w0_0=120.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3333333d0
       go to 1001
      end if
C     SkT3: NPA 420 (1984) 297
      if(index(skyrme,'Tond3').ne.0)then 
       t0=-1791.8d0
       t1=298.5d0
       t2=-99.5d0
       t3=12794.d0
       x0=0.138d0
       x1=-1.d0
       x2=1.d0
       x3=0.075d0
       w0_0=126.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3333333d0
       go to 1001
      end if
C     SkT4: NPA 420 (1984) 297
      if(index(skyrme,'Tond4').ne.0)then 
       t0=-1808.8d0
       t1=303.4d0
       t2=-303.4d0
       t3=12980.d0
       x0=-0.177d0
       x1=-0.5d0
       x2=-0.5d0
       x3=-0.5d0
       w0_0=113.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3333333d0
       go to 1001
      end if
C     SkT5: NPA 420 (1984) 297
      if(index(skyrme,'Tond5').ne.0)then 
       t0=-2917.1d0
       t1=328.2d0
       t2=-328.2d0
       t3=18584.d0
       x0=-0.295d0
       x1=-0.5d0
       x2=-0.5d0
       x3=-0.5d0
       w0_0=114.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.16666667d0
       go to 1001
      end if
C     SkT6: NPA 420 (1984) 297
      if(index(skyrme,'Tond6').ne.0)then 
       t0=-1794.2d0
       t1=294.d0
       t2=-294.d0
       t3=12817.d0
       x0=0.392d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.5d0
       w0_0=107.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3333333d0
       go to 1001
      end if
C     SkT7: NPA 420 (1984) 297
      if(index(skyrme,'Tond7').ne.0)then 
       t0=-1892.5d0
       t1=366.6d0
       t2=-21.d0
       t3=11983.d0
       x0=0.334d0
       x1=-0.359d0
       x2=6.9d0
       x3=0.366d0
       w0_0=109.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.285d0
       go to 1001
      end if
C     SkT8: NPA 420 (1984) 297
      if(index(skyrme,'Tond8').ne.0)then 
       t0=-1892.5d0
       t1=367.d0
       t2=-228.76d0
       t3=11983.d0
       x0=0.448d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.695d0
       w0_0=109.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.285d0
       go to 1001
      end if
C     SkT9: NPA 420 (1984) 297
      if(index(skyrme,'Tond9').ne.0)then 
       t0=-1891.4d0
       t1=377.4d0
       t2=-239.16d0
       t3=11982.d0
       x0=0.441d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.686d0
       w0_0=130.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.285d0
       go to 1001
      end if
C     SkP: NPA 422 (1984) 103
      if(index(skyrme,'SkP').ne.0.or.index(skyrme,'SKP').ne.0)then 
       t0=-2931.7d0
       t1=320.62d0
       t2=-337.41d0
       t3=18709.d0
       x0=0.29215d0
       x1=0.65318d0
       x2=-0.53732d0
       x3=0.18103d0
       w0_0=100.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     Sktry: NPA 422 (1984) 103
      if(index(skyrme,'Sktry').ne.0.or.index(skyrme,'SKtry').ne.0)then 
       t0=-2931.7d0
       t1=0.d0
       t2=0.d0
       t3=18709.d0
       x0=0.d0
       x1=0.d0
       x2=0.d0
       x3=0.d0
       w0_0=0.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     Esigma: PRC 33 (1986) 335
      if(index(skyrme,'SkEs').ne.0)then 
       t0=-1664.05d0
       t1=358.83d0
       t2=-137.22d0
       t3=10931.5d0
       x0=1.077d0
       x1=0.d0
       x2=0.d0
       x3=1.6918d0
       w0_0=120.14d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.35d0
       go to 1001
      end if
C     Zsigma: PRC 33 (1986) 335
      if(index(skyrme,'SkZs0').ne.0)then 
       t0=-1983.76d0
       t1=362.25d0
       t2=-104.27d0
       t3=11861.4d0
       x0=1.1717d0
       x1=0.d0
       x2=0.d0
       x3=1.762d0
       w0_0=123.69d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     Zsigma*: PRC 33 (1986) 335
      if(index(skyrme,'SkZs*').ne.0)then 
       t0=-1987.64d0
       t1=380.92d0
       t2=-109.88d0
       t3=11837.7d0
       x0=0.8897d0
       x1=0.d0
       x2=0.d0
       x3=1.278d0
       w0_0=126.13d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     Rsigma: PRC 33 (1986) 335
      if(index(skyrme,'SkRs').ne.0)then 
       t0=-1798.d0
       t1=335.97d0
       t2=-84.81d0
       t3=11083.9d0
       x0=-0.4036d0
       x1=0.d0
       x2=0.d0
       x3=-0.8705d0
       w0_0=121.59d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3d0
       go to 1001
      end if
C     Gsigma: PRC 33 (1986) 335
      if(index(skyrme,'SkGs').ne.0)then 
       t0=-1800.16d0
       t1=336.23d0
       t2=-85.74d0
       t3=11113.5d0
       x0=-0.4862d0
       x1=0.d0
       x2=0.d0
       x3=-1.0295d0
       w0_0=121.86d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3d0
       go to 1001
      end if
C     SkSC1: NPA 528 (1991) 1
c     if(index(skyrme,'SkSC1stop').ne.0)then 
c      t0=-1788.59d0
c      t1=282.623d0
c      t2=-282.623d0
c      t3=12775.3d0
c      x0=0.72d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=1.04564d0
c      w0_0=126.997d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC2: NPA 528 (1991) 1
c     if(index(skyrme,'SkSC2').ne.0)then 
c      t0=-1791.64d0
c      t1=290.601d0
c      t2=-290.601d0
c      t3=12800.6d0
c      x0=0.38d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=0.59276d0
c      w0_0=128.781d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC3: NPA 528 (1991) 1
c     if(index(skyrme,'SkSC3').ne.0)then 
c      t0=-1788.11d0
c      t1=289.901d0
c      t2=-96.6337d0
c      t3=12771.4d0
c      x0=0.65d0
c      x1=-1.d0
c      x2=1.d0
c      x3=0.960976d0
c      w0_0=143.452d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC4: NPA 549 (1992) 155
c     if(index(skyrme,'SkSC4').ne.0)then 
c      t0=-1789.42d0
c      t1=283.467d0
c      t2=-283.467d0
c      t3=12782.3d0
c      x0=0.79d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=1.13871d0
c      w0_0=124.877d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC5: PRC 50 (1994) 460
c     if(index(skyrme,'SkSC5').ne.0)then 
c      t0=-1788.17d0
c      t1=281.931d0
c      t2=-281.931d0
c      t3=12771.9d0
c      x0=0.98d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=1.138526d0
c      w0_0=126.219d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC6: PRC 50 (1994) 460
c     if(index(skyrme,'SkSC6').ne.0)then 
c      t0=-1792.47d0
c      t1=291.964d0
c      t2=-291.964d0
c      t3=12805.7d0
c      x0=0.370038d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=0.581085d0
c      w0_0=126.014d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkSC10: PRC 50 (1994) 460
c     if(index(skyrme,'SkSC10').ne.0)then 
c      t0=-1795.12d0
c      t1=298.950d0
c      t2=-298.950d0
c      t3=12827.7d0
c      x0=0.159124d0
c      x1=-0.5d0
c      x2=-0.5d0
c      x3=0.292918d0
c      w0_0=127.137d0
c      w2p_0=w0_0
c      t13=0.d0
c      alfe=0.333333d0
c      go to 1001
c     end if
C     SkI1: NPA 584 (1995) 467
      if(index(skyrme,'SkI1').ne.0.or.index(skyrme,'SKI1').ne.0)then 
       t0=-1913.62d0
       t1=439.809d0
       t2=2697.59d0
       t3=10592.3d0
       x0=-0.9545d0
       x1=-5.7824d0
       x2=-1.2874d0
       x3=-1.5614d0
       w0_0=124.260d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SkI2: NPA 584 (1995) 467
      if(index(skyrme,'SkI2').ne.0.or.index(skyrme,'SKI2').ne.0)then 
       t0=-1915.43d0
       t1=438.449d0
       t2=305.446d0
       t3=10548.9d0
       x0=-0.2108d0
       x1=-1.7376d0
       x2=-1.5336d0
       x3=-0.178d0
       w0_0=120.602d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SkI3: NPA 584 (1995) 467
      if(index(skyrme,'SkI3').ne.0.or.index(skyrme,'SKI3').ne.0)then
       t0=-1762.88d0
       t1=561.608d0
       t2=-227.09d0
       t3=8106.2d0
       x0=0.3083d0
       x1=-1.1722d0
       x2=-1.0907d0
       x3=1.2926d0
       w0_0=188.260d0
       w2p_0=0.d0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SkI4: NPA 584 (1995) 467
      if(index(skyrme,'SkI4').ne.0.or.index(skyrme,'SKI4').ne.0)then 
       t0=-1855.827d0
       t1=473.829d0
       t2=1006.855d0
       t3=9703.607d0
       x0=0.40508d0
       x1=-2.88915d0
       x2=-1.32515d0
       x3=1.1452d0
       w0_0=366.194d0
       w2p_0=-360.702
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SkI5: NPA 584 (1995) 467
      if(index(skyrme,'SkI5').ne.0.or.index(skyrme,'SKI5').ne.0)then 
       t0=-1772.91d0
       t1=550.840d0
       t2=-126.685d0
       t3=8206.25d0
       x0=-0.1171d0
       x1=-1.3088d0
       x2=-1.0487d0
       x3=0.341d0
       w0_0=123.63d0
       w2p_0=123.63d0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SLy230a: NPA 627 (1997) 710
      if(index(skyrme,'230a').ne.0)then 
       t0=-2490.23d0
       t1=489.53d0
       t2=-566.58d0
       t3=13803.d0
       x0=1.1318d0
       x1=-0.8426d0
       x2=-1.d0
       x3=1.9219d0
       w0_0=131.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     SLy4: NPA 643 (1998) 441
      if(index(skyrme,'SLy4').ne.0.or.index(skyrme,'Sly4').ne.0)then 
       t0=-2488.91d0
       t1=486.82d0
       t2=-546.390d0
       t3=13777.d0
       x0=0.834d0
       x1=-0.344d0
       x2=-1.d0
       x3=1.354d0
       w0_0=123.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     SLy5: NPA 643 (1998) 441
      if(index(skyrme,'SLy5').ne.0.or.index(skyrme,'Sly5').ne.0)then 
       t0=-2484.88d0
       t1=483.13d0
       t2=-549.40d0
       t3=13763.d0
       x0=0.778d0
       x1=-0.328d0
       x2=-1.d0
       x3=1.267d0
       w0_0=126.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     SkX: PRC 58 (1998) 220   
      if(index(skyrme,'SkX').ne.0.or.index(skyrme,'SKX').ne.0)then 
       t0=-1445.3d0
       t1=246.9d0
       t2=-131.8d0
       t3=12103.9d0
       x0=0.340d0
       x1=0.580d0
       x2=0.127d0
       x3=0.03d0
       w0_0=148.6d0
c_TEST
       w0_0=159.7d0
       w2p_0=0.0d0
       t13=0.d0
       alfe=0.5d0
       go to 1001
      end if
C     SkO: PRC 60 (1999) 014316
      if(index(skyrme,'SkO-').ne.0.or.index(skyrme,'SKO-').ne.0)then 
       t0=-2103.653d0
       t1=303.352d0
       t2=791.674d0
       t3=13553.252d0
       x0=-0.210701d0
       x1=-2.810752d0
       x2=-1.461595d0
       x3=-0.429881d0
       w0_0=353.156d0
       w2p_0=-397.498d0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     SkOprime: PRC 60 (1999) 014316
      if(index(skyrme,'SkOp').ne.0.or.index(skyrme,'SKOp').ne.0)then 
       t0=-2099.419d0
       t1=301.531d0
       t2=154.781d0
       t3=13526.464d0
       x0=-0.029503d0
       x1=-1.325732d0
       x2=-2.323439d0
       x3=-0.147404d0
       w0_0=287.79d0
       w2p_0=-165.7776d0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     MSk1: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk1').ne.0)then 
       t0=-1813.03d0
       t1=274.828d0
       t2=-274.828d0
       t3=13050.1d0
       x0=0.365395d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.449882d0
       w0_0=116.708d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk2: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk2').ne.0)then 
       t0=-1830.67d0
       t1=260.301d0
       t2=-293.742d0
       t3=13442.1d0
       x0=0.356875d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.409759d0
       w0_0=116.663d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk3: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk3').ne.0)then 
       t0=-1810.32d0
       t1=269.092d0
       t2=-269.092d0
       t3=13027.5d0
       x0=0.631485d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.903680d0
       w0_0=116.871d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk4: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk4').ne.0)then 
       t0=-1827.96d0
       t1=254.129d0
       t2=-287.569d0
       t3=13419.5d0
       x0=0.610360d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.835063d0
       w0_0=115.943d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk5: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk5').ne.0)then 
       t0=-1827.96d0
       t1=254.326d0
       t2=-287.766d0
       t3=13419.5d0
       x0=0.605152d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.827182d0
       w0_0=115.932d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk6: Phys. Rev. C 62 (2000) 024308
      if(index(skyrme,'MSk6').ne.0)then 
       t0=-1827.96d0
       t1=258.483d0
       t2=-291.924d0
       t3=13419.5d0
       x0=0.576591d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.783956d0
       w0_0=118.807d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     MSk7: At. Data Nucl. Data Tables 77 (2001) 311
      if(index(skyrme,'MSk7').ne.0)then 
       t0=-1828.23d0
       t1=259.4d0
       t2=-292.84d0
       t3=13421.7d0
       x0=0.576761d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.78529d0
       w0_0=118.807d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v110: PRC 64 (2001) 027301
      if(index(skyrme,'v110').ne.0)then 
       t0=-1827.96d0
       t1=252.949d0
       t2=-347.19d0
       t3=13419.5d0
       x0=0.627871d0
       x1=-0.5d0
       x2=-0.631341d0
       x3=0.894757d0
       w0_0=113.366d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v105: PRC 64 (2001) 027301
      if(index(skyrme,'v105').ne.0)then 
       t0=-1827.96d0
       t1=253.114d0
       t2=-286.554d0
       t3=13419.5d0
       x0=0.611956d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.837479d0
       w0_0=115.404d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v100: PRC 64 (2001) 027301
      if(index(skyrme,'v100').ne.0)then 
       t0=-1827.96d0
       t1=253.702d0
       t2=-220.261d0
       t3=13419.5d0
       x0=0.564938d0
       x1=-0.5d0
       x2=-0.272269d0
       x3=0.729808d0
       w0_0=117.008d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v090: PRC 64 (2001) 027301
      if(index(skyrme,'v090').ne.0)then 
       t0=-1827.96d0
       t1=253.565d0
       t2=-199.305d0
       t3=13419.5d0
       x0=0.559906d0
       x1=-0.1d0
       x2=-0.169977d0
       x3=0.637001d0
       w0_0=114.933d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v080: PRC 64 (2001) 027301
      if(index(skyrme,'v080').ne.0)then 
       t0=-1827.96d0
       t1=251.271d0
       t2=-168.233d0
       t3=13419.5d0
       x0=0.528552d0
       x1=0.4d0
       x2=0.019271d0
       x3=0.483059d0
       w0_0=109.448d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v075: PRC 64 (2001) 027301
      if(index(skyrme,'v075').ne.0)then 
       t0=-1828.64d0
       t1=252.349d0
       t2=-171.327d0
       t3=13425.1d0
       x0=0.521301d0
       x1=0.75d0
       x2=0.001069d0
       x3=0.408383d0
       w0_0=105.640d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     v070: PRC 64 (2001) 027301
      if(index(skyrme,'v070').ne.0)then 
       t0=-1828.64d0
       t1=248.527d0
       t2=-176.49d0
       t3=13425.1d0
       x0=0.517539d0
       x1=1.2d0
       x2=-0.051769d0
       x3=0.329698d0
       w0_0=98.0860d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     BSk1: NPA 700 (2002) 142
      if(index(skyrme,'BSk1-').ne.0)then 
       t0=-1830.45d0
       t1=262.970d0
       t2=-296.446d0
       t3=13444.7d0
       x0=0.599988d0
       x1=-0.5d0
       x2=-0.5d0
       x3=0.823074d0
       w0_0=117.97d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.33333333d0
       go to 1001
      end if
C     BSk14
      if(index(skyrme,'BSk14').ne.0)then 
       t0=-1822.67d0
       t1=377.47d0
       t2=-2.41056d0
       t3=11406.3d0
       x0=0.302096d0
       x1=-0.823575d0
       x2=61.9411d0
       x3=0.47346d0
       w0_0=135.565d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3d0
       go to 1001
      end if
C     SK255: PRC 68 (2003) 031304
      if(index(skyrme,'SK255').ne.0.or.index(skyrme,'Sk255').ne.0)then 
       t0=-1689.35d0
       t1=389.3d0
       t2=-126.07d0
       t3=10989.6d0
       x0=-0.1461d0
       x1=0.116d0
       x2=0.0012d0
       x3=-0.7449d0
       w0_0=95.39d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3563
       go to 1001
      end if
C     SK270: PRC 68 (2003) 031304
      if(index(skyrme,'SK272').ne.0.or.index(skyrme,'Sk272').ne.0)then 
       t0=-1496.84d0
       t1=397.66d0
       t2=-112.82d0
       t3=10191.64d0
       x0=0.0008d0
       x1=0.0102d0
       x2=0.002d0
       x3=-0.5519d0
       w0_0=106.58d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.4492
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy26
      if(index(skyrme,'asy26').ne.0)then 
       t0=-2453.3851d0
       t1=475.0118
       t2=-494.1886d0
       t3=13475.0482d0
       x0=2.1607d0
       x1=-0.4159d0
       x2=-1.d0
       x3=3.4314d0
       w0_0=123.4102d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy28
      if(index(skyrme,'asy28').ne.0)then 
       t0=-2461.9026d0
       t1=472.0822d0
       t2=-530.2394d0
       t3=13599.5253d0
       x0=1.6830
       x1=-0.3349d0
       x2=-1.d0
       x3=2.6602d0
       w0_0=124.7655d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy30
      if(index(skyrme,'asy30').ne.0)then 
       t0=-2477.0134d0
       t1=475.6039d0
       t2=-533.6572d0
       t3=13709.0945d0
       x0=1.1408d0
       x1=-0.3365d0
       x2=-1.d0
       x3=1.8302d0
       w0_0=117.4816d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy32
      if(index(skyrme,'asy32').ne.0)then 
       t0=-2491.8497d0
       t1=486.6144d0
       t2=-579.6609d0
       t3=13842.3513d0
       x0=0.9962d0
       x1=-0.2791d0
       x2=-1.d0
       x3=1.5763d0
       w0_0=128.0277
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy34
      if(index(skyrme,'asy34').ne.0)then 
       t0=-2503.4549
       t1=489.1119d0
       t2=-591.6478d0
       t3=13939.3951d0
       x0=0.7308d0
       x1=-0.2622d0
       x2=-1.d0
       x3=1.1527d0
       w0_0=126.8462d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy36
      if(index(skyrme,'asy36').ne.0)then 
       t0=-2513.9507d0
       t1=491.0175d0
       t2=-615.6619d0
       t3=14045.3088d0
       x0=0.5457d0
       x1=-0.2210d0
       x2=-1.d0
       x3=0.8437d0
       w0_0=128.5284d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy38
      if(index(skyrme,'asy38').ne.0)then 
       t0=-2524.5940d0
       t1=494.0606d0
       t2=-635.51d0
       t3=14142.7263d0
       x0=0.3550
       x1=-0.1916d0
       x2=-1.d0
       x3=0.5315
       w0_0=129.1997d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: asy40
      if(index(skyrme,'asy40').ne.0)then 
       t0=-2531.6883d0
       t1=489.5977d0
       t2=-625.9761d0
       t3=14203.4838d0
       x0=0.1249d0
       x1=-0.1940d0
       x2=-1.d0
       x3=0.1662d0
       w0_0=122.0453d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=26
      if(index(skyrme,'s4026').ne.0)then 
       t0=-2295.8967d0
       t1=526.0481d0
       t2=-343.6887d0
       t3=11749.5781d0
       x0=2.0777d0
       x1=-0.8554d0
       x2=-1.d0
       x3=3.7114d0
       w0_0=117.9997d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=28
      if(index(skyrme,'s4028').ne.0)then 
       t0=-2296.5337d0
       t1=515.5864d0
       t2=-336.4641d0
       t3=11786.3268d0
       x0=1.7714d0
       x1=-0.8491d0
       x2=-1.d0
       x3=3.1703d0
       w0_0=118.8025d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=30
      if(index(skyrme,'s4030').ne.0)then 
       t0=-2317.5462d0
       t1=532.8981d0
       t2=-355.0194d0
       t3=11907.0695d0
       x0=1.2924d0
       x1=-0.8449d0
       x2=-1.d0
       x3=2.3824d0
       w0_0=118.2625d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=32
      if(index(skyrme,'s4032').ne.0)then 
       t0=-2321.0874d0
       t1=520.2351d0
       t2=-491.5039d0
       t3=12164.329d0
       x0=1.2626d0
       x1=-0.5755d0
       x2=-1.d0
       x3=2.1916d0
       w0_0=136.9896d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=34
      if(index(skyrme,'s4034').ne.0)then 
       t0=-2329.7936d0
       t1=510.9318d0
       t2=-423.7083d0
       t3=12179.5398d0
       x0=0.8546d0
       x1=-0.6759d0
       x2=-1.d0
       x3=1.5447d0
       w0_0=117.9319d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=36
      if(index(skyrme,'s4036').ne.0)then 
       t0=-2337.5162d0
       t1=507.8342d0
       t2=-424.5683d0
       t3=12253.2247d0
       x0=0.7156d0
       x1=-0.6664d0
       x2=-1.d0
       x3=1.2882d0
       w0_0=118.4138d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=38
      if(index(skyrme,'s4038').ne.0)then 
       t0=-2355.7652d0
       t1=523.0509d0
       t2=-465.1452d0
       t3=12389.163d0
       x0=0.379d0
       x1=-0.6244d0
       x2=-1.d0
       x3=0.7299d0
       w0_0=118.2375d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=240, J=40
      if(index(skyrme,'s4040').ne.0)then 
       t0=-2359.0006d0
       t1=513.3245d0
       t2=-475.5631d0
       t3=12466.4566d0
       x0=0.2959d0
       x1=-0.5827d0
       x2=-1.d0
       x3=0.5558d0
       w0_0=119.4059d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=28
      if(index(skyrme,'s5028').ne.0)then 
       t0=-2157.1788d0
       t1=599.6077d0
       t2=-396.2711d0
       t3=10352.3783d0
       x0=1.4017d0
       x1=-0.9193d0
       x2=-1.d0
       x3=2.9016d0
       w0_0=138.8744d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=30
      if(index(skyrme,'s5030').ne.0)then 
       t0=-2155.1989d0
       t1=570.3777d0
       t2=-299.172d0
       t3=10324.2868d0
       x0=1.6366d0
       x1=-1.0334d0
       x2=-1.d0
       x3=3.2711d0
       w0_0=131.2176d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=32
      if(index(skyrme,'s5032').ne.0)then 
       t0=-2172.5098d0
       t1=583.5116d0
       t2=-348.9922d0
       t3=10472.9193d0
       x0=1.4649d0
       x1=-0.9688d0
       x2=-1.d0
       x3=2.914d0
       w0_0=141.5897d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=34
      if(index(skyrme,'s5034').ne.0)then 
       t0=-2179.6101d0
       t1=584.0692d0
       t2=-331.05254d0
       t3=10501.7935d0
       x0=1.2593d0
       x1=-0.9977d0
       x2=-1.d0
       x3=2.5371d0
       w0_0=141.0935d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=36
      if(index(skyrme,'s5036').ne.0)then 
       t0=-2195.4023d0
       t1=599.4421
       t2=-365.5564d0
       t3=10609.2498d0
       x0=0.9683d0
       x1=-0.9629d0
       x2=-1.d0
       x3=2.002d0
       w0_0=145.6384d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=38
      if(index(skyrme,'s5038').ne.0)then 
       t0=-2205.3575d0
       t1=602.0931d0
       t2=-375.9978d0
       t3=10690.1918d0
       x0=0.7284d0
       x1=-0.9493d0
       x2=-1.d0
       x3=1.5517d0
       w0_0=145.124d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=40
      if(index(skyrme,'s5040').ne.0)then 
       t0=-2214.2498d0
       t1=600.7144d0
       t2=-404.3826d0
       t3=10802.958d0
       x0=0.5278d0
       x1=-0.9019d0
       x2=-1.d0
       x3=1.1569d0
       w0_0=146.3581d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=26
      if(index(skyrme,'s2026').ne.0)then 
       t0=-2656.2503d0
       t1=472.3188d0
       t2=-601.7814d0
       t3=15329.0650d0
       x0=1.6352d0
       x1=-0.1424d0
       x2=-1.d0
       x3=2.4638d0
       w0_0=123.8088
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=28
      if(index(skyrme,'s2028').ne.0)then 
       t0=-2634.0562d0
       t1=418.0827d0
       t2=-603.0149d0
       t3=15353.9648d0
       x0=1.1019d0
       x1=0.0686d0
       x2=-1.d0
       x3=1.6104d0
       w0_0=96.024d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=30
      if(index(skyrme,'s2030').ne.0)then 
       t0=-2647.2077d0
       t1=435.8828d0
       t2=-643.9488d0
       t3=15442.7977d0
       x0=1.1117d0
       x1=0.0827d0
       x2=-1.d0
       x3=1.599d0
       w0_0=118.9192d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=32
      if(index(skyrme,'s2032').ne.0)then 
       t0=-2675.2552d0
       t1=462.6729d0
       t2=-715.2266d0
       t3=15654.2779d0
       x0=0.844d0
       x1=0.1228d0
       x2=-1.d0
       x3=1.1937d0
       w0_0=130.5423d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=34
      if(index(skyrme,'s2034').ne.0)then 
       t0=-2682.7194d0
       t1=454.1955d0
       t2=-721.5723d0
       t3=15752.5009d0
       x0=0.8123d0
       x1=0.1726d0
       x2=-1.d0
       x3=1.1097d0
       w0_0=132.2429d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=36
      if(index(skyrme,'s2036').ne.0)then 
       t0=-2693.3779d0
       t1=468.2162d0
       t2=-762.2838d0
       t3=15835.528d0
       x0=0.5935d0
       x1=0.1966d0
       x2=-1.d0
       x3=0.7769d0
       w0_0=142.4865d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=38
      if(index(skyrme,'s2038').ne.0)then 
       t0=-2704.1367d0
       t1=474.7342d0
       t2=-783.3956d0
       t3=15922.9639d0
       x0=0.4616d0
       x1=0.2126d0
       x2=-1.d0
       x3=0.5645d0
       w0_0=149.1597d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=220, J=40
      if(index(skyrme,'s2040').ne.0)then 
       t0=-2723.4256d0
       t1=483.567d0
       t2=-805.0598d0
       t3=16070.5815d0
       x0=0.3939d0
       x1=0.2211d0
       x2=-1.d0
       x3=0.4419d0
       w0_0=155.2334d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=28, a la Shlomo
      if(index(skyrme,'h5028').ne.0)then 
       t0=-1680.0032d0
       t1=420.2722d0
       t2=-609.118d0
       t3=10744.453d0
       x0=0.7882d0
       x1=0.0954d0
       x2=-1.d0
       x3=1.5541d0
       w0_0=106.2412d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=30, a la Shlomo
      if(index(skyrme,'h5030').ne.0)then 
       t0=-1681.3213d0
       t1=414.1716d0
       t2=-590.3460d0
       t3=10755.4218d0
       x0=0.8088d0
       x1=0.0811d0
       x2=-1.d0
       x3=1.5345d0
       w0_0=112.4526d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=32, a la Shlomo
      if(index(skyrme,'h5032').ne.0)then 
       t0=-1691.5635d0
       t1=423.1798d0
       t2=-613.1853d0
       t3=10862.7212d0
       x0=0.5807d0
       x1=0.0952d0
       x2=-1.d0
       x3=1.0546d0
       w0_0=115.0102d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=34, a la Shlomo
      if(index(skyrme,'h5034').ne.0)then 
       t0=-1702.0903d0
       t1=433.8932d0
       t2=-640.6202d0
       t3=10972.5981d0
       x0=0.3772d0
       x1=0.1112d0
       x2=-1.d0
       x3=0.621d0
       w0_0=117.7697d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=36, a la Shlomo
      if(index(skyrme,'h5036').ne.0)then 
       t0=-1705.9318d0
       t1=430.6711d0
       t2=-629.3056d0
       t3=11012.9164d0
       x0=0.3041d0
       x1=0.1024d0
       x2=-1.d0
       x3=0.4289d0
       w0_0=120.1445d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=38, a la Shlomo
      if(index(skyrme,'h5038').ne.0)then 
       t0=-1719.751d0
       t1=449.2871d0
       t2=-678.7636d0
       t3=11157.4426d0
       x0=0.0864d0
       x1=0.1311d0
       x2=-1.d0
       x3=-0.03d0
       w0_0=122.7334d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=40, a la Shlomo
      if(index(skyrme,'h5040').ne.0)then 
       t0=-1724.0371d0
       t1=449.6155d0
       t2=-678.0642d0
       t3=11202.9305d0
       x0=-0.0128d0
       x1=0.1297d0
       x2=-1.d0
       x3=-0.2694d0
       w0_0=124.7361d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=250, J=42, a la Shlomo
      if(index(skyrme,'h5042').ne.0)then 
       t0=-1732.762d0
       t1=458.7049d0
       t2=-701.2595d0
       t3=11294.5303d0
       x0=-0.1316d0
       x1=0.1416d0
       x2=-1.d0
       x3=-0.545d0
       w0_0=127.7638d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=28, a la Shlomo
      if(index(skyrme,'h6028').ne.0)then 
       t0=-1588.6915d0
       t1=490.9261d0
       t2=-433.9863d0
       t3=8873.5405d0
       x0=1.1467d0
       x1=-0.5786d0
       x2=-1.d0
       x3=2.8843d0
       w0_0=124.6211d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=30, a la Shlomo
      if(index(skyrme,'h6030').ne.0)then 
       t0=-1597.9895d0
       t1=498.5904d0
       t2=-480.1063d0
       t3=9019.3198d0
       x0=0.8924d0
       x1=-0.5069d0
       x2=-1.d0
       x3=2.2429d0
       w0_0=130.1215d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=32, a la Shlomo
      if(index(skyrme,'h6032').ne.0)then 
       t0=-1598.7848d0
       t1=480.9093d0
       t2=-519.913d0
       t3=9202.4395d0
       x0=0.7801d0
       x1=-0.3767d0
       x2=-1.d0
       x3=1.8388d0
       w0_0=130.2798d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=34, a la Shlomo
      if(index(skyrme,'h6034').ne.0)then 
       t0=-1604.0281d0
       t1=478.2814d0
       t2=-552.3588d0
       t3=9335.8041d0
       x0=0.6372d0
       x1=-0.3016d0
       x2=-1.d0
       x3=1.4406d0
       w0_0=132.8337d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=36, a la Shlomo
      if(index(skyrme,'h6036').ne.0)then 
       t0=-1611.0753d0
       t1=480.6789d0
       t2=-601.0021d0
       t3=9491.0354d0
       x0=0.4744d0
       x1=-0.21d0
       x2=-1.d0
       x3=1.01d0
       w0_0=136.0071d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=38, a la Shlomo
      if(index(skyrme,'h6038').ne.0)then 
       t0=-1617.6071d0
       t1=478.2272d0
       t2=-643.2399d0
       t3=9655.8553d0
       x0=0.3608d0
       x1=-0.1146d0
       x2=-1.d0
       x3=0.6826d0
       w0_0=137.0264d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=40, a la Shlomo
      if(index(skyrme,'h6040').ne.0)then 
       t0=-1621.5358d0
       t1=477.9287d0
       t2=-686.5383d0
       t3=9781.9373d0
       x0=0.2777d0
       x1=-0.0253d0
       x2=-1.d0
       x3=0.4285d0
       w0_0=142.6583d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=260, J=42, a la Shlomo
      if(index(skyrme,'h6042').ne.0)then 
       t0=-1626.1089d0
       t1=478.9898d0
       t2=-710.7283d0
       t3=9781.522d0
       x0=0.1639d0
       x1=0.0205d0
       x2=-1.d0
       x3=0.137d0
       w0_0=144.8878d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=28, a la Shlomo
      if(index(skyrme,'h7028').ne.0)then 
       t0=-1514.496d0
       t1=596.1779d0
       t2=-222.6384d0
       t3=6936.6843d0
       x0=1.3677d0
       x1=-1.1823d0
       x2=-1.d0
       x3=4.6994d0
       w0_0=141.2427d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=30, a la Shlomo
      if(index(skyrme,'h7030').ne.0)then 
       t0=-1521.0886d0
       t1=595.9587d0
       t2=-245.5034d0
       t3=7053.4006d0
       x0=1.07d0
       x1=-1.1435d0
       x2=-1.d0
       x3=3.7811d0
       w0_0=141.4883d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=30, a la Shlomo
      if(index(skyrme,'h7030').ne.0)then 
       t0=-1521.0886d0
       t1=595.9587d0
       t2=-245.5034d0
       t3=7053.4006d0
       x0=1.07d0
       x1=-1.1435d0
       x2=-1.d0
       x3=3.7811d0
       w0_0=141.4883d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=32, a la Shlomo
      if(index(skyrme,'h7032').ne.0)then 
       t0=-1527.799d0
       t1=594.7912d0
       t2=-272.2628d0
       t3=7184.2076d0
       x0=0.8373d0
       x1=-1.0971d0
       x2=-1.d0
       x3=3.0347d0
       w0_0=141.5227d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=34, a la Shlomo
      if(index(skyrme,'h7034').ne.0)then 
       t0=-1534.1662d0
       t1=594.0466d0
       t2=-307.964d0
       t3=7325.8157d0
       x0=0.6273d0
       x1=-1.0367d0
       x2=-1.d0
       x3=2.3551d0
       w0_0=142.1694d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=36, a la Shlomo
      if(index(skyrme,'h7036').ne.0)then 
       t0=-1539.5661d0
       t1=590.9636d0
       t2=-330.2249d0
       t3=7444.6327d0
       x0=0.4648d0
       x1=-0.9946d0
       x2=-1.d0
       x3=1.8227d0
       w0_0=141.7349d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=38, a la Shlomo
      if(index(skyrme,'h7038').ne.0)then 
       t0=-1543.3239d0
       t1=584.6841d0
       t2=-361.6233d0
       t3=7581.0304d0
       x0=0.3126d0
       x1=-0.9323d0
       x2=-1.d0
       x3=1.3121d0
       w0_0=141.6285d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=40, a la Shlomo
      if(index(skyrme,'h7040').ne.0)then 
       t0=-1549.1667d0
       t1=587.3336d0
       t2=-369.0556d0
       t3=7643.8837d0
       x0=0.1749d0
       x1=-0.923d0
       x2=-1.d0
       x3=0.8989d0
       w0_0=141.7513d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces with different J introduced by Bennaceur-Meyer: K=270, J=42, a la Shlomo
      if(index(skyrme,'h7042').ne.0)then 
       t0=-1554.5684d0
       t1=591.4131d0
       t2=-391.2705d0
       t3=7721.8707d0
       x0=0.0058d0
       x1=-0.8922d0
       x2=-1.d0
       x3=0.4072d0
       w0_0=142.0649d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.356303d0
       go to 1001
      end if
C     forces by S. Goriely: G01 (gamma=0.25 ---> K=231 MeV)    
      if(index(skyrme,'sgo01').ne.0)then 
       t0=-2040.118468d0
       t1=406.678743d0
       t2=-166.877155d0
       t3=12468.2855d0
       x0=0.488011d0
       x1=-0.962986d0
       x2=-0.204160d0
       x3=0.847586d0
       w0_0=150.516407d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.25d0
       go to 1001
      end if
C     forces by S. Goriely: G02 (gamma=1/3 ---> K=245 MeV)    
      if(index(skyrme,'sgo02').ne.0)then 
       t0=-1714.385307d0
       t1=396.496212d0
       t2=-185.180996d0
       t3=10919.1925d0
       x0=0.446199d0
       x1=-0.98549d0
       x2=-0.348774d0
       x3=0.917413d0
       w0_0=147.648992d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.333333d0
       go to 1001
      end if
C     forces by S. Goriely: G06 (gamma=0.4 ---> K=258 MeV)                     
      if(index(skyrme,'sgo06').ne.0)then 
       t0=-1551.125020d0
       t1=387.838064d0
       t2=-172.717289d0
       t3=10295.642258d0
       x0=0.396169d0
       x1=-1.0194140
       x2=-0.321336d0
       x3=0.923664d0
       w0_0=146.072549d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.4d0
       go to 1001
      end if
C     forces by S. Goriely: G07 (gamma=0.5 ---> K=276 MeV)      
      if(index(skyrme,'sgo07').ne.0)then 
       t0=-1388.217876d0
       t1=375.327679d0
       t2=-174.347874d0
       t3=9915.158007d0
       x0=0.343815d0
       x1=-1.017146d0
       x2=-0.38383838d0
       x3=0.949256d0
       w0_0=141.630025d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.5d0
       go to 1001
      end if
C     forces by S. Goriely: G05 (gamma=0.6 ---> K=294 MeV)      
      if(index(skyrme,'sgo05').ne.0)then 
       t0=-1279.768167d0
       t1=363.840536d0
       t2=-190.990332d0
       t3=9946.407055d0
       x0=0.306822d0
       x1=-1.010894d0
       x2=-0.504422d0
       x3=1.004625d0
       w0_0=136.893428d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.6d0
       go to 1001
      end if
C     forces fitted on BHF: LNS (arXiv:nucl-th/0512086)      
      if(index(skyrme,'LNS0').ne.0)then 
       t0=-2484.97d0
       t1=266.735d0
       t2=-337.135d0
       t3=14588.2d0
       x0=0.06277d0
       x1=0.65845d0
       x2=-0.95382d0
       x3=-0.03413d0
       w0_0=96.0d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     Forces by Lesinski: T44
      if(index(skyrme,'T44').ne.0)then
       t0=-2485.67d0
       t1=494.477d0
       t2=-337.961d0
       t3=13794.75d0
       x0=0.721557d0
       x1=-0.661848d0
       x2=-0.803184d0
       x3=1.175908d0
       w0_0=161.367d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     Forces by Chamel: BSk17
      if(index(skyrme,'BSk17').ne.0)then
       t0=-1837.33d0
       t1=389.102d0
       t2=-3.1742d0
       t3=11523.8d0
       x0=0.411377d0
       x1=-0.832102d0
       x2=49.4875d0
       x3=0.654962d0
       w0_0=145.885d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.3d0
       go to 1001
      end if
C     New forces fitted on BHF by Lulu Li: set IV
      if(index(skyrme,'Li4').ne.0)then 
       t0=-2652.3194d0
       t1=182.9656d0
       t2=156.3080d0
       t3=15778.1701d0
       x0=0.4705d0
       x1=-4.3052d0
       x2=-1.3445d0
       x3=0.8464d0
       w0_0=90.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     New forces fitted on BHF by Lulu Li: set V
      if(index(skyrme,'Li5').ne.0)then 
       t0=-2878.7337d0
       t1=63.1982d0
       t2=152.5453d0
       t3=18094.9337d0
       x0=0.5679d0
       x1=-12.9810d0
       x2=-1.3519d0
       x3=0.9091d0
       w0_0=90.d0
       w2p_0=w0_0
       t13=0.d0
       alfe=0.1666666d0
       go to 1001
      end if
C     New forces fitted on BHF by Danilo G.
      if(index(skyrme,'LNSp1').ne.0)then
       t0=-2215.34455d0
       t1=532.53032d0
       t2=67.76211d0
       t3=10931.60885d0
       x0=0.46291d0
       x1=0.12809d0
       x2=-2.17349d0
       x3=0.61493d0
       alfe=0.16666667d0
       w0_0=116.78885d0
       w2p_0=w0_0
       t13=0.d0
       go to 1001
      end if
C     New forces fitted on BHF by Danilo G.
      if(index(skyrme,'LNSp2').ne.0)then
       t0=-2204.56407d0
       t1=506.97225d0
       t2=104.10407d0
       t3=10854.20582d0
       x0=0.31765d0
       x1=0.16399d0
       x2=-1.64450d0
       x3=0.38179d0
       alfe=0.16666667d0
       w0_0=111.65529d0
       w2p_0=w0_0
       t13=0.d0
       go to 1001
      end if
C     New forces fitted on BHF by Danilo G.
      if(index(skyrme,'LNSp5').ne.0)then
       t0=-2194.77853d0
       t1=482.51828d0
       t2=138.13769d0
       t3=10784.16615d0
       x0=0.13451d0
       x1=-0.09735d0
       x2=-1.3992d0
       x3=0.17601d0
       alfe=0.16666667d0
       w0_0=105.67416d0
       w2p_0=w0_0
       t13=0.d0
       go to 1001
      end if
C     SAMi0
      if(index(skyrme,'SAM0').ne.0)then
       t0=-1877.746d0
       t1=475.5856d0
       t2=-85.20021d0
       t3=10219.58d0
       x0=0.3197176d0
       x1=-0.5319419d0
       x2=-0.0137857d0
       x3=0.6883226d0
       alfe=0.2561388d0
       w0_0=137.0603d0
       w2p_0=42.32571d0
       t13=0.d0
       go to 1001
      end if
C     SAMi27
      if(index(skyrme,'SAM27').ne.0)then
       t0=-1876.09978010d0
       t1=481.08876425d0
       t2=-75.70799313d0
       t3=10184.61537483d0
       x0=0.48223611d0
       x1=-0.55796649d0
       x2=0.21305272d0
       x3=1.00218922d0
       alfe=0.25463412d0
       w0_0=81.96951358d0
       w2p_0=180.29552278d0
       t13=0.d0
       go to 1001
      end if
C     SAMi28
      if(index(skyrme,'SAM28').ne.0)then
       t0=-1872.51042551d0
       t1=475.54644843d0
       t2=-83.72657601d0
       t3=10181.87247955d0
       x0=0.51173115d0
       x1=-0.53519351d0
       x2=0.01268932d0
       x3=1.03255750d0
       alfe=0.25655506d0
       w0_0=94.53763726d0
       w2p_0=155.00180995d0
       t13=0.d0
       go to 1001
      end if
C     SAMi
      if(index(skyrme,'SAM29').ne.0)then
       t0=-1862.57029141d0
       t1=471.27726753d0
       t2=-92.75726511d0
       t3=10161.12750744d0
       x0=0.23677003d0
       x1=-0.510222736d0
       x2=-0.16371424d0
       x3=0.52385309d0
       alfe=0.26095501d0
       w0_0=150.86589854d0
       w2p_0=20.17035587d0
       t13=0.d0
       go to 1001
      end if
C     SAMi
      if(index(skyrme,'SAM30').ne.0)then
       t0=-1853.89273683d0
       t1=466.11794800d0
       t2=-101.44699522d0
       t3=10139.48429553d0
       x0=0.10253334d0
       x1=-0.48490834d0
       x2=-0.30800570d0
       x3=0.26024746d0
       alfe=0.26456590d0
       w0_0=185.72231632d0
       w2p_0=-60.91326849d0
       t13=0.d0
       go to 1001
      end if
C     SAMi
      if(index(skyrme,'SAM31').ne.0)then
       t0=-1844.26557358d0
       t1=460.72266259d0
       t2=-110.19682595d0
       t3=10112.27153205d0
       x0=-0.02357191d0
       x1=-0.45861146d0
       x2=-0.43123417d0
       x3=0.00789889d0
       alfe=0.26837131d0
       w0_0=216.85878758d0
       w2p_0=-133.50977370d0
       t13=0.d0
       go to 1001
      end if
C     fake
      if(index(skyrme,'fake').ne.0)then
       t0=-2400.26557358d0
       t1=460.72266259d0
       t2=-110.19682595d0
       t3=10000.27153205d0
       x0=-0.02357191d0
       x1=-0.45861146d0
       x2=-0.43123417d0
       x3=0.00789889d0
       alfe=0.26837131d0
       w0_0=216.85878758d0
       w2p_0=-133.50977370d0
       t13=0.d0
       go to 1001
      end if

      write(2,'(a5)') skyrme
      stop '>>> Skyrme parameters not found'
1001  write(2,2) t0,t1,t2,t3,x0,x1,x2,x3
      write(2,2) t13,alfe,w0_0

c     HERE   Write the coupling constants 'a la Dobaczewski'
      write(2,*)
      write(2,*) ' Coupling constants'
      write(2,*)

      c0rho = 3.d0/8.d0 * t0 
      c0rho_dens = t3/16.d0 
      write(2,*) 'C_0^\rho = ',c0rho,' + ',c0rho_dens,' *\rho^',alfe
      c1rho = (-t0-2.d0*t0*x0)/8.d0 
      c1rho_dens = (-t3-2.d0*t3*x3)/48.d0
      write(2,*) 'C_1^\rho = ',c1rho,' + ',c1rho_dens,' *\rho^',alfe

      c0deltarho = (-9.d0*t1+5.d0*t2+4.d0*t2*x2)/64.d0 
      write(2,*) 'C_0^\delta\rho = ',c0deltarho
      c1deltarho = (3.d0*t1+6.d0*t1*x1+t2+2.d0*t2*x2)/64.d0 
      write(2,*) 'C_1^\delta\rho = ',c1deltarho

      c0tau = (12.d0*t1+20.d0*t2+16.d0*t2*x2)/64.d0 
      write(2,*) 'C_0^\tau = ',c0tau
      c1tau = (-4.d0*t1-8.d0*t1*x1+4.d0*t2+8.d0*t2*x2)/64.d0 
      write(2,*) 'C_1^\tau = ',c1tau

      c0s = (-t0+2.d0*t0*x0)/8.d0 
      c0s_dens = (-t3+2.d0*t3*x3)/48.d0
      write(2,*) 'C_0^s = ',c0s,' + ',c0s_dens,' *\rho^',alfe
      c1s = -t0/8.d0 
      c1s_dens = -t3/48.d0
      write(2,*) 'C_1^s = ',c1s,' + ',c1s_dens,' *\rho^',alfe

      c0nablaj = -w0_0/2.d0-w2p_0/4.d0
      write(2,*) 'C_0^\nabla J = ',c0nablaj
      c1nablaj = -w2p_0/4.d0
      write(2,*) 'C_1^\nabla J = ',c1nablaj

      c0deltas = (3.d0*t1-6.d0*t1*x1+t2+2.d0*t2*x2)/64.d0 
      write(2,*) 'C_0^\Delta s = ',c0deltas
      c1deltas = (3.d0*t1+t2)/64.d0 
      write(2,*) 'C_1^\Delta s = ',c1deltas

      c0t = (-4.d0*t1+8.d0*t1*x1+4.d0*t2+8.d0*t2*x2)/64.d0 
      write(2,*) 'C_0^T = ',c0t
      c1t = (-4.d0*t1+4.d0*t2)/64.d0 
      write(2,*) 'C_1^T = ',c1t

      return
      end 



      subroutine create_list_unocc
      implicit double precision(a-h,o-z)
      include 'param.qrpa'
      parameter (lmx=11,neoc=1000)
      common/b_imore_n/i_more_n
      common/hf/nmax,nocc,nunocc,norb
      common/Basis/lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &             iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      COMMON/CTR2/IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      common/b_occup/nos,npo,ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
     &ioc04,ideg4,ioc05,ideg5,ioc06,ideg6
      common/bfo/fo(neoc) 
      common/b_gr0/igr0_j 
      dimension ichosen(0:lmx,2,0:1),lastn(0:lmx,2,0:1)

c     nocc = nos
      ichex=0
	if(igr0_j.eq.1)ichex=1

	write(2,*)
	write(2,*) 'ichex = ',ichex
	write(2,*)

      open(unit=99,status='old',file='lines')

      do i=1,1000
       read(99,*,end=4) nn(i),ll(i),lj(i),iq(i)
       if(iq(i).gt.0)ipr=i
       if((lj(i)+2*ispin).gt.(2*lmx+1))then
        write(2,*) '>>> Increase LMX in CREATE_LIST_UNOCC'
        stop
       end if
      end do

4     nocc=i-1
      nos=nocc
      npo=ipr

C_TEST
      do l_un=0,lmx
c     do l_un=0,1 
       do is=1,2
        do iqq=0,1
        ichosen(l_un,is,iqq)=0
        lastn(l_un,is,iqq)=0
        end do
       end do
      end do
      do 1 i=1,nocc
       jpmin = iabs(lj(i)-2*ispin)
       jpmax = lj(i)+2*ispin
       ipar_occ = 1-2*mod(ll(i),2)
       if(lj(i).lt.(2*ll(i)))iss1=1
       if(lj(i).gt.(2*ll(i)))iss1=2
       if(nn(i).gt.lastn(ll(i),iss1,iq(i)))lastn(ll(i),iss1,iq(i))=nn(i)
C_TEST
       do 2 l_un=0,lmx        
c      do 2 l_un=0,1
        ipar_unocc = 1-2*mod(l_un,2)
        itest = ipar_occ*ipar_unocc*ipar
c       if(itest.le.0)go to 2
        do 3 is=1,2
         if(l_un.eq.0.and.is.eq.1)go to 3
         iss = -1 + 2*(is-1)
         j_un = 2*l_un + iss
c        if(j_un.lt.jpmin)go to 3
c        if(j_un.gt.jpmax)go to 3
         if(ichex.eq.0)ichosen(l_un,is,iq(i))=1
c        if(ichex.eq.1.and.iq(i).eq.1)ichosen(l_un,is,0)=1
c        if(ichex.eq.1.and.iq(i).eq.0)ichosen(l_un,is,1)=1
    3   continue
    2  continue
    1 continue
      do iqq=1,0,-1
       do l_un=0,lmx
        do is=1,2
         if(ichosen(l_un,is,iqq).eq.1)then
          jjj = 2*l_un - 1 + 2*(is-1)
          write(99,*) 
     &    l_un,jjj,iqq,lastn(l_un,is,iqq)+1,lastn(l_un,is,iqq)+i_more_n
c         write(99,*) is,iqq
c         write(99,*) 'Next'
         end if
        end do
       end do
      end do

      close(99)

c     stop
      return
      end
