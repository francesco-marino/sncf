      subroutine cbcs_new
      implicit double precision(a-h,o-z)
      double precision ma
      parameter (nc=100)
      dimension a(0:4),b(0:4),co(-nc:nc),ec(-nc:nc),r2(-nc:nc)
      common/b_chf_1/irpa_chf,n_constr_step
      common/b_chf_2/constr_step
      common/b_chf_3/constr
      common/b_chf_4/etot,r2lambda
      common/b_chf_5/am3
      common/b_chf_6/e_di_tot,epote
      common/bma/ma
      common/bunits/pi,hbc,amc2,hbdm
      data a(0),a(1),a(2),a(3),a(4) /2.d0,-16.d0,0.d0,16.d0,-2.d0/
      data b(0),b(1),b(2),b(3),b(4) /-1.d0,16.d0,-30.d0,16.d0,-1.d0/
    5 format(13x,4(f12.6))

      fpi=4.d0*pi
      if(n_constr_step.gt.nc)then
       write(2,*) '>>> INCREASE DIM. OF EC AND R2'
       stop
      end if
      relinf=0.d0
      relsup=0.5d0
      r_equi2 = 0.d0
      r2lambda0 = 0.d0
      etot0 = 0.d0
      icount = 0
      sigma = 0.001d0

      write(2,*) 'Number in steps in lambda: ',n_constr_step

      do 1 i=-n_constr_step,n_constr_step,1
 
      constr = float(i)*constr_step

c     constr=0.08

      call bcsgen

      factor = 1.d0
c     factor = 1.d0 - 1.d0/ma
c     write(*,*) ma
c     write(*,*) factor

      etot_all = etot
      etot = etot-constr*factor*r2lambda/2.d0
      e_di_tot = e_di_tot / ma
      if(i.eq.-n_constr_step)then
       write(99,*) '  i  lambda        r2            E/A
     &           E/A(integral) diff%        EPOTE'
      end if
      if(i.eq.0)then
       r2lambda0=factor*r2lambda 
       etot0=etot
C_NEW
       etot0=e_di_tot
       r2rel=0.d0
      end if
      if(i.le.0)go to 7
      r2rel=(r2lambda0-r2lambda)/r2lambda0
      write(2,*) 'Using constr.: ',constr,'   ',r2lambda,r2rel
      co(i) = constr
      r2(i)= factor*r2lambda
      ec(i)= etot
C_NEW
      ec(i)= e_di_tot
    7 write(99,6) i,constr,factor*r2lambda,etot,e_di_tot,
     &dabs(etot-e_di_tot)/dabs(etot)*100.d0,epote*ma

    1 continue
    6 format(1x,i3,6(2x,e12.5))

      d2p = (r2(1)-r2lambda0)/co(1)      
c     write(2,*) r2(1),r2lambda0,co(1),d2p
c     write(2,*) 
c    &'From the 2-point forward formula, m-1 is: '
c    &,-d2p*ma/2.d0

c     d3p = (r2(2)-r2lambda0)/co(1)/2.d0      
c     write(2,*) d3p
c     write(2,*) 
c    &'From the 3-point forward formula, m-1 is: '
c    &,-d3p*ma/2.d0

      d3p_abr = (-r2(2)+4.d0*r2(1)-3.d0*r2lambda0)/co(1)/2.d0      
c     write(2,*) 
c    &'From the 3-point forward formula (ABR), der. and m-1 are: '
c    &,d3p_abr,-d3p_abr*ma/2.d0

c     d4p = (r2(3)/3.d0-r2(2)/2.d0+r2(1)-5.d0*r2lambda0/6.d0)/co(1)      
c     write(2,*) d4p
c     write(2,*) 
c    &'From the 4-point forward formula, m-1 is: '
c    &,-d4p*ma/2.d0

      d4p_abr = (2.d0*r2(3)-9.d0*r2(2)+18.d0*r2(1)
     &-11.d0*r2lambda0)/co(1)/6.d0      
c     write(2,*) 
c    &'From the 4-point forward formula (ABR), der. and m-1 are: '
c    &,d4p_abr,-d4p_abr*ma/2.d0

      d5p_abr = (-6.d0*r2(4)+32.d0*r2(3)-72.d0*r2(2)
     &+96.d0*r2(1)
     &-50.d0*r2lambda0)/co(1)/24.d0      
c     write(2,*) 
c    &'From the 5-point forward formula (ABR), der. and m-1 are: '
c    &,d5p_abr,-d5p_abr*ma/2.d0

      d6p_abr = (24.d0*r2(5)-150.d0*r2(4)+400.d0*r2(3)-600.d0*r2(2)
     &+600.d0*r2(1)
     &-274.d0*r2lambda0)/co(1)/120.d0      
c     write(2,*) 
c    &'From the 6-point forward formula (ABR), der. and m-1 are: '
c    &,d6p_abr,-d6p_abr*ma/2.d0

      e3p_abr = (ec(2)-2.d0*ec(1)+etot0)/co(1)**2      
c     write(2,*)
c    &'From the 3-point formula applied to the s.d. of E: ',
c    &e3p_abr,e3p_abr*ma/2.d0

      e5p_abr = (11.d0*ec(4)-56.d0*ec(3)+114.d0*ec(2)-104.d0*ec(1)
     &+35.d0*etot0)
     &/12.d0/co(1)**2      
c     write(2,*)
c    &'From the 5-point formula applied to the s.d. of E: ',
c    &e5p_abr,e5p_abr*ma/2.d0

      d3p=-d3p_abr*ma/2.d0
      d5p=-d5p_abr*ma/2.d0
      e3p=e3p_abr*ma/2.d0
      e5p=e5p_abr*ma/2.d0

      write(2,*)
      write(99,*) 
     &'m(-1) from:  <r^2> 3p    <r^2> 5p     E 3p        E 5p'
      write(99,5) -d3p_abr*ma/2.d0,-d5p_abr*ma/2.d0,
     &e3p_abr*ma/2.d0,e5p_abr*ma/2.d0 
      write(99,5) -d3p_abr*ma/2.d0/fpi,-d5p_abr*ma/2.d0/fpi,
     &e3p_abr*ma/2.d0/fpi,e5p_abr*ma/2.d0/fpi

c     amean=(d3p+d5p+e3p+e5p)/4.d0
      amean=(d3p+d5p)/2.d0
      err_d3p=(d3p-amean)/amean
      err_d5p=(d5p-amean)/amean
      err_e3p=(e3p-amean)/amean
      err_e5p=(e5p-amean)/amean
      write(99,*) amean
      write(99,5) err_d3p*100.d0,err_d5p*100.d0,
     &err_e3p*100.d0,err_e5p*100.d0

      write(2,*)
      write(2,*) '*** MONOPOLE ***'
      write(2,*)
      ewsr = 20.73d0*ma*4.d0*r2lambda0/16.d0/datan(1.d0)
      amean_to_comp_RPA = amean/16.d0/datan(1.d0)
      write(2,*) 'm(1),m(-1),sqrt[m(1)/m(-1)] = '
      write(2,*) ewsr,amean_to_comp_RPA,dsqrt(ewsr/amean_to_comp_RPA)
      akappaa = 939.d0*r2lambda0*(ewsr/amean_to_comp_RPA)/197.3d0**2
      write(2,*) 'A**(-1/3), K_c(A) = ',ma**(-1.d0/3.d0),akappaa
      write(2,*) 'm(3),m(1),sqrt[m(3)/m(1)] = '
      write(2,*) am3/16.d0/datan(1.d0),ewsr,
     1dsqrt(am3/16.d0/datan(1.d0)/ewsr)
      akappaas = 939.d0*r2lambda0*(am3/16.d0/datan(1.d0)/ewsr)
     1/197.3d0**2
      write(2,*) 'A**(-1/3), K_s(A) = ',ma**(-1.d0/3.d0),akappaas

      return
      end
