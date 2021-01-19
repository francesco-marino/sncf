      subroutine bcsgen
ccc   Program HARTREE-FOCK generalized to accept a pairing interaction
ccc   Generalized to have any number of points and step

c      use custom_interaction
      use ieee_arithmetic

      PARAMETER (nnp=10000,NEOC=1000)
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      double precision MA,MESMP(NNP),MESMN(NNP)                  
      common/b_to_do/ichf,ihf,irpa
      common/bdeb/ideb0,ideb
      COMMON/ENE/EMIN0,EMAX0,ECUT,NITER2
      COMMON/DER/FCT(NNP),DF(NNP),H                             
      COMMON RM(NNP),DU(NNP),VC(NNP),U(NNP),VS(NNP),TAUN(NNP),   
     +TAUP(NNP),TAUT(NNP),DJN1(NNP),DJP1(NNP),DJT1(NNP),VN(NNP),VP(NNP),
     +VSON(NNP),VSOP(NNP),ECOULB(NNP),VQRE(NNP),RM2(NNP),RM3(NNP),      
     +HMEN(NNP),HMEP(NNP),DHMEN(NNP),DHMEP(NNP),D2HMEN(NNP),D2HMEP(NNP),
     +DN(NNP),DP(NNP),DT(NNP),DN1(NNP),DP1(NNP),DT1(NNP),DN2(NNP),      
     +DP2(NNP),DT2(NNP),DJN(NNP),DJP(NNP),DJT(NNP),DC(NNP),AVP(NNP,     
     +NEOC),DUNL(NNP,NEOC),UNL(NNP,NEOC),VCD(NNP),DVC(NNP)   
C_GL_11_4_05
C     DVC introduced
      common/bmax/nmaxt
      common/b_box/ibox
      common/b_chf_3/constr
      common/b_chf_4/etot,r2lambda
      common/b_chf_5/am3
      common/b_chf_6/e_di_tot,epote
      common/b_mesh_r/step
      common/b_Skyrme/t0,t1,t2,t3,x0,x1,x2,x3,w0_0,w2p_0,t13,alfe
cc    HERE coefff. a la Dobaceski 
      common/b_Doba/c0rho,c0rho_dens,c1rho,c1rho_dens,c0deltarho,
     + c1deltarho,c0tau,c1tau,c0s,c0s_dens,c1s_dens,c0nablaj,c1nablaj,
     + c0deltas,c1deltas,c0t,c1t

      common/b_occup/nos,NPO,ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
     &ioc04,ideg4,ioc05,ideg5,ioc06,ideg6
      common/b_potential/icoul_d,icoul_e,iso,isj,icm1
      common/b_iter/nhf,niter
      common/b_extra/npplus,nnplus
      common/b_bcs/bcsp,bcsn
      common/b_vpair/vzero,xx,densc,gamma
      common/b_radii/r2l(12)
      common/bma/ma
      common/b_isotope/aa,zz
C_GL_17_6
      common/becutp/ecutp_p,ecutp_n
      common/bfo/fo(neoc)
      common/b_tensor/stc_t,stc_u
      common/b_tensor_i/itens
      common/b_dipo/DipoKap

ccc   constants
      common/bunits/pi,hbc,amc2,hbdm


      real*8 lampold,lamnold
      real*8 delpi,delni,lamp,lamn,G,gp,gpp,gn,gpn,z,a,SUMN,SUMP
      real*8 VZERO,densc,gamma,FAC,GGG,elrid1,posiz
      integer bcsp,bcsn,nhf     !variable pour calcul bcs
      dimension Gelem(neoc,neoc)
      DIMENSION densr(nnp),densrp(nnp),densrn(nnp)
      DIMENSION VPAIR(nnp),SOVR(NEOC,NEOC)     
      DIMENSION WF(20,nnp) 
      DIMENSION RMSM(nnp)                
      DIMENSION DELPOLD(NEOC),DELNOLD(NEOC)
      dimension deltag(neoc)
      DIMENSION deltap(NEOC),deltan(NEOC),dprep(NEOC),dpren(NEOC)
      DIMENSION elmar(NEOC,NEOC),VPRO(NEOC,NEOC),VNEU(NEOC,NEOC)
C_GL_22_4

      real(8), dimension(NNP) :: BETA(NNP)    !Beta: asymmetry (FRANCESCO)i
      real(8), dimension(NNP) :: so_dens      ! Spin-orbit density
      ! densities of direct and exchange coulomb  
      real(8), dimension(NNP) :: d_coul_dens
      real(8), dimension(NNP) :: e_coul_dens
      real(8), dimension(NNP) :: en_dens       ! Energy density
      real(8), dimension(NNP) :: kin_dens      ! kinetic energy density
      real(8), dimension(NNP) :: bulk_arr       ! potential energy density
      real(8), dimension(NNP) :: surf_dens     ! gradient terms
      real(8), dimension(NNP) :: cen_dens      ! central (t0123) energy density

ccc   isovector densities (DOba)
      DIMENSION D_IV(NNP),D_IV1(NNP),D_IV2(NNP), 
     +  DJ_IV(NNP),DJ_IV1(NNP), TAU_IV(NNP)
      DIMENSION NN(NEOC),LL(NEOC),LJ(NEOC),LT(NEOC),DEG(NEOC),EHF(NEOC),
     +EVSO(NEOC),DEG2(NEOC),EHF_OLD(NEOC)                     
      DIMENSION OV(NEOC,NEOC)                   
      DIMENSION REC(8),UNLOSC(NNP,NEOC)       
      DIMENSION DAL(3),MASH(3)               
      DIMENSION NX(231),LX(231),JX(231)
      dimension rr_ft(512),dens1_ft(512),dens2_ft(512)
      CHARACTER*6 ST(NEOC),STRY
      CHARACTER*6 SX(231) /' 1S1/2',' 1P3/2',' 1P1/2',' 1D5/2',' 2S1/2',
     &' 1D3/2',' 1F7/2',' 2P3/2',' 1F5/2',' 2P1/2',' 1G9/2',' 2D5/2',
     &' 1G7/2',' 3S1/2',' 2D3/2','1H11/2',' 2F7/2',' 1H9/2','1I13/2',
     &' 3P3/2',' 2F5/2',' 3P1/2','1I11/2',' 2G9/2','1J15/2',' 3D5/2',
     &' 4S1/2',' 2G7/2',' 3D3/2',' 3F5/2',' 3F7/2','2H11/2','1J13/2',
     &'1K17/2',' 2H9/2',' 4P3/2',' 4P1/2','1K15/2','2I13/2','2I11/2',
     &' 3G9/2',' 3G7/2',' 4D5/2',' 4D3/2',' 5S1/2','1L19/2','1L17/2',
     &'2J15/2','2J13/2','3H11/2',' 3H9/2',' 4F7/2',' 4F5/2',' 5P3/2',
     &' 5P1/2','1M21/2','1M19/2','2K17/2','2K15/2','3I13/2','3I11/2',
     &' 4G9/2',' 4G7/2',' 5D5/2',' 5D3/2',' 6S1/2','1N23/2','1N21/2',
     &'2L19/2','2L17/2','3J15/2','3J13/2','4H11/2',' 4H9/2',' 5F7/2',
     &' 5F5/2',' 6P3/2',' 6P1/2','1O25/2','1O23/2','2M21/2','2M19/2',
     &'3K17/2','3K15/2','4I13/2','4I11/2',' 5G9/2',' 5G7/2',' 6D5/2',
     &' 6D3/2',' 7S1/2',
     &'1Q27/2','1Q25/2','2N23/2','2N21/2','3L19/2',
     &'3L17/2','4J15/2','4J13/2','5H11/2',' 5H9/2',' 6F7/2',' 6F5/2',
     &' 7P3/2',' 7P1/2',
     &'1R29/2','1R27/2','2O25/2','2O23/2','3M21/2',
     &'3M19/2','4K17/2','4K15/2','5I13/2','5I11/2',' 6G9/2',' 6G7/2',
     &' 7D5/2',' 7D3/2',' 8S1/2',
     &'1T31/2','1T29/2',
     &'2Q27/2','2Q25/2','3N23/2','3N21/2','4L19/2',
     &'4L17/2','5J15/2','5J13/2','6H11/2',' 6H9/2',' 7F7/2',' 7F5/2',
     &' 8P3/2',' 8P1/2',
     &'1U33/2','1U31/2',
     &'2R29/2','2R27/2','3O25/2','3O23/2','4M21/2',
     &'4M19/2','5K17/2','5K15/2','6I13/2','6I11/2',' 7G9/2',' 7G7/2',
     &' 8D5/2',' 8D3/2',' 9S1/2',
     &'1V35/2','1V33/2',
     &'2T31/2','2T29/2',
     &'3Q27/2','3Q25/2','4N23/2','4N21/2','5L19/2',
     &'5L17/2','6J15/2','6J13/2','7H11/2',' 7H9/2',' 8F7/2',' 8F5/2',
     &' 9P3/2',' 9P1/2',
     &'1W37/2','1W35/2',
     &'2U33/2','2U31/2',
     &'3R29/2','3R27/2','4O25/2','4O23/2','5M21/2',
     &'5M19/2','6K17/2','6K15/2','7I13/2','7I11/2',' 8G9/2',' 8G7/2',
     &' 9D5/2',' 9D3/2','10S1/2',
     &'1Y39/2','1Y37/2',
     &'2V35/2','2V33/2',
     &'3T31/2','3T29/2',
     &'4Q27/2','4Q25/2','5N23/2','5N21/2','6L19/2',
     &'6L17/2','7J15/2','7J13/2','8H11/2',' 8H9/2',' 9F7/2',' 9F5/2',
     &'10P3/2','10P1/2',
     &'1Z41/2','1Z39/2',
     &'2W37/2','2W35/2',
     &'3U33/2','3U31/2',
     &'4R29/2','4R27/2','5O25/2','5O23/2','6M21/2',
     &'6M19/2','7K17/2','7K15/2','8I13/2','8I11/2',' 9G9/2',' 9G7/2',
     &'10D5/2','10D3/2','11S1/2'/
      DIMENSION VQ(nnp),UQ(nnp)                        
      DATA CP/.79577471d-01/,QP/12.566371/               
      DATA NX/1,1,1,1,2,1,1,2,1,2,1,2,1,3,2,1,2,1,1,3,2,3,1,2,1,3,4,2,3,
     &3,3,2,1,1,2,4,4,1,2,2,3,3,4,4,5,1,1,2,2,3,3,4,4,5,5,1,1,2,2,3,3,4,
     &4,5,5,6,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,7,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
     &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11
     &/
      DATA LX/0,1,1,2,0,2,3,1,3,1,4,2,4,0,2,5,3,5,6,1,3,1,6,4,7,2,0,4,2,
     &3,3,5,7,8,5,1,1,8,6,6,4,4,2,2,0,9,9,7,7,5,5,3,3,1,1,10,10,8,8,6,6,
     &4,4,2,2,0,11,11,9,9,7,7,5,5,3,3,1,1,12,12,10,10,8,8,6,6,4,4,2,2,0,
     &13,13,11,11,9,9,7,7,5,5,3,3,1,1,
     &14,14,12,12,10,10,8,8,6,6,4,4,2,2,0,
     &15,15,13,13,11,11,9,9,7,7,5,5,3,3,1,1,
     &16,16,14,14,12,12,10,10,8,8,6,6,4,4,2,2,0,
     &17,17,15,15,13,13,11,11,9,9,7,7,5,5,3,3,1,1,
     &18,18,16,16,14,14,12,12,10,10,8,8,6,6,4,4,2,2,0,
     &19,19,17,17,15,15,13,13,11,11,9,9,7,7,5,5,3,3,1,1,
     &20,20,18,18,16,16,14,14,12,12,10,10,8,8,6,6,4,4,2,2,0
     &/
      DATA JX/1,3,1,5,1,3,7,3,5,1,9,5,7,1,3,11,7,9,13,3,5,1,11,9,
     &15,5,1,7,3,5,7,11,13,17,9,3,1,15,13,11,9,7,5,3,1,19,17,15,13,11,9,
     &7,5,3,1,21,19,17,15,13,11,9,7,5,3,1,23,21,19,17,15,13,11,9,7,5,3,1
     &,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1,
     &41,39,37,35,33,31,29,27,25,23,21,19,17,15,13,11,9,7,5,3,1
     &/


c     FLAG   CUSTOM interaction (FRANCESCO)
      integer :: n_terms
      real(8), dimension(:), allocatable:: coeff_SM, coeff_NM, coeff_sym
      real(8), dimension(:), allocatable:: exponents

      logical :: FRANCESCO= .false.
      logical :: debug =  .true.   ! .false.
      character(Len = 90) :: dummy_string = ' '
      logical :: verbose = .false.   !.true. 


      real(8) :: i_ext = 0.        ! flag external potential
      real(8) :: omega = 3.        ! frequency, omega = [Energy]=MeV



  200 FORMAT(  //' ITERATION ', I3,'      XMU=',F5.2//)                 
  205 FORMAT(/'NOYAU DE MASSE A=',F9.1,3X,'ET DE CHARGE Z=',I5/)       
  209 FORMAT(8E10.3)                                                    
  211 FORMAT(10E12.5)                                                   
  212 FORMAT(10E12.5)                                                   
  213 FORMAT(//'  DENSITE NEUTRONS'/)                                  
  214 FORMAT(//'  DENSITE PROTONS '/)                                  
  216 FORMAT(//'  DENSITE CHARGE  '/)                                  
  217 FORMAT
     1(/' EC/A=',E15.8,' EHF/A=',E15.8,' EREA/A=',E15.8,' ETOT/A=',E15.8
     1//)       
  221 FORMAT(//' RN=',E13.5,' RP=',E13.5,' RC=',E13.5//)          
  350 format(75('*'))
 1200 FORMAT(//2X,'DENSITE DE MASSE'/)       
  408 format('ITER:',I3,3X,'I',I3,3X,'Energie=',E13.6,3X,
     &  'NI=',E13.6)
  409 format('lamp=',E13.6,3X,'lamn=',E13.6,3X,'Z=',E13.6,3X,'A=',E13.6)
 4090 format(75('*'))

      accur_bcs = 1e-6
      icount_chf=1
      if(ichf.eq.1.and.constr.gt.1e-10)icount_chf=2

      open (unit=10, file='fort.10', status='unknown')
      open (unit=16, file='dens.out', status='unknown')
      open (unit=56, file='dens1.out', status='unknown')
      open (unit=17, file='pot.out', status='unknown')
      open (unit=28, file='spur_dens.dat', status='unknown') 
      open (unit=29, file='dens_for_ft.dat',status='unknown')
      if(icount_chf.eq.1)then
      open (unit=12, file='lecpot.dat', status='unknown')       
      open (unit=27, file='facoc.in', status='unknown')                
      open (unit=50, file='occupied.dat', status='unknown')
      open (unit=51, file='unoccupied.dat', status='unknown')

      end if


ccc   Write the energy density
      open(unit=1997, file='energy_density.out', status='unknown')
      open(unit=1998, file='energy_density_full.out', status='unknown')
      open(unit=1897, file='contributions.out', status='unknown')


cc    Reading parameters          FRANCESCO
      open (unit=1996, file='interaction.in', status='old') 
      read(1996,'(A)') dummy_string
      read(1996,*) n_terms
      if (n_terms.ge.1) then
          FRANCESCO = .true.
          print *, "Using custom interaction "
      end if


      if (FRANCESCO) then

        allocate( exponents(1:n_terms) )
        allocate( coeff_SM(1:n_terms) )
        allocate( coeff_NM(1:n_terms) )
        allocate( coeff_sym(1:n_terms) )

        read(1996,'(A)') dummy_string
c       read exponents & coefficients
        do i=1,n_terms
           read(1996,*) exponents(i), coeff_SM(i), coeff_NM(i)
           exponents(i) = exponents(i)/3.0d0
           coeff_sym(i) = coeff_NM(i) - coeff_SM(i)
        end do
c       gradient terms 
        read(1996,'(A)') dummy_string 
        read(1996,*) c0deltarho
        read(1996,*) c1deltarho

c       nucleus (A, Z)
        read(1996,'(A)') dummy_string
        read(1996,*) aa, zz

ccc     energy harmonic trap
        read(1996,'(A)') dummy_string
        read(1996,*) omega

        if (omega.eq.0.) then
           i_ext=0.
        end if

c       spin orbit (yes or no)
        read(1996,'(A)') dummy_string
        read(1996,'(A)') dummy_string
        if (index(dummy_string,'no spin orbit').ne.0) then
           iso = 0
        else
           iso = 1
           
        end if
cc      Coulomb        
        if (index(dummy_string,'no coul').ne.0) then
           icoul_d=0
           icoul_e=0
           print *, "No Coulomb"
        else
           icoul_d=1
           icoul_e=1
        end if

      end if 


      if (i_ext.ne.0) then
         print *, "Harmonic trap of energy E=", omega," MeV"
      end if 



cc    read from command line
      if (command_argument_count().ge.2) then   ! read A and Z from command line
          call get_command_argument(1, dummy_string)
          read(dummy_string,*) aa    !convert ro real
          call get_command_argument(2, dummy_string)
          read(dummy_string,*) zz        

          if (command_argument_count().ge.4) then
             call get_command_argument(3, dummy_string)
             read(dummy_string,*) c0deltarho
             call get_command_argument(4, dummy_string)
             read(dummy_string,*) c1deltarho        
          end if
ccc       spin-orbit
          if (command_argument_count().ge.5) then
             call get_command_argument(5, dummy_string)
             read(dummy_string,*) w0_0
             w2p_0=w0_0        ! assuming the 2 spin strengths are equal
             c0nablaj = -0.75d0 * w0_0
             c1nablaj = -0.25d0 * w0_0
             if (w0_0.ne.0) then   
                iso= 1
             end if
          end if  

        end if
        
        if (FRANCESCO) then
          print *, "C0 = ", c0deltarho
          print *, "C1 = ", c1deltarho
        end if
        if (iso.eq.1) then
           print '(3(A8,F8.2))',"W",w0_0,"cj0",c0nablaj,"cj1",c1nablaj 
        else         
           print *, "No spin orbit "
        end if


CCC
CCC    WRITING PARAMETERS                                              
CCC
      i_grnp0_box=1
      write(2,*) 'Number of points and mesh: ',nmaxt,step
      w2 = w0_0 
      w1 = 2.d0*w0_0
      w2p = w2p_0
      w1p = 2.d0*w2p_0
      IF(ABS(ALFE).LT.1.E-04) ALFE=1.           
      ni = niter
      kop1 = 2
      ipchpo = 1 
      ipchwf = 0
      idcc = 0
      ipunch = 0
      lamp = -7.d0
      lamn = -7.d0
      delpi = 0.d0
      delni = 0.d0

      write(2,*)
      write(2,*) 'Initial guess for lambda_p,lambda_n, 
     &parameters bcsp,bcsn: '
      WRITE(2,*) lamp,lamn,bcsp,bcsn
      if(bcsp.ne.0)then 
       delpi=2.d0
       write(2,*) 'Initial guess for delta_p: ',delpi
      end if
      if(bcsn.ne.0)then 
       delni=2.d0
       write(2,*) 'Initial guess for delta_n: ',delni
      end if
      write(2,*) 'Pairing force: V0, x, rho0 and exponent: '
      WRITE(2,*) VZERO,xx,densc,gamma
      write(2,*) 'Max n. iterations: ',niter
      write(2,*) 

c     HERE some options      
      if(icoul_d.eq.0)write(2,*) 'NO COULOMB'
      if(icoul_e.eq.0)write(2,*) 'NO COULOMB EXCHANGE'
      if(iso.eq.0)write(2,*) 'NO SPIN-ORBIT'
      if(isj.eq.0)write(2,*) 'NO J2 TERMS'
      if(icm1.eq.0)write(2,*) 'NO CENTER OF MASS CORRECTION'
      if(icoul_d.eq.1.and.icoul_e.eq.1.and.iso.eq.1.and.isj.eq.1.and.
     &icm1.eq.1)write(2,*) 'FULL HF POTENTIAL INCLUDED'
      if(itens.eq.1)write(2,*) 'TENSOR INCLUDED; T, U = ',stc_t,stc_u
      NNN=NMAXT
      if(icount_chf.eq.1)nnn2=nnn-2
      DO J=1,NNN                                                   
       U(J)=0.    
      end do        
      H=STEP

      do ir=1,NNN
       densr(ir)=densrp(ir)+densrn(ir)
      enddo
      
      IPRE = 2                              
      DO J=1,NNN                        
       RM(J) = H*FLOAT(J)                    
       RM2(J) = RM(J) * RM(J)                
       RM3(J) = RM(J)*RM2(J)                 
      end do                             
      UNS4PI = CP                         
      NN1 = NNN- 1                      
      NN2 = NNN- 2                     
      NN3 = NNN- 3                    
      NN4 = NNN- 4                   
c_CHF
      if(ichf.eq.1)iadd=1
      if(ichf.eq.1.and.irpa_chf.eq.3.and.iter.ne.itmax)iadd = 0
c_CHF
      W3 = W2                       
      W0 = W1
c0995
      W3P = W2P
      W0P = W1P
c0995
      CRO = .25*(T1*(1.+.5*X1)+T2*(1.+.5*X2))                      
      CROQ = .125*(T2*(1.+2.*X2) -T1*(1.+2.*X1))                   
      TXRO = T0*(1. + .5*X0)                                       
      TXROQ = -T0*(X0 +.5)                                         
      TXDRO = .5*(T2*(1.+.5*X2)-T1*(1.+.5*X1))                     
      TXD2RO = .125*(T2*(1.+.5*X2) - 3.*T1*(1.+.5*X1))             
      TXDROQ = .25*(T1*(1.+2.*X1) + T2*(1.+2.*X2))                 
      TXD2RQ = .0625*(3.*T1*(1. +2.*X1) + T2*(1.+2.*X2))           
      TXTAU = CRO                                                 
      TXTAUQ = CROQ                                              
      C13=T13/24.                                               

c_order_levels
      ma = aa
      jz = zz
      ntp = 0
      do 1011 icp=1,231
      deg(icp) = 0.d0
      nn(icp) = nx(icp)
      ll(icp) = lx(icp)
      lj(icp) = jx(icp)
      lt(icp) = 1
      st(icp) = sx(icp)
      ntp = ntp + lj(icp) + 1
      deg(icp) = lj(icp) + 1
      if(ntp.ge.int(jz))go to 1012
 1011 continue
 1244 stop '>>> THE REQUESTED NUCLEUS HAS TOO MANY PROTON LEVELS'
 1012 NPO = icp
      if(ntp.gt.int(jz))deg(icp)=deg(icp)-(float(ntp)-jz)
      if(npplus.ne.0)then
       NPO=NPO+npplus
       do iextra=icp+1,NPO
       if(iextra.gt.231)go to 1244
       deg(iextra) = 0.d0
       nn(iextra) = nx(iextra)
       ll(iextra) = lx(iextra)
       lj(iextra) = jx(iextra)
       lt(iextra) = 1
       st(iextra) = sx(iextra)
       end do
      end if

      ntn = 0
      do 1013 icn=1,231
      icn1 = icn + NPO
      deg(icn1) = 0.d0
      nn(icn1) = nx(icn)
      ll(icn1) = lx(icn)
      lj(icn1) = jx(icn)
      lt(icn1) = 0
      st(icn1) = sx(icn)
      ntn = ntn + lj(icn1) + 1
      deg(icn1) = lj(icn1) + 1
      if(ntn.ge.int(ma-jz))go to 1014
 1013 continue
 1245 stop '>>> THE REQUESTED NUCLEUS HAS TOO MANY NEUTRON LEVELS'
 1014 nos = icn1
      if(ntn.gt.int(ma-jz))deg(icn1)=deg(icn1)-(float(ntn)-ma+jz)
      if(nnplus.ne.0)then
       nos=nos+nnplus
       do iextra=icn+1,icn+nnplus
       if(iextra.gt.231)go to 1245
       iextra1 = iextra + NPO
       deg(iextra1) = 0.d0
       nn(iextra1) = nx(iextra)
       ll(iextra1) = lx(iextra)
       lj(iextra1) = jx(iextra)
       lt(iextra1) = 0
       st(iextra1) = sx(iextra)
       end do
      end if

      if(icount_chf.eq.1)then
      write(2,*)
      write(2,*) 'Added (total) proton and neutron levels: '
      write(2,*) npplus,NPO
      write(2,*) nnplus,nos
C_GL_17_6
      write(2,*) 'BCS cutoff energies: '
      write(2,*) ecutp_p,ecutp_n
      write(2,*)
      end if

      z_ver = 0.
      do i=1,NPO
      z_ver = z_ver + deg(i)
      end do

      amz_ver = 0.
      do i=NPO+1,nos
      amz_ver = amz_ver + deg(i)
      end do
    
      if(dabs(z_ver-jz).gt.1e-3.
     1   or.dabs((z_ver+amz_ver)-ma).gt.1e-3)then
       write(2,*)
       write(2,*) '>>> WRONG COUNT OF PARTICLES: '
       write(2,*) jz,z_ver,dabs(z_ver-jz)
       write(2,*) ma,z_ver+amz_ver,dabs((z_ver+amz_ver)-ma)
       write(2,*)
       stop
      end if

      do I=1,NPO
      deltap(I)=delpi
      dprep(I)=delpi
      enddo
      do I=NPO+1,NOS
      deltan(I)=delni
      dpren(I)=delni
      enddo


      do I=1,NPO
      do k=1,NPO
      GGG=0.d0
C_GL_2_4_05
      lcount_min=iabs(lj(i)-lj(k))/2 
      lcount_max=(lj(i)+lj(k))/2
c     do lcount=0,13
      do lcount=lcount_min,lcount_max
      MULTIP=lcount 
      elrid1=yl(LL(I),LJ(I),LL(K),LJ(K),MULTIP)
      GGG=GGG+elrid1**2
      enddo
      elmar(I,K)=GGG
      enddo
      enddo
C     parte di neutroni
      do I=NPO+1,NOS
         do k=NPO+1,NOS
         GGG=0.d0
C_GL_2_4_05
         lcount_min=iabs(lj(i)-lj(k))/2 
         lcount_max=(lj(i)+lj(k))/2
c     do lcount=0,13
         do lcount=lcount_min,lcount_max
            MULTIP=lcount
            elrid1=yl(LL(I),LJ(I),LL(K),LJ(K),MULTIP)
            GGG=GGG+elrid1**2
         enddo
         elmar(I,K)=GGG
         enddo
      enddo


c     if(icount_chf.eq.1) then
       write(2,4090) 
       write(2,205) MA,JZ                        
c      write(2,*) NPO,nos
       write(2,4090) 
c     end if

      if(icm1.eq.0) DMSHB = .04823    
      hb0=1./.04823
      if(icm1.eq.1) DMSHB = .04823*MA/(MA-1)                    
cc    hbar^2/(2m)      
      hb = 1./DMSHB                               
      if(icm1.eq.2) then
       vt_butler = (3.d0*ma/2.d0)**(1./3.)
       ff_butler = 2.d0/(vt_butler+1.d0/3.d0/vt_butler) 
       write(2,*)
       write(2,*) 'f(A) in the paper by Butler et al. : ',ff_butler
       write(2,*)
       DMSHB = .04823*(1-ff_butler/ma)
       hb = 1./DMSHB
      end if
CCC                                                                             
CCC    PUITS DE SAXON WOODS                                                     
CCC                                                                             
      R0 = 1.21
      A1 = 55.5
      A2 = 33.2                                
      A3 = .36                                
      A5 = .68                               
      CWPI = 1.41                           
      CW2 = CWPI*CWPI                      
      AL = A5                             
      R = R0*MA**(1./3.)    

CCC                                                                             
CCC   STARTING  VALUE                                                           
CCC                                                                             
      ITER = 0                                                        
      GO TO(137,341,601), KOP1                                       
  137 if(icount_chf.eq.1)read 209, ((AVP(J,I),J=1,NNN),I=1,NOS)
     &,(EHF(I),I=1,NOS)    
      GO TO 107                                                      
  341 CONTINUE     

      DO J=1,NNN                                               
        HMEN(J) = HB                                                
        HMEP(J) = HB                                               
        DHMEN(J) = 0.d0                                                 
        DHMEP(J) = 0.d0  
        D2HMEN (J) = 0.d0                                              
        D2HMEP (J) = 0.d0                                             
        MESMN(J) = 1.d0     
        MESMP(J) = 1.d0                                              
      end do   
                                                   
      KX = 1                                                     
      E = 0.d0
ccc   Electron charge^2 [MeV fm]   here 
      qe2 =  1.43986d0     
      E22 = 1.43986d0*H
      F = 0.d0                                                 
      RC = R*1.09d0/R0                                        
      SYM = (MA-2*JZ)/MA                                     
      Z = 0.d0                                              
      DO K = 1,NNN
        X = RM(K)                                            
        DP(K) = 1.d0/(1.d0 + EXP((X-RC)/0.55d0)) 
        Z = Z + X*X*DP(K)
      end do                    
      Y = FLOAT(JZ)/(QP*Z*H)                             
      DO K = 1,NNN                                  
        DP(K) = Y*DP(K) 
      end do                                
      DO K=1,NNN                                  
       Y = QP*RM(K)*RM(K)                            
       E = E + Y*DP(K)                              
       F = F + Y*DP(K)/RM(K)                       
       VC(K) = E/RM(K) - F
      end do                        
      FF = F                                    
      DO J=1,NNN                          
       VC(J) = E22*(VC(J) + FF)                
       X = RM(J)                              
       PE = EXP((X-R)/AL)                    
       FPE = -1./(1.+PE)                    
       VN(J) = FPE                         
       VS(J) = - PE*FPE*FPE/AL            
       DO I=1,NOS                    
        N = NN(I)                        
        L = LL(I)                       
        JM = LJ(I)                     
        FL = (JM*(JM+2) - 4*L*(L+1) - 3)/8.                          
        ISIG = 2*LT(I) - 1                                           
        V = A1 + ISIG*A2*SYM                                         
        AVP(J,I)=V*VN(J)+A3*CW2*VS(J)*2.*FL*float(iso)
        IF (I.LE.NPO) then 
         AVP(J,I) = AVP(J,I) + float(icoul_d)*VC(J)
        end if
       end do  
      end do

      DO I=1,NOS                                           
       EHF(I) = .3 * AVP(1,I)                                  
       IF(I.LE.NPO) then 
        EHF(I) = EHF(I) +.3*float(icoul_d)*VC(1) 
       end if
      end do              
  601 CONTINUE                                              
  107 CONTINUE                                            
      ITMAX = NI  


CCC                                                                             
CCC   STARTING ITERATION                                                         
CCC
      XMU=0.5
c     THIS IS THE FAST ITERATION PROCEDURE
c     XMU=0.15
c     THIS IS THE SLOW ITERATION PROCEDURE
      XMU=0.99     

  104 ITER = ITER + 1                                  

c     THIS IS THE FAST ITERATION PROCEDURE
c     IF(ITER.EQ.5) XMU=0.3
c     IF(ITER.EQ.10) XMU=0.5
c     IF(ITER.EQ.15) XMU=0.7
c     IF(ITER.EQ.20) XMU=0.85

c     THIS IS THE SLOW ITERATION PROCEDURE
      IF(ITER.EQ.5)  XMU=0.95
      IF(ITER.EQ.10) XMU=0.93                             
      IF(ITER.EQ.15) XMU=0.90                             
      IF(ITER.EQ.20) XMU=0.85
      IF(ITER.EQ.30) XMU=0.80
      IF(ITER.EQ.40) XMU=0.70
      IF(ITER.EQ.50) XMU=0.60
      IF(ITER.EQ.70) XMU=0.50
      IF(ITER.EQ.90) XMU=0.40
      IF(ITER.EQ.100) XMU=0.30
      IF(ITER.EQ.110) XMU=0.20
      IF(ITER.EQ.120) XMU=0.10
          
      IF(ITER.EQ.ITMAX) WRITE(2,200) ITER,XMU

      IF(ichf.eq.1.and.ITER.EQ.1) then
       write(2,*) 
       write(2,*) 'Constr. ',constr
       write(2,*)
      end if


CCC                                                                             
CCC   SOLVING SCHROEDINGER EQUATION                                  
CCC   
      
      EPOTE=0.         ! sum of energy eigenvalues                                 
      DO 110  I=1,NOS                                  
c     write(2,*) i,nn(i),ll(i),lj(i),lt(i),ehf(i)
      N = NN(I)                                         
      L = LL(I)                                         
      LL1 = L*(L+1)                                    
      EI = EHF(I)                                     
      K = 0                                          
      S2 = 0.                                       
    
      ESO=(LJ(I)*(LJ(I)+2.)-4.*L*(L+1.)-3.)/8.     ! see Vautherin&Brink

      DO 300 J=1,NNN   

       IF(LT(I).EQ.0) then
        RMSM(J)=(1.-MESMN(J))
        VQRE(J) = MESMN(J) * (AVP(J,I) - .25 * DHMEN(J) *         
     1  DHMEN(J) / HMEN(J) + .5*D2HMEN(J) )                    
       endif
       IF(LT(I).EQ.1) then
        RMSM(J)=(1.-MESMP(J))
     
        VQRE(J) = MESMP(J) * (AVP(J,I) - .25 * DHMEP (J)*         
     1  DHMEP(J) / HMEP(J) + .5*D2HMEP(J) )                    
       endif

c_CHF
       VQRE(J) = VQRE(J)+float(iadd)*CONSTR*RM(J)*RM(J)
       write(98,*) 
     1 float(j)*h,VQRE(J)-float(iadd)*CONSTR*RM(J)*RM(J),vqre(j)
c_CHF
       VQRE(J) = VQRE(J)/HB + LL1/(RM(J)*RM(J))                    

  300 CONTINUE                                                     
      GO TO 112       

  113 K=K+2                                                      
      IF (K-nnn) 114,114,115                                 
  115 WRITE(2,*) I,NN(I),LL(I),LJ(I),LT(I)                 
      GO TO 111     

  114 EI = -.5*FLOAT(K)                                   
  112 continue
      E = EI * DMSHB                                      
      
      NO = N - 1                                         

      ideb=ideb0

ccc   solve SCHRODINGER ???
      CALL NUM2B(NNN-1,H,E,S2,VQRE,RMSM,U,NO,HB)
      IF (NO.LT.0) GO TO 113

  111 CONTINUE                                
      U(NNN)=0.
      DO 12 J=1,NNN                          
   

      FCT(J) = U(J)
   12 continue                         
      CALL DERIV                           
      S = 0.                              
      DO J= 1,NNN                   
        DU(J)  = DF(J)                    
        S = S + U(J)*U(J)
      end do

      S = SQRT(S*H)                   
      SIG = U(1)/ABS(U(1))           
      SIG = 1.
      DO J=1,NNN               
        DU(J) = SIG*DU(J)/S        
        U(J) = SIG*U(J)/S    
      end do       
      E = E/DMSHB
      EHF(I)=E
      DO J=1,NNN             
       DUNL(J,I) = DU(J)         
       UNL(J,I) = U(J)
      end do                            
      XKZ=0.d0                                                  
      DO J=1,NNN                                            
        ZY=(-UNL(J,I)/RM(J)+DUNL(J,I))/RM(J)                      
        ZY1=UNL(J,I)/RM2(J)                                       
        XKZ=XKZ+(ZY*ZY+LL1*ZY1*ZY1)*RM2(J) 
      end do                       
      XKZ=SQRT(H*XKZ)                                           
      X = (E-EI)/E                                              
      X = ABS(X)                                               


 2221 FORMAT(2X,A6,I2,5X,'E=',E11.4,5X,'DEG=',E11.4,5X,' D=', 
     1    F8.4,5X,' NI=',E11.4,' s.o.(e.m.)=',E11.4)                  
      EPOTE=EPOTE+E*DEG(I)                                  
  110 CONTINUE                                             

      if ((bcsp .eq. 0) .and. (bcsn .eq. 0 )
     & .or. nhf .ge. ITER )goto 1234  !calcul bcs demande 
      
c     open(unit=22,status='unknown',file='vpair.dat')
      itbcs=itbcs+1
      do ir=1,NNN
        VPAIR(ir)=VZERO*(1.d0-xx*(DT(ir)/densc)**gamma)
        potpair=-vpair(ir)
c     write(22,*)ir,potpair
      enddo
c     close(22)


C     parte di protoni
      DO i=1,NPO
       do k=1,NPO
        SOVR(i,k)=0.d0
        do ir=1,NNN
         SOVR(i,k)=SOVR(i,k)+VPAIR(ir)*(unl(ir,i)**2)*
     1          (unl(ir,k)**2)/ir/ir/step
        enddo
       enddo
      enddo
      
C     parte di neutroni
      DO i=NPO+1,NOS
       do k=NPO+1,NOS
        SOVR(i,k)=0.d0
        do ir=1,NNN
         SOVR(i,k)=SOVR(i,k)+VPAIR(ir)*(unl(ir,i)**2)*
     1          (unl(ir,k)**2)/ir/ir/step
        enddo
       enddo
      enddo
   
C     fine elemento matrice radiale

      do i=NPO+1,NOS
        do j=NPO+1,NOS
C_GL_17_6_05
        vneu(i,j)=0.d0
        if(ehf(i).lt.ecutp_n.and.ehf(j).lt.ecutp_n)then
          VNEU(i,j)=SOVR(i,j)*elmar(i,j)
        end if

         write(85,350)
         write(85,*) 'Neutrons'
         write(85,*)ll(i),lj(i),ehf(i)
         write(85,*)ll(j),lj(j),ehf(j)
         write(85,*)sovr(i,j),elmar(i,j),2.*vneu(i,j)
   
        enddo
      enddo


      do i=1,NPO
        do j=1,NPO
C_GL_17_6_05
        vpro(i,j)=0.d0
        if(ehf(i).lt.ecutp_p.and.ehf(j).lt.ecutp_p)then
           VPRO(i,j)=SOVR(i,j)*elmar(i,j)/dsqrt(LJ(i)+1.d0)/
     1          dsqrt(LJ(j)+1.d0)
        end if

        write(85,350)
        write(85,*) 'Protons'
        write(85,*)ll(i),lj(i),ehf(i)
        write(85,*)ll(j),lj(j),ehf(j)
        write(85,*)sovr(i,j),elmar(i,j),2.*vpro(i,j)

        enddo
      enddo


CCC   PARTIE BCS
      
      if (bcsp .eq. 1) then
        do I=1,NPO
        DEG(I)=LJ(I)+1
        enddo
      endif
      
      if (bcsn .eq. 1) then
        do I=NPO+1,NOS
        DEG(I)=LJ(I)+1
        enddo
      endif


CCC  1) CHEMICAL POTENTIAL

      if (abs(lamn).ge.0.0000001) goto  743

      test0p=-5.
      test0n=-5.

c     premier test pour voir si ca diverge avec les valeurs de lambda


c     calcul de gp,gpp
 874  testp=test0p
      testn=test0n
      do icount=1,2   !nbe iterations pour teste si le lambda est bon
      gn =0.
      gp =0.
      gpp =0.
      gpn =0.
      do I=1,NPO
        gp=gp+((DEG(I)/2.)*(1.-((EHF(I)-testp)/((EHF(I)-testp)**2
     &    +dprep(I)**2)**0.5)))
        gpp=gpp+DEG(I)/2.*(dprep(I)**2/((EHF(I)-testp)**2
     &     +dprep(I)**2)**1.5)
      enddo
c     calcul de gn,gpn
      do I=NPO+1,NOS
        gn=gn+((DEG(I)/2.)*(1.-((EHF(I)-testn)/((EHF(I)-testn)**2
     &    +dpren(I)**2)**0.5)))
        gpn=gpn+DEG(I)/2.*(dpren(I)**2/((EHF(I)-testn)**2
     &     +dpren(I)**2)**1.5)
      enddo  

      testp=testp+(float(JZ)-gp)/gpp
      testn=testn+(MA-float(JZ)-gn)/gpn

      enddo

c     write(6,*) test0p,testp,test0n,testn

      if((testp.le.-20).or.(testp.ge.0)) test0p=test0p-0.1
      if((testn.le.-20).or.(testn.ge.0)) test0n=test0n-0.1
      if((testp.le.-20).or.(testp.ge.0).or.
     &(testn.le.-20).or.(testn.ge.0))  goto 874 

      lamp=testp    !valeur de lambda acceptee pour algo
      lamn=testn

 743  continue
      write(2,*) 'lambdap et n initiaux:',lamp,lamn

CC    ciclo per la soluzione BCS

      do 3007 icountbcs=1,30

CC    ciclo per la soluzione della GAP EQ. per i lambda attuali

cc    LA GAP EQUATION è PER ELEMENTO DI MATRICE ANTISIMMETRIZZATO
   
cc    delta des etats de proton
      
        if (bcsp.eq.0) goto 3001
    
        do i=1,NPO
           delpold(i)=dprep(i)
        enddo
 
      do 3004 igappro=1,10000

      do I=1,NPO
        deltap(I)=0

        do J=1,NPO
     
        SUMP=dprep(J)*VPRO(I,J)/(2.d0*((EHF(J)-lamp)
     1     **2+dprep(J)**2)**0.5)*dsqrt(LJ(j)+1.d0)/
     1     dsqrt(LJ(i)+1.d0)
        deltap(I)=deltap(I)+SUMP
      
        enddo
      enddo

cc    check di convergenza gap eq.
      chdel=0
      do i=1,NPO
        diffdel=abs(deltap(i)-dprep(i))
        if(diffdel.gt.accur_bcs) chdel=chdel+1
      enddo
      if(chdel.eq.0) goto 3003      

cc    nuova assegnazione
      do i=1,NPO
         dprep(i)=deltap(i)
      enddo

 3004 continue

 3003 continue

cc    nuova assegnazione
      do i=1,NPO
        dprep(i)=deltap(i)
      enddo
 
cc    calcolo delta medio protonico
      somdeg=0
      delmedp=0
      do i=1,NPO
C_GL_17_6
        if(ehf(i).gt.(lamp-5.).and.ehf(i).lt.(lamp+5.))then
          somdeg=somdeg+LJ(I)+1.
          delmedp=delmedp+deltap(I)*(LJ(i)+1.)
        end if
      enddo
      deltamediop=delmedp/somdeg
 
 3001 continue


c     delta des etats de neutron      
      if (bcsn.eq.0) goto 3002
    
      do i=NPO+1,nos
        delnold(i)=dpren(i)
      enddo

      do 3006 igapneu=1,10000
     
      do I=NPO+1,NOS
        deltan(I)=0

        do J=NPO+1,NOS

        SUMN=dpren(J)
     1     *VNEU(I,J)/float(lj(i)+1) 
     1     /(2.d0*((EHF(J)-lamn)**2+dpren(J)**2)**0.5)
        deltan(I)=deltan(I)+SUMN
 
        enddo
      enddo    


cc    check di convergenza gap eq.
      chdel=0
      do i=NPO+1,nos
        diffdel=abs(deltan(i)-dpren(i))
        if(diffdel.gt.accur_bcs) chdel=chdel+1
        enddo
      if(chdel.eq.0) goto 3005      

cc    nuova assegnazione
      do i=NPO+1,nos
        dpren(i)=deltan(i)
      enddo

 3006 continue

 3005 continue

cc    nuova assegazione
      do i=NPO+1,nos
        dpren(i)=deltan(i)
      enddo

cc    calcolo delta medio neutronico
      somdeg=0
      delmedn=0

      do i=NPO+1,nos
C_GL_17_6
c     write(2,*) i,ehf(i),lamn
      if(ehf(i).gt.(lamn-5.).and.ehf(i).lt.(lamn+5.))then
       somdeg=somdeg+LJ(I)+1.
       delmedn=delmedn+deltan(I)*(LJ(i)+1.)
c      if(iter.eq.itmax-1)then
c        write(2,*) i,ehf(i),lj(i),deltan(i)
c        write(2,*) somdeg,delmedn
c      end if 
      end if
      enddo
      deltamedion=delmedn/somdeg
      do I=NPO+1,NOS
        dpren(I)=deltan(I)
      enddo
 3002 continue

cc    fine ciclo gap eq. neutr.


CC    Ciclo per la NUMBER EQ.

      lampold=lamp
      lamnold=lamn

      do icount=1,10000          !iterations algo de Newton-Raph
      
      gn =0.
      gp =0.
      gpp =0.
      gpn =0.
      
      if(bcsp.eq.0) goto 3009
c     calcul de gp,gpp
    
        do I=1,NPO
      
        gp=gp+((DEG(I)/2.0d0)*(1.0d0-((EHF(I)-lamp)/
     &     ((EHF(I)-lamp)**2+deltap(I)**2)**0.5)))
        gpp=gpp+(DEG(I)/2./(((EHF(I)-lamp)**2
     &     +deltap(I)**2)**0.5))

        enddo

      lamp=lamp+(float(JZ)-gp)/gpp

 3009 continue

        if(bcsn.eq.0) goto 3010
c     calcul de gn,gpn
        do I=NPO+1,NOS
         gn=gn+((DEG(I)/2.0d0)*(1.0d0-((EHF(I)-lamn)/
     &     ((EHF(I)-lamn)**2+deltan(I)**2)**0.5)))
         gpn=gpn+(DEG(I)/2./(((EHF(I)-lamn)**2
     &      +deltan(I)**2)**0.5))
       
         enddo  
     
       lamn=lamn+(MA-float(JZ)-gn)/gpn

 3010 continue
         
      enddo
CC    fine ciclo number equation    

cc    test per chiudere iterazione interna bcs
      chdel=0
      do i=1,NPO
      diffdel=abs(delpold(i)-deltap(i))
      if(diffdel.gt.accur_bcs) chdel=chdel+1
      enddo
      do i=NPO+1,nos
      diffdel=abs(delnold(i)-deltan(i))
      if(diffdel.gt.accur_bcs) chdel=chdel+1
      enddo
      difflamp=abs(lampold-lamp)
      difflamn=abs(lamnold-lamn) 
      if(difflamn.lt.accur_bcs.and.difflamp
     1   .lt.accur_bcs.and.chdel.eq.0) goto 3008 

 3007 continue   
 3008 continue

CC    fine ciclo bcs


      write(2,*)'iter=',iter
      write(2,*)'lamp,lamn',lamp,lamn
      do i=1,NPO
      write(2,*)'delta  ',i,deltap(i)
      enddo
      do i=NPO+1,nos
      write(2,*)'delta  ',i,deltan(i)
      enddo
c     fin du calcul des deltas
     
      if(itbcs.eq.1) then   
        do i=1,NPO
        eqp=((ehf(i)-lamp)**2+deltap(i)**2)**.5
        enddo
        do i=NPO+1,nos
        eqp=((ehf(i)-lamn)**2+deltan(i)**2)**.5
        enddo
        do i=NPO+1,nos
          do j=NPO+1,nos
          Gelem(i,j)=vneu(i,j)*2/sqrt(lj(i)+1.)/sqrt(lj(j)+1.)
          enddo
        enddo
      endif


c     2) OCCUPATION FACTORS

      if (bcsp .eq. 1) then
      do I=1,NPO
        FO(I)=0.5*(1-(EHF(I)-lamp)/((EHF(I)-lamp)**2
     &  +deltap(I)**2)**0.5)
        DEG(I)=DEG(I)*FO(I)
      enddo
      else
        do i=1,NPO
          fo(i)=deg(i)/(lj(i)+1)  !si pas de BCS facteur d'occupation
         enddo 
      endif
      if (bcsn .eq. 1) then
      do I=NPO+1,NOS
        FO(I)=0.5*(1-(EHF(I)-lamn)/((EHF(I)-lamn)**2 +
     &    deltan(I)**2)**0.5)
        DEG(I)=DEG(I)*FO(I)
      enddo
      else
        do i=NPO+1,nos
          fo(i)=deg(i)/(lj(i)+1)
         enddo 
      endif
c     Recalcul du nombre de nucleons

      z=0.
      a=0.
      
      do I=1,NPO
        z=z+(LJ(I)+1)*FO(I)
      enddo  

      do I=1,NOS
        a=a+(LJ(I)+1)*FO(I)
      enddo  

      DO I=1,NOS
        write(23,408) ITER,I,EHF(I),FO(I)
        enddo
        write(23,350) 
        write(23,409)lamp,lamn,z,a 
        write(23,350) 
      
 1234   continue

c_order_levels

c     Rearrangement des niveaux par energie

      ibench=0
      if(ibench.eq.1)go to 1241

      DO 1235 I=1,NPO
      DO 1236 J=I+1,NPO
      IF(EHF(J).LE.EHF(I))THEN
       P=EHF(I)
       EHF(I)=EHF(J)
       EHF(J)=P
       NTRY=NN(I)
       NN(I)=NN(J)
       NN(J)=NTRY
       LTRY=LL(I)
       LL(I)=LL(J)
       LL(J)=LTRY
       JTRY=LJ(I)
       LJ(I)=LJ(J)
       LJ(J)=JTRY
       if(bcsp.eq.1) then
        P=DEG(I)
        DEG(I)=DEG(J)
        DEG(J)=P
        P=FO(I)
        FO(I)=FO(J)
        FO(J)=P
       end if
       EVTRY=EVSO(I)
       EVSO(I)=EVSO(J)
       EVSO(J)=EVTRY
       STRY=ST(I)
       ST(I)=ST(J)
       ST(J)=STRY
       p=deltap(i)
       deltap(i)=deltap(j)
       deltap(j)=p
       p=dprep(i)
       dprep(i)=dprep(j)
       dprep(j)=p
       DO 1237 J1=1,NNN
         P=UNL(J1,I)
         UNL(J1,I)=UNL(J1,J)
         UNL(J1,J)=P
         P=DUNL(J1,I)
         DUNL(J1,I)=DUNL(J1,J)
         DUNL(J1,J)=P
         P=AVP(J1,I)
         AVP(J1,I)=AVP(J1,J)
         AVP(J1,J)=P
 1237  CONTINUE
      ENDIF
 1236 CONTINUE
 1235 CONTINUE

      if(bcsp.eq.0)then
       ntp = 0
       do icp=1,NPO
        deg(icp) = 0.d0
        ntp = ntp + lj(icp) + 1
        deg(icp) = lj(icp) + 1
        if(ntp.ge.int(jz))go to 1242
        end do
        stop '>>> ERROR IN REATTRIBUTING THE OCCUPANCIES'
 1242   if(ntp.gt.int(jz))deg(icp)=deg(icp)-(float(ntp)-jz)
         do iextra=icp+1,NPO
          deg(iextra) = 0.d0
         end do
      end if

      DO I=NPO+1,NOS
       DO J=I+1,NOS
        IF(EHF(J).LE.EHF(I))THEN
         P=EHF(I)
         EHF(I)=EHF(J)
         EHF(J)=P
         NTRY=NN(I)
         NN(I)=NN(J)
         NN(J)=NTRY
         LTRY=LL(I)
         LL(I)=LL(J)
         LL(J)=LTRY
         JTRY=LJ(I)
         LJ(I)=LJ(J)
         LJ(J)=JTRY
         if(bcsn.eq.1)then
          P=DEG(I)
          DEG(I)=DEG(J)
          DEG(J)=P
          P=FO(I)
          FO(I)=FO(J)
          FO(J)=P
         end if
         EVTRY=EVSO(I)
         EVSO(I)=EVSO(J)
         EVSO(J)=EVTRY
         STRY=ST(I)
         ST(I)=ST(J)
         ST(J)=STRY
         p=deltan(i)
         deltan(i)=deltan(j)
         deltan(j)=p
         p=dpren(i)
         dpren(i)=dpren(j)
         dpren(j)=p
         DO J1=1,NNN
          P=UNL(J1,I)
          UNL(J1,I)=UNL(J1,J)
          UNL(J1,J)=P
          P=DUNL(J1,I)
          DUNL(J1,I)=DUNL(J1,J)
          DUNL(J1,J)=P
          P=AVP(J1,I)
          AVP(J1,I)=AVP(J1,J)
          AVP(J1,J)=P
         end do
        ENDIF
       end do
      end do

      if(bcsn.eq.0)then
      ntn = 0
      do icn=NPO+1,nos
      deg(icn) = 0.d0
      ntn = ntn + lj(icn) + 1
      deg(icn) = lj(icn) + 1
      if(ntn.ge.int(ma-jz))go to 1243
      end do
      stop '>>> ERROR IN REATTRIBUTING THE OCCUPANCIES'
 1243 if(ntn.gt.int(ma-jz))deg(icn)=deg(icn)-(float(ntn)-ma+jz)
      do iextra=icn+1,nos
      deg(iextra) = 0.d0
      end do
      end if

      z_ver = 0.
      do i=1,NPO
       z_ver = z_ver + deg(i)
      end do

      amz_ver = 0.
      do i=NPO+1,nos
       amz_ver = amz_ver + deg(i)
      end do
    
      if(dabs(z_ver-jz).gt.1e-3.
     1   or.dabs((z_ver+amz_ver)-ma).gt.1e-3)then
       write(2,*)
       write(2,*) '>>> WRONG COUNT OF PARTICLES: '
       write(2,*) jz,z_ver,dabs(z_ver-jz)
       write(2,*) ma,z_ver+amz_ver,dabs((z_ver+amz_ver)-ma)
       write(2,*)
       stop
      end if

      do i=1,nos
       fo(i)=deg(i)/(lj(i)+1)
      enddo 

      do I=1,NPO
       do k=1,NPO
        GGG=0.d0
        do lcount=0,13
         MULTIP=lcount 
         elrid1=yl(LL(I),LJ(I),LL(K),LJ(K),MULTIP)
         GGG=GGG+elrid1**2
        enddo
        elmar(I,K)=GGG
       enddo
      enddo
C     parte di neutroni
      do I=NPO+1,NOS
       do k=NPO+1,NOS
        GGG=0.d0
        do lcount=0,13
         MULTIP=lcount
         elrid1=yl(LL(I),LJ(I),LL(K),LJ(K),MULTIP)
         GGG=GGG+elrid1**2
        enddo
        elmar(I,K)=GGG
       enddo
      enddo

 1241 continue

c_order_levels_END


CCC                                                                             
CCC   DENSITiIES                                                                  
CCC
      check_val_int = 0.d0
      partial_dp=0.d0
      partial_dn=0.d0   

      DO  J=1,NNN                    
        X = RM(J)                          
        R2 = RM2(J)                        
        R3 = RM3(J)                        
        E = 0.      
        F = 0.                             
        E2= 0.                             
        F2 = 0.                           
        TE = 0.                            
        TF= 0.

        DO I= 1,NOS                 
          L = LL(I)                      
          N = NN(I)                     
          JM = LJ(I)                                                 
          LL1 = L*(L+1)                                              
          FL = (JM*(JM+2)-4 *L*(L+1) - 3.)/8.                        
          EVSO (I) = FL                                             
          Y = UNL(J,I)*UNL(J,I)*DEG(I)                               
          E = E + Y*(1-LT(I))
          F = F + Y*LT(I)                                            
          EE = 2*FL*Y                                                
          E2 = E2 + EE*(1-LT(I))                                     
          F2 = F2 + EE*LT(I)                                         
          Y = (-UNL(J,I)/X + DUNL(J,I))/X                           
          Y1 = UNL(J,I)/R2                                          
          Y2 = DEG(I)*(Y*Y + LL1*Y1*Y1)                                
          TE = TE + Y2*(1-LT(I))                                       
          TF = TF + Y2*LT(I)                                         
        end do

        DP(J) = F*UNS4PI/R2                                       
        DN(J) = E*UNS4PI/R2                                     
        DT(J) = DN(J) + DP(J)       
        D_IV(J) = DN(J)-DP(J) 
        BETA(J) = (DN(J)-DP(J))/DT(J)      ! asymmetry   

        DJN(J) = E2*UNS4PI/R3                                  
        DJP(J) = F2*UNS4PI/R3                                 
        DJT(J) = DJN(J) + DJP(J)
        DJ_IV(J)= DJN(J) - DJP(J)
     
        TAUN(J) = TE*UNS4PI                                 
        TAUP(J) = TF*UNS4PI                              
        TAUT(J) = TAUN(J) + TAUP(J)     
        TAU_IV(J) = TAUN(J) - TAUP(J)  

      end do  

CCC                                                                             
CCC   Density derivatives                                                  
CCC                                                                             
      DO J=1,NNN                                                    
        FCT(J) = DN(J)    ! load into fct and get its der. by calling  DERIV
      end do                                                  
      CALL DERIV                                                     
      DO J=1,NNN                                                
        DN1(J) = DF(J)                                               
        FCT(J) = DP(J)
      end do                                              
      CALL DERIV                                                 
      DO J=1,NNN                                             
        DP1(J) = DF(J)                                           
        DT1(J) = DN1(J) + DP1(J)   
        D_IV1(J) = DN1(J) - DP1(J)    
        FCT(J) = DN1(J)
      end do                                        
      CALL DERIV                                            
      DO J=1,NNN                                        
        DN2(J) = DF(J)                                      
        FCT(J) = DP1(J)
      end do                                   
      CALL DERIV                                        
      DO J=1,NNN                                   
        DP2(J) = DF(J)                                  
        DT2(J) = DP2(J) + DN2(J)              
        D_IV2(J) = DN2(J) - DP2(J)
        FCT(J) = DJN(J)
      end do                               
      CALL DERIV                                   
      DO J=1,NNN                              
        DJN1(J) = DF(J)                            
        FCT(J) = DJP(J)
      end do      
      CALL DERIV                                                    
      DO J=1,NNN                                               
        DJP1(J) = DF(J)                                             
        DJT1(J) = DJP1(J) + DJN1(J)
        DJ_IV1(J) = DJN1(J) - DJP1(J)
      end do 

CCC                                                                             
CCC   CALCULATING FIELDS                                                     
CCC
CCC   COULOMB TERM                                                          
CCC
      E = 0.                                                    
      F = 0.                                                  
      DO J=1,NNN                                          
        X = RM(J)                                              
        Y = QP*X*X                                            
        E = E + Y*DP(J)                                      
        F = F + Y*DP(J)/X                                          
        VCD(J) = E/X - F
      end do       

      E222=((12./QP)**0.333333)*E22/H                                  
      DO J=1,NNN                                                  
        VCD(J) = E22*(VCD(J)+F)
        VC(J)=VCD(J)-float(icoul_e)*E222*(DP(J)**0.333333)
      end do                               
C_GL_11_4_05
      DO J=1,NNN              
        FCT(J) = VC(J)       
      end do      
      CALL DERIV                 
      DO J=1,NNN           
        DVC(J) = DF(J)
      end do


c     Starts defining SKYRME MEAN FIELD         
c     STARTING CYCLE
      DO 167  J = 1,NNN               

cccc  eq. 21 Vautherin-Brink:   HB: bare; HN, HP: con massa efficace

      HN=HB+CRO*DT(J)+CROQ*DN(J)+C13*(DT(J)**2-DN(J)**2)               
      HP=HB+CRO*DT(J)+CROQ*DP(J)+C13*(DT(J)**2-DP(J)**2)


ccc   For now, BARE mass
      if (FRANCESCO) then
          HN = HB
          HP = HB
      end if 


c     update masses
      HMEN(J)=HN*(1.-XMU)+HMEN(J)*XMU                                  
      HMEP(J)=HP*(1.-XMU)+HMEP(J)*XMU                                  
      MESMN(J) = HB/HMEN(J)                                            
      MESMP(J) = HB/HMEP(J)


      HN1=CRO*DT1(J)+CROQ*DN1(J)+C13*2.*(DT(J)*DT1(J)-DN(J)*DN1(J))    
      HP1=CRO*DT1(J)+CROQ*DP1(J)+C13*2.*(DT(J)*DT1(J)-DP(J)*DP1(J))   
      HN2=CRO*DT2(J)+CROQ*DN2(J)+C13*2.*(DT1(J)**2+DT(J)*DT2(J)        
     1-DN1(J)**2-DN(J)*DN2(J))                                         
      HP2=CRO*DT2(J)+CROQ*DP2(J)+C13*2.*(DT1(J)**2+DT(J)*DT2(J)        
     1-DP1(J)**2-DP(J)*DP2(J))


ccc  Bare mass => derivatives are set to zero
      if (FRANCESCO) then
             HN1 = 0.
             HP1 = 0.
             HN2 = 0.
             HP2 = 0.
      end if 

      DHMEN(J)=HN1*(1.-XMU)+DHMEN(J)*XMU                               
      DHMEP(J)=HP1*(1.-XMU)+DHMEP(J)*XMU                               
      D2HMEN(J)=HN2*(1.-XMU)+D2HMEN(J)*XMU                             
      D2HMEP(J)=HP2*(1.-XMU)+D2HMEP(J)*XMU 


cc    Defining the potential  
  167 CONTINUE                                                         
      DO 1500  J=1,NNN                                                 
        VNT33=(1.-X3)*(ALFE+2.)*DT(J)*DT(J)                              
        VNT3=(VNT33+(2.+4.*X3)*DP(J)*(DT(J)+ALFE*DN(J)))/24.             
        VPT3=(VNT33+(2.+4.*X3)*DN(J)*(DT(J)+ALFE*DP(J)))/24.             
      IF(ABS(ALFE-1.).LT.1.E-06) GO TO 8888                            
      if (dt(j).eq.0) then  
        vnt=0.
      else 
        VNT=DT(J)**(ALFE-1.)                                             
      endif

      VNT3=VNT3*VNT                                                    
      VPT3=VPT3*VNT                                                    
 8888 CONTINUE                                                         
      VNT3=-VNT3                                                       
      VPT3=-VPT3                          
                             
      VN(J) = TXRO*DT(J) + TXDRO *DT1(J)/RM(J) + TXD2RO*DT2(J) +       
     1 TXTAU*TAUT(J)-W3*(DJT(J)/RM(J)+DJT1(J)/2.)                       

      if(iter.eq.1.and.j.eq.1)then
        write(2,*)
        write(2,*) isj,'  beta_C = ',-.125*(T1*X1 + T2*X2)
        write(2,*) itens,'  beta_T = ',5.d0*(stc_t+stc_u)/24.d0
        write(2,*)
      end if


      VN(J)=VN(J)-C13*(2.5*DT(J)*(DT2(J)+2.*DT1(J)/RM(J))              
     1+1.25*(DT1(J)**2)-2.*DT(J)*TAUT(J))       
     
      
      if(iter.eq.1)then
        write(96,*) float(j)*h,W0* DT1(J)/(2.* RM(J)), 
     1  -.125*(T1*X1 + T2*X2)*DJT(J)/RM(J),
     1  5.d0*(stc_t+stc_u)/24.d0*DJT(J)/RM(J),w0
      end if    
   
      VP(J) = VN(J) + TXROQ*DP(J) + TXDROQ* DP1(J)/RM(J) +TXD2RQ *    
     1  DP2(J) +TXTAUQ *TAUP(J) -W3P*(DJP(J)/RM(J) +.5*DJP1(J)) - .25*  
     2  T3*4.*VPT3    

      if(iter.eq.1)then
       write(96,*) float(j)*h,VN(J),TXROQ*DP(J),TXDROQ* DP1(J)/RM(J), 
     1 TXD2RQ*DP2(J),TXTAUQ *TAUP(J),-W3P*(DJP(J)/RM(J) +.5*DJP1(J))
     1 ,- .25*T3*4.*VPT3,vp(j)
      end if

      VP(J)=VP(J)-C13*(-2.5*DP(J)*(DP2(J)+2.*DP1(J)/RM(J))-1.25*(DP1(J)
     1 **2)+2.*DP(J)*TAUP(J)-0.25*(DJP(J)**2)-0.50*(DJN(J)**2))         
 

      if(iter.eq.1.and.j.eq.1)then
       write(2,*)
       write(2,*) isj,'  alpha_C = ',.125*(T1-T2)-.125*(T1*X1 + T2*X2)
       write(2,*) itens,'  alpha_T = ',5.d0*stc_u/12.d0
       write(2,*)
      end if


cc     common block of the s.o. potential (total densities)
       VSON(J) = W0* DT1(J)/(2.* RM(J)) 
     1-2.d0*float(isj)*.125*(T1*X1 + T2*X2)*DJT(J)/RM(J)      
     1+2.d0*float(itens)*5.d0*(stc_t+stc_u)/24.d0*DJT(J)/RM(J)  
    

ccc   In input there is f lag to turn SO on/off
      VSOP(J) = VSON(J) + W0P*DP1(J)/(2.*RM(J))
     1+2.d0*float(isj)*.125*(T1-T2)*DJP(J)/RM(J) 
     1+2.d0*float(itens)*5.d0*(stc_u-stc_t)/24.d0*DJP(J)/RM(J)


ccc   Change here PROTONS    
      if (FRANCESCO) then
          tauz = -1.d0 
          VP(J) = mean_field(DT(J),BETA(J),tauz)
     1  + 2.d0*c0deltarho*(dt2(j)  +2.d0*dt1(j)/rm(j))
     1  + 2.d0*c1deltarho*(d_iv2(j)+2.d0*d_iv1(j)/rm(j))*tauz
ccc       Spin-orbit (nabla j + nabla j_p)
     1  - w0_0/2.* (  (djt1(j) + 2.d0/rm(j)*djt(j)) + 
     1    (djp1(j) +2.d0/rm(j)*djp(j))  ) 
c          VSON(J) = W0* DT1(J)/(2.* RM(J))
c          VSOP(J) = VSON(J) + W0P*DP1(J)/(2.*RM(J))
      end if 
      

      if(iter.eq.1)then
       write(96,*) float(j)*h,vson(j),vsop(j)
      end if   
    
 
ccc   I: loop over protons;   J: position (external loop)
ccc   Update average potential PROTONS
      do  I=1,NPO                                  
      AVP(J,I) = (VP(J) + 
     1            float(iso)*VSOP(J)*EVSO(I) + 
     1            float(icoul_d)*VC(J) +
     1            i_ext*harmonic_potential(rm(j),omega) )        ! external
     1          *(1. -XMU) + XMU*AVP(J,I)
      end do




ccc   NEUTRONS      
      VN(J) = VN(J) + TXROQ*DN(J) + TXDROQ*DN1(J)/RM(J) + TXD2RQ* DN2(J
     1 ) + TXTAUQ * TAUN(J) -W3P*( DJN(J)/RM(J) + .5*DJN1(J)) -.25*T3   
     2 *4.*VNT3

      VN(J)=VN(J)-C13*(-2.5*DN(J)*(DN2(J)+2.*DN1(J)/RM(J))-1.25*(DN1(J)
     1 **2)+2.*DN(J)*TAUN(J)-0.25*(DJN(J)**2)-0.50*(DJP(J)**2))      
 
      VSON(J) = VSON(J) + W0P*DN1(J)/(2.*RM(J))
     1 +2.d0*float(isj)*.125*(T1-T2)*DJN(J)/RM(J)
     1 +2.d0*float(itens)*5.d0*(stc_u-stc_t)/24.d0* DJN(J)/RM(J)
     

      if (FRANCESCO) then
ccc   Define VN(J) here     (att: in VN +1-beta;  in VP, -1-beta)
        tauz=+1.d0
        VN(J) = mean_field(DT(J),BETA(J),tauz)
cc    Gradients
     1  + 2.d0*c0deltarho*(dt2(j)  +2.d0*dt1(j)/rm(j))
     1  + 2.d0*c1deltarho*(d_iv2(j)+2.d0*d_iv1(j)/rm(j))*tauz
ccc   Spin-orbit (nabla j + nabla j_n)
     1  - w0_0/2.d0 * (  (djt1(j) + 2.d0/rm(j)*djt(j)) +
     1    (djn1(j) +2.d0/rm(j)*djn(j))  )
cc        VSON(J) = W0* DT1(J)/(2.* RM(J))
cc        VSON(J) = VSON(J) + W0P*DN1(J)/(2.*RM(J))
      end if 



      NNEUT  = NOS - NPO

cc    Update average potential NEUTRON
      DO I2=1,NNEUT                                              
        I = I2 + NPO                                                    
        AVP(J,I) = ( VN(J) 
     1         + float(iso)*VSON(J)*EVSO(I) 
     1         + i_ext*harmonic_potential(rm(j),omega) )   ! external       
     1   *(1.-XMU) + XMU*AVP(J,I)

        if(i.eq.13)then
         write(85,*) dfloat(j)*h,
     1   W0P*DN1(J)/(2.*RM(J))
         write(85,*) dfloat(j)*h,vn(j),vson(j)
        end if
      end do


 1500 CONTINUE
      diffmax=0.d0
      do i=1,nos
       if(iter.eq.1)go to 1505
       diff=dabs(ehf(i)-ehf_old(i))
       if(diff.gt.diffmax) diffmax=diff
 1505  ehf_old(i)=ehf(i)
      end do

c     Convergence    
      if(iter.gt.1.and.diffmax.lt.1e-3)go to 116  ! def. value: 1e-5
      IF(iter. GE.ITMAX)go to 1506
      GO TO 104     

CCC                                                                             
CCC   END of iteration
CCC



ccc   Write stuff to file
 1506 write(2,*)
      write(2,*) '>>> convergence reached only at the level: '
      write(2,*) diffmax
      if (debug) then
          print *, "Max. difference = ", diffmax
      end if
      stop '>>> convergence not reached'
      go to 1507

  116 continue
      write(2,*) 'reached desired convergence: '
      write(2,*) diffmax
 1507 write(2,*)
c     stop '>>> max number of iterations reached' 
      do i=1,NPO
       deltag(i)=deltap(i)
      enddo
      do i=NPO+1,nos
       deltag(i)=deltan(i)
      enddo


C_GL_11_4_05
      fac_em_so = 0.022075d0

      if(ichf.eq.0.and.ibox.eq.1)
     &open(unit=99,status='unknown',file='lines')
      do i=1,nos
       fac_mu = 2.39d0
       if(i.gt.NPO)fac_mu = -1.91d0
       fac_l = -ll(i)-1
       ll2=2*ll(i)
       if(lj(i).gt.ll2)fac_l = ll(i)
       aint=0.d0
       write(7,*) ehf(i),nn(i),ll(i),lj(i),lt(i)
       do iprint=1,nnn2
       r=float(iprint)*h
       write(7,*) r,unl(iprint,i)
       aint=aint+h*unl(iprint,i)**2/r*dvc(iprint)
       end do
       aint = fac_em_so*fac_mu*aint*fac_l
      write(2,2221) st(i),lt(i),ehf(i),deg(i),deltag(i),fo(i),aint    
      if(ichf.eq.0.and.ibox.eq.1)then
       if(fo(i).ge.0.00001d0) then 
        write(99,*) nn(i),ll(i),lj(i),lt(i)
       end if
      end if
      if(ichf.eq.0.and.i_grnp0_box.eq.1)then
       if(fo(i).ge.0.00001d0) then 
        write(50,*) ehf(i),nn(i),ll(i),lj(i),lt(i)
        write(50,*)(unl(j,i),j=1,nnn2)
        write(50,*)(dunl(j,i),j=1,nnn2)
       end if
      end if
CCC   <r^2> for each orbit
c     RO=0.d0                                                                     
c     DO J=1,NNN                                                            
c     X2=RM2(J)
c     RO=RO+X2*UNL(J,I)**2                                                         
c     END DO                                                                    
c     RO=RO*H
      enddo

      if(ichf.eq.0.and.ibox.eq.1)then
       close(99)
       call create_list_unocc
      end if

      emaxp = -1000.d0
      do 1160 i=1,NPO
      if(ehf(i).gt.emaxp)emaxp=ehf(i)
 1160 continue      

      emaxn = -1000.d0
      do 1161 i=NPO+1,nos
      if(ehf(i).gt.emaxn)emaxn=ehf(i)
 1161 continue      

c     write(2,*)
c     write(2,*) 'emaxp,emaxn= ',emaxp,emaxn
c     write(2,*) 


cc    doing something with Coulomb
      ecbd=0.                                                          
      ecbe=0.                                                          
      DO J=1,NNN                                                  
        d_coul_dens(j) = vcd(j)*dp(j)/2.d0 
        cost_ = ((3.d0/pi)**(1./3.)) * qe2 
        e_coul_dens(j) = -0.75 * (dp(j)**(4./3.)) * cost_

        ecbd = ecbd + QP*RM2(J)*h * VCD(J)*DP(J)*0.5     
        ecbe = ecbe + QP*RM2(J)*h * (DP(J)**(4./3.))                 
      end do
      costalpha=E22*((12./QP)**(1./3.))/H
c     write(23,*) 'C,ecbe=',costalpha,ecbe
      ecbe=ecbe*costalpha
      ecbe= -ecbe*0.75      
      if(icoul_d.eq.0)ECBD=0.d0
      if(icoul_e.eq.0)ECBE=0.d0
      WRITE(2,5556) ecbd, ecbe                                      
 5556 FORMAT(2X,'ECB(DIRECT)=',E12.5,5X,'ECB(EXCHANGE)=',E12.5)       
      IF(IBOX.NE.1) GO TO 461                                         
  453 FORMAT(2X,'N=',I4,5X,'NOT CONVERGED')                           
      nmax = nmaxt
      eps = 1.e-6
      em = 50.d0
c     read(14,*) NMAX,EPS,EM                                         
      EM1=EM*DMSHB                                                  
      open(unit=99,status='old',file='lines')
      do i=1,nos
      read(99,*) i_dum1,i_dum2,i_dum3,i_dum4
      end do
  460 read(99,*,end=461) L,J,NT,IND,INF                                
      INN=INF-IND+1                                              
      N=IND-1                                                   
      JK=0                                                     
  463 N=N+1                                                  
      JK=JK+1                                               
      IF(N.GT.INF) GO TO 4460                              
      LL1=L*(L+1)                                         
      K=0                                                
      S2=0.                                             
      ESO=(J*(J+2)-4*L*(L+1)-3.)/8.                    
      NMAX1=NMAX+1                                    
      DO464 KK=1,NMAX1                                                 
      IF(KK-NNN) 465,465,466                                        
  465 IF(NT.EQ.0) then 
      VQ(KK)=MESMN(KK)*(VN(KK)+float(iso)*ESO*VSON(KK)
     1-0.25*DHMEN(KK)*DHMEN(KK)/HMEN(KK)+0.5*D2HMEN(KK))
      else IF(NT.EQ.1) then
      VQ(KK)=MESMP(KK)*(VP(KK)+float(icoul_d)*VC(KK)
     1+float(iso)*ESO*VSOP(KK)-0.25*   
     1DHMEP(KK)*DHMEP(KK)/HMEP(KK)+0.5*D2HMEP(KK))       
      end if       
      VQ(KK)=VQ(KK)/HB+LL1/(RM(KK)*RM(KK))                             
      RMSM(KK)=(1.-MESMN(KK))*(1.-NT)+(1.-MESMP(KK))*NT                
      GO TO 464                                                        
  466 VQ(KK)=LL1/(H*H*FLOAT(KK*KK))                                    
      RMSM(KK)=0.                                                      
  464 CONTINUE                                                         
      GO TO 480                                                        
  470 WRITE(2,453) N                                                 
      GO TO 463                                                        
  480 NO=N-1                                                           
      ideb=0
      EMIN0=-500.
      emax0=500.
      ecut=500.
      niter2=420
      emin=emin0
      nprim=0
      eprim=0.
      IF (NPRIM.EQ.(N-1)) EMIN=EPRIM

      nprim=n
      eprim=e

      ideb=1
      CALL NUM2B(NMAX,H,E,S2,VQ,RMSM,UQ,NO,HB)                          
      IF(NO.LT.0) GO TO 470                                            
      E=E/DMSHB
      NOL=NO+1                                                         

c     write(6,*) NMAX,EPS,EM

      DO  KX=1,NMAX                                                
       VQ(KX)=VQ(KX)/DMSHB                                              
       WF(NOL,KX)=UQ(KX)                                                
      end do
c     WRITE (*,454) NOL,E                                             
c     WRITE (*,211) (UQ(IK),IK=1,NMAX)                                

      if(i_grnp0_box.eq.1)then
c     write(7,*) nol,l,j,nt,e  
      write(51,*) e,nol,l,j,nt 
      write(51,510)(UQ(IK),IK=1,nmax-2)
  510 format(6(e12.5,1x)) 
      write(7,*) e,nol,l,j,nt
      do iprint=1,nmax-2
      r=float(iprint)*h
      write(7,*) r,uq(iprint)
      end do
      DO ik=1,nmax                                                     
       FCT(ik) = UQ(ik)
      end do 
      CALL DERIV                                                       
      write(51,510)(df(ik),ik=1,nmax-2)
      end if

9998  FORMAT(6E12.5)                                                   
      GO TO 463                                                        
 4460 continue
c     WRITE (*,456)                                                   
      write(2,*)
      write(2,*) 'Overlap matrix for L,J,NT: ',l,j,nt
      DO4464 K1=IND,INF                                                
      DO4462 K2=K1,INF                                                 
      OV(K1,K2)=0.                                                     
      DO 4463 KX=1,NMAX                                                
      OV(K1,K2)=OV(K1,K2)+WF(K1,KX)*WF(K2,KX)*H                        
 4463 continue
      OV(K2,K1)=OV(K1,K2)                                              
 4462 CONTINUE                                                         
      WRITE (2,211) (OV(K1,K2),K2=IND,INF)           
                
 4464 CONTINUE                                                         
      GO TO 460                                                        
  461 CONTINUE                                                         
      N64=64                                                          
 1002 FORMAT(A6,4I5,E12.5)                                             
      MASH(1)=12                                                       
      MASH(2)=24                                                       
      MASH(3)=36                                                       
      NX1=2                                                            
      NX2=2                                                            
      NX3=2                                                           
      DAL(1)=0.20                                                      
      DAL(2)=0.20                                                      
      DAL(3)=0.20                                                      
      if(i_grnp0_box.eq.1)then
      write(10,*) nnn2
      write(10,9997) T0,T1,T2,T3,T13,X0                                
      write(10,9997) ALFE,X3,X1,X2,W2,W2P                              
      FJZ=FLOAT(JZ)                                                    
      write(10,9997) MA,FJZ,H,(DAL(I),I=1,3)                           
      write(10,9996) (MASH(I),I=1,3)                                   
      NXXX=1                                                           
      ND=1                                                             
      write(10,9996) NXXX,ND,NX1,NX2,NX3                              
      write(10,9998) (HMEN(J),J=1,NNN2),(HMEP(J),J=1,NNN2)             
      write(10,9998) (VN(J),J=1,NNN2),(VP(J),J=1,NNN2)                 
      write(10,9998) (DHMEN(J),J=1,NNN2),(DHMEP(J),J=1,NNN2)           
      write(10,9998) (VSON(J),J=1,NNN2),(VSOP(J),J=1,NNN2)             
      write(10,9998) (VC(J),J=1,NNN2)                                  
      write(10,9998) (DT(J),J=1,NNN2),(TAUT(J),J=1,NNN2)               
      write(10,9998) (DT1(J),J=1,NNN2),(DT2(J),J=1,NNN2)               
      write(10,9998) (DP(J),J=1,NNN2),(DN(J),J=1,NNN2)                 
      go to 1791
      end if

      IF(IPCHPO.EQ.0) GO TO 1789                                       
      write(12,*) nnn
      write(12,9997) T0,T1,T2,T3,T13,X0                                
      write(12,9997) ALFE,X3,X1,X2,W2                                  
      FJZ=FLOAT(JZ)                                                    
      write(12,9997) MA,FJZ,H,(DAL(I),I=1,3)                           
      write(12,9996) (MASH(I),I=1,3)                                   
      NXXX=1                                                           
      ND=1                                                             
      write(12,9996) NXXX,ND,NX1,NX2,NX3                               
      write(12,9998) (HMEN(J),J=1,NNN),(HMEP(J),J=1,NNN)               
      write(12,9998) (VN(J),J=1,NNN),(VP(J),J=1,NNN)                   
      write(12,9998) (DHMEN(J),J=1,NNN),(DHMEP(J),J=1,NNN)             
      write(12,9998) (VSON(J),J=1,NNN),(VSOP(J),J=1,NNN)               
      write(12,9998) (VC(J),J=1,NNN)                                   
      write(12,9998) (DT(J),J=1,NNN),(TAUT(J),J=1,NNN)                 
      write(12,9998) (DT1(J),J=1,NNN),(DT2(J),J=1,NNN)                 
      write(12,9998) (DP(J),J=1,NNN),(DN(J),J=1,NNN)            

 1791 continue


c     Writing densities
      do J=1,NNN
        write(16,*) H*J,DP(J),DN(J)
        write(56,*) H*J,DP(J)/FJZ,DN(J)/(MA-FJZ)
c       write(90,*) h*j,dn_valence(j),dn(j)-dn_valence(j)
        write(28,280) H*J,DP1(J),DN1(J),DP1(J)+DN1(J)
        write(29,*) h*j,djp(j),djn(j) 
      enddo
  280 format(4(2x,e12.5))


c density for the FT

      rhomin_ft = dlog(0.0001d0)
      rhomax_ft = dlog(30.d0)
      akmax_ft = 500.d0
      nexp_ft = 8
      nr_ft = 2**nexp_ft
      dr_ft = (rhomax_ft-rhomin_ft)/dfloat(nr_ft-1)
      akpmin_ft = dlog(akmax_ft)-rhomax_ft+rhomin_ft

      rr_ft(1) = dexp(rhomin_ft)
      cf_ft = dexp(dr_ft)
      do i = 2,nr_ft
       rr_ft(i) = cf_ft*rr_ft(i-1)
      enddo

      do i=1,nr_ft
       ind_ft=(rr_ft(i)-h)/h+1.0001 
       if(ind_ft.gt.nnn) then
        dens1_ft(i)=0.d0
        dens2_ft(i)=0.d0
       end if
      end do

cc       end if
      resto_ft=(rm(ind+1)-rr_ft(i))/h 
      dens1_ft(i)=dp(ind)+dp(ind+1)
      dens1_ft(i)=dens1_ft(i)+(dp(ind+1)-dp(ind))*
     &         (1.d0-2.d0*resto_ft)
      dens1_ft(i)=dens1_ft(i)/2.d0
      dens2_ft(i)=dn(ind)+dn(ind+1)
      dens2_ft(i)=dens2_ft(i)+(dn(ind+1)-dn(ind))*
     &         (1.d0-2.d0*resto_ft)
      dens2_ft(i)=dens2_ft(i)/2.d0
400   continue
      write(29,*) rr_ft(i),dens1_ft(i),dens2_ft(i)
3     continue


c     Writing potentials
      write(17,*) 'r     V(p)     V(n)     m*(p)/m     m*(n)/m    VC '
      do J=1,NNN
        write(17,1792) H*J,VP(J),VN(J),VSOP(J),VSON(J),VC(J)  !,VSOP(J),VSON(J)

      enddo
1792  format(6(2x,e12.5))

c     Writing occupation factors
      if(bcsp.eq.0)lamp=0.d0
      if(bcsn.eq.0)lamn=0.d0
      write(27,*) lamp,lamn,bcsp,bcsn 
      do i=1,nos
      write(27,*) fo(i)
      enddo


 9997 FORMAT(6E12.5)                                                  
 9996 FORMAT(15I3)
      do j=1,nnn
         write(11,*) dfloat(j)*h,djp(j),djn(j)
      end do 
      espo = 0.d0
      aint1 = 0.d0
      aint2 = 0.d0
ccc   spin-orbit energy
      do j=1,nnn
        espo = espo + w2 *dt(j)*(2.d0*rm(j)*djt(j)+rm2(j)*djt1(j)) 
     2              + w2p*dp(j)*(2.d0*rm(j)*djp(j)+rm2(j)*djp1(j)) 
     2              + w2p*dn(j)*(2.d0*rm(j)*djn(j)+rm2(j)*djn1(j) )
ccc     spin-orbit  density
        so_dens(j) = 
     1  (  w2*dt(j)*( 2.d0*djt(j)/rm(j) + djt1(j) ) 
     1  + w2p*dp(j)*( 2.d0*djp(j)/rm(j) + djp1(j) ) 
     1  + w2p*dn(j)*( 2.d0*djn(j)/rm(j) + djn1(j) ) ) *(-0.5d0)


        aint1 = aint1 + h*qp*rm2(j)*(djp(j)**2+djn(j)**2)
        aint2 = aint2 + h*qp*rm2(j)*djp(j)*djn(j) 
      end do
      espo = espo*qp*(-0.5d0)*h

      prefac = -4.d0 / 3.d0 / (ma-2*jz) 
      fac1 = 5.d0/24.d0 * (stc_t+3.d0*stc_u)
      fac2 = 5.d0/12.d0 * stc_u 
      e_gt_shift = prefac*(aint2*fac1+aint1*fac2) 
      prefac_old = 1.d0 / 9.d0 / (ma-2*jz) 
      facold1 = (stc_t - 9.d0*stc_u ) 
      facold2 = -5.d0*stc_u 
      e_old = prefac_old*(aint2*facold1+aint1*facold2)
c     e_gt_shift = -10.d0*qp/36.d0/(ma-jz) * e_gt_shift 
      
      write(2,*)
      write(2,*) 'prefac,fac1,fac2,aint1,aint2: ',prefac,fac1,fac2
     1,aint1,aint2
      write(2,*) 'e_gt_shift: ',e_gt_shift
      write(2,*) 'old WRONG result was: ',e_old
c0995x 
      if(icount_chf.eq.1)then
      write(2,*) 
      write(2,*) 'E(s.o.) = ',espo 
      end if


 1789 IF(IPCHWF.EQ.0) GO TO 1790                                       
      DO9999 I=1,NOS                                                   
      IF(DEG(I).LT.0.1) GO TO 9999                                     
      write(12,1002) ST(I),NN(I),LL(I),LJ(I),LT(I),EHF(I)              
      I1I=ND-1+NX1                                                     
      I1F=MASH(1)*NX1                                                  
      I2I=I1F+NX2                                                      
      I2F=I1F+(MASH(2)-MASH(1))*NX2                                    
      I3I=I2F+NX3                                                      
      I3F=I2F+(MASH(3)-MASH(2))*NX3                                    
c     if(i.eq.11)then
                                       
 9999 CONTINUE                                                         
 1790 CONTINUE                                                         


ccc   scrittura dei delta di ogni stato per input qrpa
      open(unit=20,status='unknown',file='delta.dat')
      do i=1,NPO
      write(20,*) deltap(i)
      enddo
      do i=NPO+1,nos
      write(20,*) deltan(i)
      enddo
      close(20)





CCC
CCC   Calculating energies                                                       
CCC

ccc   Initialize several quantities to zero
      erea =0.                                                          
      e_di_t0=0.d0
      e_di_constr=0.d0
      e_di_t12=0.d0
      e_di_t3=0.d0
c      e_di_so=0.d0
c     Custom potential
      e_bulk = 0.d0
      e_ext = 0.d0
C_GL_19_04
      e_d_coul=0.d0


cc    debug   Doba     

      e_rho0 =0.d0
      e_rho1 =0.d0
      e_rho  =0.d0

      e_deltarho0=0.d0
      e_deltarho1=0.d0
      e_deltarho =0.d0

      e_tau0=0.d0
      e_tau1=0.d0
      e_tau =0.d0

      e_nablaj0=0.d0
      e_nablaj1=0.d0
      e_nablaj =0.d0

cc    tensor, J2
      e_t0 =0.d0
      e_t1 =0.d0
      e_t  =0.d0

      e_bulk0= 0.d0
      e_bulk1= 0.d0
      e_bulk = 0.d0

      e_surf0= 0.d0
      e_surf1= 0.d0
      e_surf = 0.d0


ccc   Energy density
cc    R, energy density, kin., pot., central, so, coul_d, coul_e
cc    Terms: kin, pot, spin-orbit, Coul (d/e), ext       
      bulk_arr = 0.d0
      surf_dens=0.d0
      do j = 1,nnn
        kin_dens(j) = taut(j)*hb
!        bulk_arr(j) = 
!     1  float(icoul_d)*d_coul_dens(j) + float(icoul_e)*e_coul_dens(j) 
        surf_dens(j) = 
     3  so_dens(j)*float(iso) + 
     3  i_ext *harmonic_potential(rm(J),omega) *dt(j)
        if (FRANCESCO) then
          bulk_arr(j) = bulk_arr(j) + pot_density(dt(j),beta(j))
          surf_dens(j) = surf_dens(j)
     2    - c0deltarho* dt1(j)**2 - c1deltarho* d_iv1(j)**2
c     2    c0deltarho*dt(j)*   (dt2(j)+dt1(j)/rm(j)) +
c     2    c1deltarho*d_iv(j)* (d_iv2(j)+d_iv1(j)/rm(j))
        else 
          cen_dens(j) =           ! central energy density (t0123)
     2  t0/4.*( (2.d0+x0)*dt(j)**2-(2.d0*x0+1.d0)*( dp(j)**2+dn(j)**2) )!t0
     2  + (1.d0/24.d0)*t3* dt(j)**alfe *                                !t3    
     2  ( (2.d0+x3)*dt(j)**2-(2.d0*x3+1.d0)*(dp(j)**2+dn(j)**2 ) )  
cc      H(eff)    effective mass contribution                           !t12
     2  + (1.d0/8.d0)*(t1*(2.d0+x1)+t2*(2.d0+x2))*taut(j)*dt(j) +
     2  (1.d0/8.d0)*(t2*(2.d0*x2+1.d0)-t1*(2.d0*x1+1.d0))*
     2  ( taup(j)*dp(j)+taun(j)*dn(j) ) +
cc      H(fin)    finite range (gradient, t1-t2) contribution
     2  (1.d0/32.d0)*( 3.d0*t1*(2.d0+x1)-t2*(2.d0+x2) )*dt1(j)**2 -
     2  (1.d0/32.d0)*( 3.d0*t1*(2.d0*x1+1.d0)+t2*(2.d0*x2+1.d0) )*
     2  ( dp1(j)**2+dn1(j)**2 )
          bulk_arr(j) = bulk_arr(j) + cen_dens(j)
        end if
        en_dens(j) = kin_dens(j) + bulk_arr(j) + surf_dens(j)
cc      todo qua è pieno di cose da sistemare
cc      R, energy density, kin. density, pot. d
        write(1997,'(8F10.4)'),rm(j),en_dens(j),kin_dens(j),bulk_arr(j)
cc      R, energy density, kin., pot., central, so, coul_d, coul_e
        write(1998,'(8F10.4)'), rm(j), en_dens(j), kin_dens(j),
     3   bulk_arr(j), cen_dens(j), so_dens(j),
     3  float(icoul_d)*d_coul_dens(j), float(icoul_e)*e_coul_dens(j)

      end do

cc    ectot    bare kinetic energy
      ectot=0.                                                         
      do J=1,NNN                                                  
        ectot=ectot + taut(J)*rm2(J)   
      end do                                            
      ectot=ectot*qp*h*hb          


c     Iterating over positions   -> computing integral
      do  J=1,NNN                           

c     Integrale potenziale "personalizzato" (in simm. sferica -> somma r^2 * integrando)
      if (FRANCESCO) then
        e_bulk = e_bulk + rm2(J) * pot_density(DT(J),BETA(J)) 
cc        e_surf = e_surf + rm2(j) * surf_dens(j)   
      end if


cc    c13 is usually 0
      ERE=-0.5*C13*1.25*(DP(J)*(DP1(J)**2)                             
     1+DN(J)*(DN1(J)**2)-DT(J)*(DT1(J)**2))                            
      ERE=ERE-0.5*C13*(TAUN(J)*(DP(J)**2)+TAUP(J)*(DN(J)**2)+2.*DT(J)  
     1*(TAUP(J)*DP(J)+TAUN(J)*DN(J)-TAUT(J)*DT(J)))                    
      ERE=ERE-C13*0.125*((DP(J)-2.*DT(J))*(DJP(J)**2)                  
     1+(DN(J)-2.*DT(J))*(DJN(J)**2))                                 
     


      if (FRANCESCO) then
        erea=erea+rm2(J)*rearrangement_density(DT(J),BETA(J))   
      else
        erea = erea-RM2(J)*ERE                                       
        erea = erea -RM2(J)*(0.5*(1.-X3)*DT(J)*DT(J)+(2.*X3+1.)*DN(J)*
     2  DP(J))*(DT(J)**ALFE)*T3*ALFE/24.
      end if




      e_di_t0=e_di_t0+0.25d0*t0*((2.d0+x0)*dt(j)**2-
     1(2.d0*x0+1.d0)*(dp(j)**2+dn(j)**2))*rm2(j)

      e_di_constr=e_di_constr+constr*dt(j)*rm2(j)*rm2(j)

      e_di_t3=e_di_t3+(1.d0/24.d0)*t3*dt(j)**alfe*
     1((2.d0+x3)*dt(j)**2-(2.d0*x3+1.d0)*(dp(j)**2+dn(j)**2))*
     1rm2(j)


ccc   E_ext(rho) = 4pi int( r^2 rho(r) v_ext(r) )      
      e_ext= e_ext+harmonic_potential(RM(J),omega)*dt(j)* rm2(j)   !external
cc      print *, "E(ext)   ", e_ext



      ctr=
cc   H(eff)    effective mass contribution
     1(1.d0/8.d0)*(t1*(2.d0+x1)+t2*(2.d0+x2))*taut(j)*dt(j)+
     1(1.d0/8.d0)*(t2*(2.d0*x2+1.d0)-t1*(2.d0*x1+1.d0))*
     1( taup(j)*dp(j)+taun(j)*dn(j) )+
cc   H(fin)    finite range (gradient, t1-t2) contribution
     1(1.d0/32.d0)*(3.d0*t1*(2.d0+x1)-t2*(2.d0+x2))*dt1(j)**2-
     1(1.d0/32.d0)*(3.d0*t1*(2.d0*x1+1.d0)+t2*(2.d0*x2+1.d0))*
     1( dp1(j)**2+dn1(j)**2 )

      e_di_t12 = e_di_t12+ctr*rm2(j)



cc    Displacement Coulomb energy 
      dint_cde=0.d0
      do jj=1,nnn
         if(jj.le.j)      fac=rm2(jj)/rm(j)
         if(jj.gt.j)      fac=rm(jj)
         dint_cde = dint_cde+(dn(jj)-dp(jj))*fac*h
      end do

      e_d_coul = e_d_coul + dint_cde*dp(j)*rm2(j)*h


cc    Doba energies
      if (debug) then

ccc       rho *lapl(rho) = - (grad rho)^2 = - (drho/dr)^2
          e_deltarho0 = e_deltarho0 - c0deltarho * DT1(J)**2 *rm2(J)
          e_deltarho1 = e_deltarho1 - c1deltarho *D_IV1(J)**2 *rm2(J)

ccc       rho * div(J)
          div = rm2(J)*DJT1(J)   + 2.*rm(J)*DJT(J)         ! div(J) * r^2
          e_nablaj0 =  e_nablaj0 + c0nablaj*DT(J)* div
          div = rm2(J)*DJ_IV1(J) + 2.*rm(J)*DJ_IV(J)
          e_nablaj1 =  e_nablaj1 + c1nablaj*D_IV(J)* div


        if (FRANCESCO.eqv..false.) then

ccc       rho^2 
          cc0 = c0rho + c0rho_dens * DT(J)**alfe
          cc1 = c1rho + c1rho_dens * DT(J)**alfe
          e_rho0 = e_rho0 + cc0 * DT(J)**2   * rm2(J)
          e_rho1 = e_rho1 + cc1 * D_IV(J)**2 * rm2(J)

ccc       rho*tau
          e_tau0 = e_tau0 + c0tau*DT(J)*TAUT(J)* rm2(J)
          e_tau1 = e_tau1 + c1tau*D_IV(J)*TAU_IV(J)* rm2(J)


ccc       J*J  tensor
          e_t0 = e_t0 + c0t*DJT(J)*  DJT(J)   * rm2(J)
          e_t1 = e_t1 + c1t*DJ_IV(J)*DJ_IV(J) * rm2(J)
         
         end if


        end if



      end do
ccc   Finished with the integral loop



ccc   Multiply by step (h) and 4pi (qp)
    
      e_di_t0=e_di_t0*h*qp
      e_di_constr=e_di_constr*h*qp
      e_di_t3=e_di_t3*h*qp
      e_di_t12=e_di_t12*h*qp
      e_d_coul=e_d_coul*e22/h*qp**2/(ma-2*jz)    
      
      
ccc   block energies doba
      e_rho0 =e_rho0*h*qp 
      e_rho1 =e_rho1*h*qp
      e_rho  =e_rho0+e_rho1

      e_deltarho0=e_deltarho0*h*qp
      e_deltarho1=e_deltarho1*h*qp
cc      print *, "E(Delta 1) = ", e_deltarho1
      e_deltarho =e_deltarho0+ e_deltarho1


      e_tau0=e_tau0*h*qp
      e_tau1=e_tau1*h*qp
      e_tau = e_tau0 + e_tau1

      e_nablaj0 = e_nablaj0*h*qp
      e_nablaj1 = e_nablaj1*h*qp
      e_nablaj  = e_nablaj0+e_nablaj1

      e_t0 = e_t0*h*qp
      e_t1 = e_t1*h*qp
      e_t  = e_t0+e_t1
        
      e_ext = e_ext *h*qp

      if (FRANCESCO) then
        e_bulk = e_bulk  * h * qp
!        e_surf = e_surf  * h * qp
      else
        e_bulk0=e_rho0+e_tau0
        e_bulk1=e_rho1+e_tau1
        e_bulk =e_rho +e_tau
      endif

      e_surf0=e_nablaj0*float(iso)+e_deltarho0+float(isj)*e_t0
      e_surf1=e_nablaj1*float(iso)+e_deltarho1+float(isj)*e_t1
!      e_surf =e_surf0+e_surf1
      e_surf = e_nablaj*float(iso)  + e_deltarho + float(isj)* e_t
  



      if(icount_chf.eq.1)then
      write(2,*) 'E(t0) = ',e_di_t0
      write(2,*) 'E(t3) = ',e_di_t3
      write(2,*) 'E(t1,t2) = ',e_di_t12
      write(2,*)
      end if


cc    complete E(rea)
      erea=erea* h*qp                                                   
      erea=erea+ ecbe/3.d0
      erea=erea/ma 
 

      EPOTE=EPOTE/MA                                                   
      ectot=ectot/ma        

cc    H0 = 1/2 [ E(kin)+ sum(eneergy eignevalues) + ext. pot. ] / A     
      h0 = 0.5*(EPOTE + ectot + i_ext*e_ext/ma )


ccc   total energy/A
ccc   H0:  sum of eigenvalues and kinetic energy; EREA: rearrangement energy      
      etot = h0 + erea                                                  
      write(2,*) 'h0,erea,etot=',h0,erea,etot


cc    pairing
      epair=0.d0
      sumuv=0.d0
      epair2=0.d0
c     write(2,*) 'NPO,nos=',NPO,nos
      do i=1,nos
       ctr=0.d0
       ctr2=0.d0
       if(dabs(fo(i)).le.1e-10)fo(i)=0.d0
       if(fo(i).le.-1e-5)stop '>>> ERROR IN THE OCCUPATION FACTORS'
       if(bcsp.eq.1.and.i.le.NPO)then 
        ctr=-deltap(i)*dsqrt(fo(i))*dsqrt(1.d0-fo(i))/2.d0
     1  *float(lj(i)+1)
       end if
       if(bcsn.eq.1.and.i.gt.NPO)then  
        ctr=-deltan(i)*dsqrt(fo(i))*dsqrt(1.d0-fo(i))/2.d0
     1  *float(lj(i)+1)
        ctr2=0.d0
        do j=NPO+1,nos
        ctr2=ctr2-vneu(i,j)*dsqrt(fo(i))*dsqrt(fo(j))*
     1  dsqrt(1.d0-fo(i))*dsqrt(1.d0-fo(j))/2.d0
        end do
       end if
       epair=epair+ctr
       epair2=epair2+ctr2
       if(bcsp.eq.1.and.i.le.NPO)sumuv=sumuv+dsqrt(fo(i))*
     1 dsqrt(1.d0-fo(i))*float(lj(i)+1)/2.d0
       if(bcsn.eq.1.and.i.gt.NPO)sumuv=sumuv+dsqrt(fo(i))*
     1 dsqrt(1.d0-fo(i))*float(lj(i)+1)/2.d0
      end do

      write(2,*) 'sumuv=',sumuv
      do i=NPO+1,nos
       delta_try=0.d0
       do j=NPO+1,nos
        delta_try = delta_try + dsqrt(fo(j))*dsqrt(1.d0-fo(j))
     1  /float(lj(i)+1)*vneu(i,j)
       end do
       write(2,*) i,delta_try
      end do


      WRITE(2,217) ectot,H0, erea, etot 
      etot = etot - e_di_constr/ma      
      write(2,*) 'epair=',epair,epair2
      write(2,*) 'EPAIR/A=',epair/ma,'ETOT/A=',etot+epair/ma
      write(2,*)
      write(2,*) 'ETOT=', etot*ma+epair
      write(2,*)



ccc   Total energy computed as an integral   (work HERE)
      e_di_tot = ectot*ma +e_di_t0+e_di_t3+e_di_t12 +  ! t0,t12,t3
     1 espo*float(iso) +      ! spin-orbit
     1 ecbd*float(icoul_d)+ ecbe*float(icoul_e) +      ! coulomb
     1 i_ext*e_ext                 ! external

cc    e_bulk: contribution of (custom) nuclear potential; ectot*ma: kin. en;
cc    esurf: gradient+spin-orbit;  ecbd+ecbe: Coulomb
      

      if (FRANCESCO) then
          e_di_tot = e_bulk + e_surf + ectot*ma +   
     1    ecbd*float(icoul_d) + ecbe*float(icoul_e) + i_ext*e_ext

          write(2,*) 'ETOT(francesco) = ',e_di_tot
          write(2,*) 'ETOT/A(francesco) = ', e_di_tot/ma

          print '(2A5)', "A","Z"
          print '(2I5)', int(aa), int(zz)
          print '(A10,F10.2)', "E(tot)", e_di_tot
          print '(A10,F10.2)', "E(tot)/A", e_di_tot/ma
          print '(A10,F10.2)', "E(kin)", ectot*ma
          print '(A10,F10.2)', "E(bulk)", e_bulk
          print '(A10,F10.2)', "E(surf)", e_surf
          print '(A10,F10.2)', "E(Coul)",
     2      ecbd*float(icoul_d) +ecbe*float(icoul_e)
          print '(A10,F10.2)', "E(so)", e_nablaj*float(iso)      
c          print '(A10,F10.2)', "Diff", espo - e_nablaj  
          if (i_ext.ne.0) then
             print *, "E(ext) = ", i_ext*e_ext
          end if

      end if


      if (debug) then
        if (FRANCESCO.eqv..false.) then
          print '(2A5)', "A","Z"
          print '(2I5)', int(aa), int(zz)
          if (iso.eq.0) then
             print *, "No spin orbit included"
          end if
          print '(A10,F10.2)', "E(tot)", e_di_tot 
c          print *, "ETOT(alternativo) = ", ectot*ma + ecbd
c     1     +ecbe + e_bulk + e_surf
          print '(A10,F10.2)', "E(tot)/A", e_di_tot/ma
          print '(A10,F10.2)', "E(kin)", ectot*ma
          print '(A10,F10.2)', "E(bulk)", e_bulk
          print '(A10,F10.2)', "E(surf)", e_surf
          print '(A10,F10.2)', "E(coul)",
     2      ecbd*float(icoul_d)+ecbe*float(icoul_e)
          print '(A10,F10.2)', "E(so)", espo*float(iso)  
cc          print '(A10,F10.2)',"E(so)-E(nablaJ)",espo*float(iso)-e_nablaj
          if (i_ext.ne.0) then      
           print *, "E(ext) = ", i_ext*e_ext
          end if
          print *, "E(rho) =", e_rho, "    E(tau) =", e_tau
          print *, "E(delta_rho) =", e_deltarho, "   E(nablaJ) =",
     1    e_nablaj,"      E(J2) =", e_t


          write (1897,'(2A5)') "A", "Z" 
          write (1897,'(2I5)') int(aa), int(zz)
          write (1897,'(4A10)'), "#", "is", "iv", "tot"
          write (1897,'(A10,2F10.2)') "Rho", e_rho0, e_rho1
          write (1897,'(A10,2F10.2)') "Tau", e_tau0, e_tau1
          write (1897,'(A10,2F10.2)') "Delta", e_deltarho0,e_deltarho1
          write (1897,'(A10,2F10.2)') "So",
     2      e_nablaj0*float(iso),e_nablaj1*float(iso)
          write (1897,'(A10,2F10.2)') "t",
     3       e_t0*float(isj),e_t1*float(isj)

          write (1897,'(A)')
          write (1897,'(A10,3F10.2)') "Bulk",e_bulk0,e_bulk1,e_bulk
          write (1897,'(A10,3F10.2)') "Surf",e_surf0,e_surf1,e_surf
          write (1897,'(A10,F10.2)') "Kin", ma*ectot
          write (1897,'(A10,F10.2)') "E(Coul)",
     3       ecbd*float(icoul_d) +ecbe*float(icoul_e)         
          write (1897,'(A10,F10.2)') "Tot",  
     2       e_bulk+e_surf+ectot*ma +
     2       ecbd*float(icoul_d) +ecbe*float(icoul_e)

        end if
      end if


c         Writing contributions to file
         open(unit=3020,file='temp2.out')
         write(3020,'(13A10)') "tot","kin","rho0","rho1",
     4       "tau0","tau1","delta0","delta1","so0","so1","coul",
     4       "bulk", "surf"
         write(3020,'(13F10.2)') 
     5    e_di_tot, ma*ectot,  e_rho0, e_rho1,
     5    e_tau0,e_tau1,e_deltarho0,e_deltarho1,
     5    e_nablaj0*float(iso),e_nablaj1*float(iso),
     5    ecbd*float(icoul_d) +ecbe*float(icoul_e),
     5    e_bulk,e_surf
         close(3020)


            
c     Comparing direct integration and sum of eigenvalues +rearrangement
      write(2,*) 'E (direct integration) = ',e_di_tot 
      diffperc=dabs(e_di_tot-etot*ma)/dabs(etot*ma)
      write(2,*) 'Delta (%) = ',diffperc*100.d0



      virial_th=2.d0*ectot*ma+3.d0*e_di_t0+
     15.d0*(e_di_t12+espo)+(3.d0*alfe+3.d0)*e_di_t3+
     1ecbd+ecbe

      write(2,*) 'virial theorem = ',virial_th
      am3=0.5d0*6874.4372d0*(2.d0*ectot*ma+
     16.d0*e_di_t0+20.d0*(e_di_t12+espo)+
     1(3.d0*alfe+3.d0)*(3.d0*alfe+2.d0)*e_di_t3)
c    1/16.d0/datan(1.d0)
      write(2,*) 'm(3) for the monopole = ',am3
      write(2,*)
      write(2,*) 'Direct Coulomb displacement energy',e_d_coul
c     end if

      etot = etot + epair/ma 
      e_di_tot = e_di_tot + epair 
c     - lamn*(ma-zz)/2.d0/ma - lamp*zz/2.d0/ma





CCC                                                                             
CCC   Calculating charge density                                        
CCC
      DO 1503 J=1,NNN                                                  
      X=RM(J)                                                          
      DC(J)=DSC(X)
 1503 continue                                                         
CCC
CCC   Calculating RADII                                        
CCC
      RN=0.                                                            
      RT=0.                                                            
      RP=0.                                                            
      RC=0.                                                            
      DO J=1,NNN                                                   
      X4=RM2(J)*RM2(J)                                                 
      RN=RN+DN(J)*X4                                                   
      RT=RT+DT(J)*X4                                                   
      RP=RP+DP(J)*X4                                                   
      RC=RC+DC(J)*X4                                                   
      end do

      RN=SQRT(QP*(H*RN/(MA-JZ)))                                       
      RT=SQRT(QP*(H*RT/MA))                                            
      if(jz.ne.0)RP=SQRT(QP*(H*RP/JZ))                                 
      if(jz.ne.0)RC=SQRT(QP*(H*RC/JZ))                                 
      r2lambda=rt**2

      WRITE(2,221) RN,RP,RC          
      if (FRANCESCO.or.debug) then
        print *, "Rp=",rp,"     Rn=",rn,"     Rc=", rc
      end if

      write(2,*) 'radii squared',rn**2,rp**2,rc**2
      write(2,*) 'SD sum rule',9.d0/QP*((MA-JZ)*rn**2-jz*rp**2)
      WRITE(2,213)                                                   
      WRITE(2,212) (DN(J),J=1,NNN)                                   
      WRITE(2,214)                                                   
      WRITE(2,212) (DP(J),J=1,NNN)                                   
      WRITE(2,216)                                                   
      WRITE(2,212)(DC(J),J=1,NNN)                                    
      WRITE(2,1200)                                                  
      WRITE(2,212) (DT(J),J=1,NNN)                


      open(unit=2020,file='temp.out')
      write(2020,'(2I6,2E15.5)') int(aa),int(zz),e_di_tot,rc
      close(2020)

      open(unit=13,position='append',file='be.dat')
c     open(unit=13,file='be.dat')
      write(13,1210) ma,ma*etot+epair,rc
      close(13)
      IF(IDCC.EQ.1) write(12,1209)(DC(J),J=1,NNN)                      
 1209 FORMAT(5E16.9)                                                   
 1210 format(1x,f5.1,2(3x,e12.5))
      LJ1=1                                                            
 7001 RP=0.                                                            
      RPN=0.                                                           
      DO7000 J=1,NNN                                                   
      X4=RM2(J)**LJ1  
      x4pn=rm(j)**(2+lj1)                                              
      RP=RP+DT(J)*X4                                                   
      RPN=RPN+DN(J)*DP(J)*x4pn                                          
 7000 CONTINUE                                                         
      RP=QP*H*RP/MA                                                    
      WRITE(2,7010) LJ1,RP                                           
      r2l(lj1)=rp
      RPN=QP*H*RPN                                                     
      WRITE(2,7011) LJ1,RPN                                            
 7011 FORMAT(2X,'INTEGRALE DE RO(N)*RO(P) AVEC R**',I3,'=',E13.5)
      DipoKap=0.d0
      skyterm=(t1*(1.d0+(x1/2.d0)))+(t2*(1.d0+(x2/2.d0)))
      denomin=hb0*ma
      DO 149 J=1,NNN
      cuca=QP*RM2(J)*DP(J)*DN(J)
      DipoKap=DipoKap+skyterm/denomin*cuca*H
  149 CONTINUE
c     write(*,*)'Kdip'
c     write(*,*)DipoKap
 7010 FORMAT(/2X,'R**(2*',I3,'-2)=',E13.5)                             
      LJ1=LJ1+1                                                        
      IF(LJ1.LE.12) GO TO 7001       

CCC                                                                             
CCC   RECOUVREMENTS                                                             
CCC

c modif ek le 22/09/97

c     go to 1000

      SUM = 0.                                                         
      DO 150 J=1,NOS                                                  
      SUM=SUM+DEG(J)*(2.*NN(J)+LL(J)-.5)
  150 continue                                        
c     write(6,*) MA,SUM
      SUM = SQRT(SUM/MA)                                              
      B = RT/SUM                                                      
      ECM = .75*41.47/(B*B*MA)                                        
      WRITE(2,231) ECM,B,RT                                         
  231 FORMAT(/,' ECM/A=', E12.5, '   B=',E12.5, '   RT=',E12.5/)      
      DO 153  J=1,NOS                                                 
      L = LL(J)                                                       
      DO 152 JJ=1,8                                                   
      REC(JJ) = 0.d0
  152 continue                                                       
      E = 0.d0                                                      
      DO 155  K = 1,8                                                  
      N = K - 1                                                        
      DO 154 NREC = 1,NNN                                              
      X = H * NREC                                                     
      UNLOSC(NREC,J) = HWF(N,L,B,X)                                    
      REC(K) = REC(K) + H*UNL(NREC,J)*UNLOSC(NREC,J)
  154 continue                            
      E = E + REC(K)*REC(K)
  155 continue                                                     
      WRITE(2,229) ST(J) , LT(J)                                     
      WRITE(2,211) (REC(M),M=1,8),E                                  
  153 CONTINUE                                                         
  229 FORMAT(5X,A6,I2,'  REC=',E12.5)                                  
      DO 156 I = 1,NOS                                                 
      L = LL(I)                                                        
      N = NN(I) - 1                                                   
      DO 1560  J=1,NNN                                                 
      X = H*J                                                          
      UNLOSC(J,I) = HWF(N,L,B,X)
 1560 continue                                                
  156 CONTINUE                                                         
      DO4500 J=1,NNN                                                   
      XR=H*FLOAT(J)                                                   
      DF(J)=(3.*DT(J)+XR*DT1(J))*XR*XR
 4500 continue                                          
      write(2,*)
      WRITE(2,4501)                                                  
 4501 FORMAT(2X,' TASSIE MONOPOLE TRANSITION DENSITY TIMES R**2'/)    
      WRITE(2,211) (DF(J),J=1,NNN)                                 
CCC
CCC   CALCUL DE (D**2)                                                          
CCC
      DSQ=RT*RT*MA/12.                                              
      DO1003 IA=1,NOS                                               
      LTA=LT(IA)                                                    
      LA=LL(IA)                                                     
      JA=LJ(IA)                                                     
      DO 1004 IB=1,NOS                                               
      IF(LT(IB)-LTA) 1004,1008,1004                                  
1008  LAB=LL(IB)+LA+1                                               
      IF(MOD(LAB,2)) 1004,1007,1004                                   
1007  JB=LJ(IB)                                                       
      RAB=0.                                                          
      DO1005 IX=1,NNN                                                 
      RAB=RAB+UNL(IX,IA)*UNL(IX,IB)*IX*H*H                            
 1005 continue
c     ANS=F3J(JA,JB,2,1,-1,0)
      ans=0.                                                   
      DSQ=DSQ-(JA+1.)*(JB+1.)*RAB*RAB*ANS*ANS/12.                    
 1004 CONTINUE                                                      
 1003 CONTINUE                                                     
      WRITE(2,1006) DSQ                                         
 1006 FORMAT(75('*')/2X,'DSQUARE=',E12.5/75('*'))              
  321 format(2X,F5.2,2X,F6.3,2X,E11.4,2X,E11.4)
C_GL_17_6
      write(2,*)
      write(2,*)'delta medio finale protonico',deltamediop
      write(2,*)'delta medio finale neutronico',deltamedio
      write(2,*)
c     write(21,*)'en.fermi:prot  neutr'
c     write(21,*)'       ',lamp,lamn
c     write(21,*)'posiz','   ','densita','   '  
c    1           ,'unl21','   ','unl20'

      open(unit=22,status='unknown',file='vpair.dat')
      do ir=1,nnn
cc SCRIVO VPAIR(R)=(VZERO/2)*(1-RHO(R)/DENSC)**GAMMA
cc (a meno del segno !!)
      write(22,*)ir,-vpair(ir)
      enddo
      close(22)

      close(29) 
      close(51)
      close(50)
      close(28)
      close(27)
      close(17)
      close(56)
      close(16)
      close(12)
      close(10)
      return



cc    FRANCESCO, CUSTOM interaction
      contains

       

       function mean_field(rho,beta,q)
       real(8) rho,beta, q, mean_field
       real(8) c_SM, c_sym, ee
       mean_field=0.0d0

       if (rho.ge.1.0E-12) then 
         do k=1,n_terms
           c_SM = coeff_SM(k)
           c_sym = coeff_sym(k)
           eee = exponents(k)

           coeff = (eee+1.)*c_SM
     &     + 2.*beta*(q-beta)*c_sym 
     &     + (eee+1.)*beta**2 *c_sym

           mean_field = mean_field + coeff * rho**eee
         end do
       end if
c       print *, "Rho    ", rho, "           Mean field    ",  mean_field
       end function mean_field


       function pot_density(rho,beta)
       real(8) rho,beta,pot_density
       pot_density=0.
       if (rho.ge.1.0E-12) then 
c       integer g
         do k=1,n_terms
           ccc = coeff_SM(k) + beta**2 *coeff_sym(k)
           eee = exponents(k)
           pot_density=pot_density + ccc* rho**(eee+1)
         end do
       end if
       end function pot_density


ccc   Contribution   c_k(beta) rho^{gamma_k +1) 
       function one_term(rho,beta,k)
         real(8) rho,beta,one_term
         integer k
         if (rho.ge.1.0E-12) then 
           ccc = coeff_SM(k) + beta**2 *coeff_sym(k)
           eee = exponents(k)
           one_term=  ccc* rho**(eee+1)
         end if
       end function one_term
      


       function rearrangement_density(rho,beta)
       real(8) rho,beta, rearrangement_density
       rearrangement_density = 0.0d0
       if (rho.ge.1.0E-12) then
         do k=1, n_terms
           eee = exponents(k)
           coeff = coeff_sm(k) + (beta**2) * coeff_sym(k)
           rearrangement_density=rearrangement_density + 
     &     (1.-eee)/2. * coeff * rho**(eee+1)
         end do
       end if

       end function rearrangement_density


       function rearrangement_term(rho,beta,q)
       real(8) rearrangement_term,rho,beta, q
       rearrangement_term = 0.
       if (rho.ge.1.0E-12) then
         do k=1, n_terms
           c_SM = coeff_SM(k)
           c_sym = coeff_sym(k)
           eee = exponents(k)

           coeff = (eee-1.)*c_SM
     &     + 2.*beta*(q-beta)*c_sym
     &     + (eee-1.)*beta**2 *c_sym

           rearrangement_term = rearrangement_term + coeff * rho**eee
         end do
       end if
       end  function rearrangement_term


cc     V(r) = 1/2 m omega^2 r^2 
       function harmonic_potential(r,omega)
         real(8) harmonic_potential, r,omega    ! omega = [MeV]
         harmonic_potential = amc2/2. * (omega/hbc*r)**2 
       end function harmonic_potential


      END subroutine bcsgen  
     
      


      FUNCTION DSC(Z1)                                                 
CCC                                                                             
CCC   DENSITE DE CHARGE                                                         
CCC   FACTEUR DE FORME GAUSSIEN                                                 
CCC   MU=0.65 FM                                                                
CCC                                                                             
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      DIMENSION YK(10),PI(10),SIG(2)                                            
      DATA (YK(I),I=1,10)/0.24534071d0,0.73747373d0,1.2340762d0,
     &   1.7385377d0,2.2549740d0,2.7888061d0,3.3478546d0,
     &  3.9447640d0,4.6036824d0,5.3874809d0/                  
      DATA(PI(I),I=1,10)/0.49092150d0,0.49384339d0,0.49992087d0,
     &  0.50967903d0,0.52408035d0,0.54485174d0,0.57526244d0,
     &  0.62227870d0,0.70433296d0,0.89859196d0/        
      DATA SIG/1.,-1./                                                          
      X=Z1                                                                      
      XMU=0.65                                                                  
      ZW=1./(0.4431125*XMU**3)                                                  
      SG=0.                                                                     
      XSM=X/XMU                                                                 
      DO 100 I=1,10                                                             
      DO 1000 K=1,2           
      Y=X+XMU*YK(I)*SIG(K)                                                      
      IF(Y)100,100,101                                           
  101 F=ZW*GUGUS(Y)*Y*Y                                                         
      YSM=Y/XMU                                                                 
      SG=SG+XMU*PI(I)*F*VKG(0,XSM,YSM)    
 1000 continue                                      
  100 CONTINUE                                                                  
      DSC=SG                                                                    
      RETURN                                                                    
      END  FUNCTION DSC                                                                     


      subroutine DERIV                                                          
C     FIVE POINT DIFFERENTATION FORMULA (ABRAMOWITZ, PAG 914)                   
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      PARAMETER (nnp=10000,NEOC=100)
      COMMON/DER/FCT(NNP),DF(NNP),H                                             
      COMMON/BMAX/NMAXT
      DIMENSION A(5,5)                                                          
      DATA A(1,1),A(1,2),A(1,3),A(1,4),A(1,5)                                   
     +                  /-50.,96.,-72.,32.,-6./                                 
      DATA A(2,1),A(2,2),A(2,3),A(2,4),A(2,5)                                   
     +                  /-6.,-20.,36.,-12.,2./                                  
      DATA A(3,1),A(3,2),A(3,3),A(3,4),A(3,5)                                   
     +                  /2.,-16.,0.,16.,-2./                                    
      DATA A(4,1),A(4,2),A(4,3),A(4,4),A(4,5)                                   
     +                     /-2.,12.,-36.,20.,6./                                
      DATA A(5,1),A(5,2),A(5,3),A(5,4),A(5,5)                                   
     +                     /6.,-32.,72.,-96.,50./                               
      DATA EMFACT/24./                                                          
      
      NNN=NMAXT                                                             
      NMX=NNN-2                                                                 
      DO J=1,NNN                                                              
       K=3                                                                       
       IF(J.LT.3)K=J                                                             
       IF(J.GT.NMX)K=J-NNN+5                                                     
       SUM=0.                                                                    
       do I=1,5                                                                
        JJ=J+I-K                                                                  
        SUM=SUM+A(K,I)*FCT(JJ)
       end do                                                    
       DF(J)=SUM/(H*EMFACT)
      end do                                                    
      return                                                                    
      end subroutine DERIV


      FUNCTION GUGUS(X)                                                         
C---- INTERPOLATION QUADRATIQUE DE FONCTIONS PAIRES SUR 0-RMAX                  
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      PARAMETER (nnp=10000,NEOC=100)
      COMMON/DER/FCT(NNP),DF(NNP),H                                             
      COMMON RM(NNP),DU(NNP),VC(NNP),U(NNP),VS(NNP),TAUN(NNP),                  
     +TAUP(NNP),TAUT(NNP),DJN1(NNP),DJP1(NNP),DJT1(NNP),VN(NNP),VP(NNP),        
     +VSON(NNP),VSOP(NNP),ECOULB(NNP),VQRE(NNP),RM2(NNP),RM3(NNP),              
     +HMEN(NNP),HMEP(NNP),DHMEN(NNP),DHMEP(NNP),D2HMEN(NNP),D2HMEP(NNP),        
     +DN(NNP),DP(NNP),DT(NNP),DN1(NNP),DP1(NNP),DT1(NNP),DN2(NNP),              
     +DP2(NNP),DT2(NNP),DJN(NNP),DJP(NNP),DJT(NNP),DC(NNP),AVP(NNP,             
     +NEOC),DUNL(NNP,NEOC),UNL(NNP,NEOC),DVC(NNP)   
      COMMON/BMAX/NMAXT
      NNN=NMAXT                                                                 
      X2=2.*H                                                                   
      IF(X.GT.X2) GO TO 2                                                       
      GUGUS=((4.*DP(1)-DP(2))+(DP(2)-DP(1))*X*X)/3.                             
      RETURN                                                                    
    2 I=X/H                                                                     
      IF(I.GE.NNN) GO TO 4                                                      
      Y=(X-I*H)/H                                                               
      S=DP(I)+.5*Y*(DP(I+1)-DP(I-1))                                            
      GUGUS=S+.5*Y*Y*(DP(I+1)+DP(I-1)-2.*DP(I))                                 
      RETURN                                                                    
    4 GUGUS=.0                                                                  
      RETURN                                                                    
      END FUNCTION GUGUS


      SUBROUTINE NUM1L(N,H,E,S2,U,S,NO,EPS)                                     
CCC   VERSION CORRIGEE LE 21 NOV 72                                             
C*****INTEGRATION DE L'EQUATION DE SCHROEDINGER PAR LA METHODE DE NUMERO        
C*****POUR E NEGATIF                                                            
C*****RECHERCHE DE L'ENERGIE PROPRE PAR LA METHODE DE RAPHSON-NEWTON            
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      DIMENSION U(1),S(1)                                                       
      DATA RAP1,RAP2/0.,0./                                                     
      H12=H*H/12.                                                               
C*****CONTROLE DES CONDITIONS ASYMPTOTIQUES                                     
      IF(E.GT..0) E=.0                                                          
      DEI=.0                                                                    
      EPSS=.1E-10                                                               
      IF(U(N-1).GT.EPSS) GO TO 10                                               
      DEI=U(N-1)-EPSS                                                           
      DO 8 K=1,N                                                                
      U(K)=U(K)-DEI
    8 continue     
   10 U(N)=U(N-1)                                                               
   

cc    Calcul du nombre d etats lies par integration a energie nulle

      S(1)=1.E-10                                                               
      B0=0.                                                                     
      AA=H12*U(1)                                                               
      IF (S2) 16,18,16                                                     
   16 B0=-S(1)*AA                                                               
   18 B1=S(1)*(1.-AA)                                                           
      DO 38 K=2,N                                                               
      B2=12.*S(K-1)-10.*B1-B0                                                   
      IF (ABS(B2).LT.1.E+10) GO TO 22                                           
      B2=B2*1.E-20                                                              
      B1=B1*1.E-20                                                              
   22 AA=H12*U(K)                                                               
      S(K)=B2/(1.-AA)                                                           
      B0=B1                                                                     
      B1=B2
   38 continue  
      DO 42 K=5,N                                                               
      N0=K                                                                      
      IF(U(K).LT.0.) GO TO 44                                                   
   42 CONTINUE                                                                  
   44 NEL=0                                                                     
      DO 52 K=N0,N                                                              
      IF (S(K-1)*S(K)) 46,50,52                                                 
   46 NEL=NEL+2                                                                 
      GO TO 52                                                                  
   50 NEL=NEL+1                                                                 
   52 CONTINUE                                                                  
      NEL=NEL/2                                                                 
      IF(NEL.GT.NO) GO TO 64                                                    
      IF(NEL.EQ.NO) GO TO 60                                                    
   62 NO=-1                                                                     
      RETURN                                                                    
   60 RAP1=S(N-1)/S(N)                                                          
      RAP2=EXP(H*SQRT(U(N-1)-E))                                                
      IF(RAP1.LT.RAP2) GO TO 62                                                 

C*****CALCUL DE EMIN ET EMAX ENTRE LESQUELLES SE TROUVE L ENERGIE PROPRE        
64    UMIN=U(1)                                                                 
      DO 70 K=2,N                                                               
      IF(U(K).LT.UMIN) UMIN=U(K)                                                
   70 CONTINUE                                                                  
      EMIN=UMIN                                                                 
      EMAX=0.                                                                   
C*****DEBUT DE LA RECHERCHE DE L'ENERGIE PROPRE DANS L'INTERVALLE MAXIMU        
      TE=EMAX-EMIN                                                              
C*****REJET DE L'ENERGIE D'ESSAI E PROPOSEE SI ELLE EST A L EXTERIEUR DE        
C*****BORNES (EMIN,EMAX)                                                        
      IF((E.LT.EMIN).OR.(E.GT.EMAX)) E=EMIN+.5*TE                               
      E1=EMIN                                                                   
      E2=EMAX                                                                   
      J=2                                                                       
      I=1                                                                       
      GO TO 102                                                                 
C*****REDUCTION DES BORNES EMIN ET EMAX                                         
   90 EMIN=E1                                                                   
      EMAX=E2                                                                   
      TE=EMAX-EMIN                                                              
      J=2                                                                       
   98 I=1                                                                       
  100 E=EMIN+TE*FLOAT(I)/FLOAT(J)                                               
  102 DE=0.                                                                     
  104 E=E+DE                                                                    
      IF(E.GT.0.) GO TO 204                                                     
      S(N)=1.E-10                                                               
      N1=N-1                                                                    
      RAP2=EXP(H*SQRT((U(N-1)+U(N))/2.-E))                                      
      S(N1)=S(N)*RAP2                                                           
      AA=H12*(U(N1)-E)                                                          
      B0=S(N)*(1.-AA)                                                           
      B1=S(N1)*(1.-AA)                                                          
      N1=N-2                                                                    
      DO 138 KAUX=1,N1                                                          
      K=N1-KAUX+1                                                               
      B2=12.*S(K+1)-10.*B1-B0                                                   
      AA=H12*(U(K)-E)                                                           
      S(K)=B2/(1.-AA)                                                           
      B0=B1                                                                     
      B1=B2                                                                     
      IF(U(K).LT.E) GO TO 140                                                   
  138 CONTINUE                                                                  
  140 N1=K                                                                      
C*****NORMALISATION DE LA FONCTION D ONDE A S(N1)                               
      DO 146 KAUX=N1,N                                                          
      K=N-KAUX+N1                                                               
      S(K)=S(K)/S(N1)
  146 continue                                                          
C*****DEBUT DE L INTEGRATION VERS L'EXTERIEUR JUSQU'A N1                        
      S(1)=1.E-10                                                               
      B0=0.                                                                     
      AA=H12*(U(1)-E)                                                           
      IF(S2) 156,158,156                                                        
  156 B0=-S(1)*AA                                                               
  158 B1=S(1)*(1.-AA)                                                           
      DO 170 K=2,N1                                                             
      B2=12.*S(K-1)-10.*B1-B0                                                   
      AA=H12*(U(K)-E)                                                           
      S(K)=B2/(1.-AA)                                                           
      B0=B1                                                                     
      B1=B2
  170 continue          
C*****NORMALISATION DE LA FONCTION A S(N1)                                      
      DO 174 K=1,N1                                                             
      S(K)=S(K)/S(N1)
  174 continue                                                           
C*****CALCUL DE LA CORECTION D ENERGIE                                          
      SOM=0.                                                                    
      DO 180 K=1,N                                                              
      SOM=SOM+S(K)*S(K)
  180 continue                                                         
      DE=((-S(N1-1)+2.-S(N1+1))/(H*H)+U(N1)-E)/SOM                              
      IF(ABS(DE).GT.EPS) GO TO 104                                              
C*****CALCUL DU NOMBRE DE NOEUDS DE L ETAT PROPRE TROUVE                        
      DO 182 K=5,N                                                              
      IF(U(K).LT.E) GO TO 184                                                   
  182 CONTINUE                                                                  
  184 N0=K                                                                      
      NEL=0                                                                     
      DO 192 K=N0,N1                                                            
      IF(S(K-1)*S(K)) 186,190,192                                               
  186 NEL=NEL+2                                                                 
      GO TO 192                                                                 
  190 NEL=NEL+1                                                                 
  192 CONTINUE                                                                  
      NEL=NEL/2                                                                 
C*****L ETAT PROPRE TROUVE EST-IL LE BON                                        
      IF(NEL-NO) 198,214,202                                                    
  198 IF(E.GT.E1) E1=E                                                          
      GO TO 204                                                                 
  202 IF(E.LT.E2) E2=E                                                          
  204 I=I+2                                                                     
      IF (I.LE.J)  GO TO 100                                                    
      J=2*J                                                                     
      IF(ABS(E1-EMIN).GT.EPS.OR.ABS(EMAX-E2).GT.EPS) GO TO 90                   
      GO TO 98                                                                  
C*****NORMALISATION DE LA FONCTION PROPRE                                       
  214 SOM=1./SQRT(SOM*H)                                                        
      DO 218 K=1,N                                                              
      S(K)=S(K)*SOM
  218 continue                                                             
      E=E+DEI                                                                   
      RETURN                                                                    
C*****DEBUT FORMATS                                                             
 2000 FORMAT(/,35X,56HL'ETAT DEMANDE N'EST PAS LIE. RETOUR DE NUM1L AVEC        
     1 NO=-1,/)                                                                 
C*****FIN FORMATS                                                               
      END SUBROUTINE NUM1L



      SUBROUTINE NUM2B(N,H,E,S2,U,RMS,S,NO,HB)                           
      parameter (nnp=10000)
CCC   VERSION DU 13 FEVRIER 1979                                                
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      DIMENSION U(nnp),RMS(nnp),S(nnp),P(nnp),T(nnp)                            
      common/bdeb/ideb0,ideb

      N1=N+1                                                                    
      H12=H*H/12.d0                              
      GO TO 1                                                                   
  100 NO=-1                                                                     
      RETURN                                                                    
c   1 E0=-2.8938d0                                                              
    1 E0=-4.3407d0
      DE=0.24115d0                                                             
    2 E0=E0+DE                                                                  
      DO 3 K=1,N1     
      P(K)=U(K)+RMS(K)*E0
    3 continue                                                       
      UMIN=P(1)                                                                 
      DO4 K=2,N1                                                                
      IF(P(K).LT.UMIN) UMIN=P(K)                                                
    4 CONTINUE                                                                  
      DO63 IE=1,2                                                               
      KK=0                                                                      
      GO TO(64,65),IE                                                           
   64 E0=UMIN-DE                                                                
      GO TO 5                                                                   
   65 DE=0.24115                                                                
      E0=E1-DE                                                                  
    5 E0=E0+DE                                                                  
      KK=KK+1                                                                   
      IF(KK.GT.50) GO TO 100                                                    
      DO 6 K=1,N1    
      P(K)=U(K)+RMS(K)*E0
    6 continue                                                       
      S(1)=1.E-10                                                               
      B0=0.d0                                                                   
      AA=H12*(P(1)-E0)                                                          
      IF(S2) 7,8,7                                                              
    7 B0=-S(1)*AA                                                               
    8 B1=S(1)*(1.d0-AA)                                                           
      DO 9 K=2,N1 

      B2=12.d0*S(K-1)-10.d0*B1-B0                                                   
      AA=H12*(P(K)-E0)                                                          
      S(K)=B2/(1.d0-AA)                                                       
      B0=B1                                                                     
      B1=B2
    9 continue                                                                     
      NEL=0                                                                     

c      write(6,*) S(N),S(N-1)
      DO10 K=8,N                                                                


c      IF((S(K-1) .gt. 1e18) .and. (S(K) .gt. 1e18 )) goto 10
      if (S(K-1)*S(K)) 11,12,10



c      else if ((S(K-1)*S(K)) .eq. 0) then
c        go to 12
c        endif 

   11 NEL=NEL+2                                                                 
      GO TO 10                                                                  
   12 NEL=NEL+1                                                                 
   10 CONTINUE                                                                  
      NEL=NEL/2                                                                 
      GO TO(66,67),IE                                                           
   66 IF(NEL-NO) 5,13,14                                                        
   13 E1=E0                                                                     
      SE1=S(N)                                                                  
      NEL1=NEL                                                                  
      GO TO 63                                                                  
   14 E0=E0-DE                                                                  
      DE=0.5*DE                                                                 
      GO TO 5                                                                   
   67 IF(NEL-NO-1) 5,15,14                                                      
   15 E2=E0                                                                     
      SE2=S(N)                                                                  
      NEL2=NEL                                                                  
   63 CONTINUE                                                                  
      KK=0                                                                      
      DO 25 K=1,N1    
      T(K)=S(K)
   25 continue                                                                 
   36 IF(NEL2-NEL1-1) 18,17,18                                                  
   18 ET1=E1*HB                                                                 
      ET2=E2*HB                                                                 
      WRITE(2,19) ET1,SE1,NEL1,ET2,SE2,NEL2                                         
   19 FORMAT(2X,2E12.5,I5,2E12.5,I5)                                            
      GO TO 100                                                                 
   17 IF(SE1*SE2) 21,20,18                                                      
   20 IF(SE1.EQ.0.d0) GO TO 18                                                    
   23 E=E2                                                                      
      NEL=NEL2                                                                  
      GO TO 101                                                                 
   21 DE=0.5d0*(E2-E1)                                                            
      DEB=DE*HB                                                                 
      IF(DEB.GT.1.E-12) GO TO 22                                                
c     IF(DEB.GT.1.E-6) GO TO 22                                                
      GO TO 23                                                                  
   22 E3=E1+DE                                                                  
      KK=KK+1                                                                   
      IF(KK.GT.50) GO TO 18                                                     
      DO 24 K=1,N1     
      P(K)=U(K)+RMS(K)*E3
   24 continue                                                       
      S(1)=1.E-10                                                               
      B0=0.d0                                                                   
      AA=H12*(P(1)-E3)                                                          
      IF(S2) 26,27,26                                                           
   26 B0=-S(1)*AA                                                               
   27 B1=S(1)*(1.d0-AA)                                                           
      DO 28 K=2,N1      
      B2=12.d0*S(K-1)-10.d0*B1-B0     
      AA=H12*(P(K)-E3)                                                          
      S(K)=B2/(1.d0-AA)                                                           
      B0=B1                                                                     
      B1=B2
   28 continue                                                       
      SE3=S(N)                                                                  
      NEL3=0                                                                    
      DO29 K=8,N                                                                
      IF(S(K-1)*S(K)) 31,30,29                                                  
   31 NEL3=NEL3+2                                                               
      GO TO 29                                                                  
   30 NEL3=NEL3+1                                                               
   29 CONTINUE                                                                  
      NEL3=NEL3/2                                                               

c      write(6,*) SE1,SE3
c      if ((SE1 .gt. 1e10) .and. (SE3 .gt. 1e10)) go to 34
      IF(SE1*SE3) 32,33,34                                                      
   32 NEL2=NEL3                                                                 
      SE2=SE3                                                                   
      E2=E3                                                                     
      DO 35 K=1,N1     
      T(K)=S(K)
   35 continue                                                                 
      GO TO 36                                                                  
   33 E=E3                                                                      
      NEL=NEL3                                                                  
      DO 37 K=1,N1        
      T(K)=S(K)
   37 continue                                                                 
      GO TO 101                                                                 
   34 NEL1=NEL3                                                                 
      SE1=SE3                                                                   
      E1=E3                                                                     
      GO TO 36                                                                  
  101 NO=NEL-1                                                                  
      DO43 K=N,4,-1                                                             
      IF(T(K-1)*T(K)) 44,45,43                                                  
   45 KK=K-2                                                                    
      GO TO 46                                                                  
   44 KK=K-1                                                                    
      GO TO 46                                                                  
   43 CONTINUE                                                                  
   46 IF(E) 47,48,48                                                            
   47 DO53 K=KK,4,-1                                                            
      IF(dABS(T(K-1))-dABS(T(K))) 54,54,53                                        
   53 CONTINUE                                                                  
   54 KL=(K+KK)/2                                                               
      K1=N-1                                                                    
      K2=N-2                                                                    
      S(N)=0.d0                                                         
      S(K1)=1.E-10                                                              
      AA=H12*(P(K1)-E)                                                          
      B0=0.d0                                                                  
      B1=S(K1)*(1.-AA)                                                          
      K4=K1-KL                                                                  
      DO55 K3=1,K4                                                              
      K=K2-K3+1                                                                 
      B2=12.*S(K+1)-10.*B1-B0                                                   
      AA=H12*(P(K)-E)                                                           
      S(K)=B2/(1.d0-AA)                                                           
      B0=B1                                                                     
      B1=B2                                                                     
   55 CONTINUE                                                                  
      FAC=T(KL)/S(KL)                                                           
      DO56 K=KL,N                                                               
      T(K)=S(K)*FAC
   56 continue                                                             
      GO TO 50                                                                  
   48 DO 51 K=KK,N         
      T(K)=T(K)-FLOAT(K-KK)*T(N)/FLOAT(N-KK)
   51 continue                                    
   50 CONTINUE      !modif ek le 25 09 97
      SOM=0.d0                                                              
      DO38 K=1,N                                                                
      Y=1.d0-RMS(K)                                                               
      T(K)=T(K)*dSQRT(Y)                                                         
      SOM=SOM+T(K)*T(K)
   38 continue                                                         
      SOM=SQRT(SOM*H)                                                           
      SIG=T(10)/dABS(T(10))                                                      
      DO39 K=1,N                                                                
      S(K)=SIG*T(K)/SOM
   39 continue                                                         
      RETURN                                                                    
      END SUBROUTINE NUM2B



      FUNCTION VKG(I1,Z1,Z2)                                                    
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      DIMENSION FJ(50)                                                          
      L=I1                                                                      
      X=Z1                                                                      
      Y=Z2                                                                      
      DELTA=Y-X                                                                 
      A=2.*X*Y                                                                  
      IF(A-15.)101,101,113                                                      
  101 IF(A-0.01)102,102,104                                                     
  102 ARG=EXP(-(X*X+Y*Y))                                                       
      FJ0=ARG*(1.+A*A/6.)                                                       
      IF(L-1)103,103,105                                                        
  103 FJ(1)=FJ0                                                                 
      FJ(2)=ARG*A*(10.+A*A)/30.                                                 
      GO TO 107                                                                 
  104 U=DELTA*DELTA                                                             
      V=(X+Y)*(X+Y)                                                             
      U=EXP(-U)                                                                 
      V=EXP(-V)                                                                 
      FJ0=(U-V)/(2.*A)                                                          
  105 L2=L+5                                                                    
      L2=L+10                                                                   
      FJ(L2)=1.E-10*FLOAT(2*L2+1)/A                                             
      FJ(L2+1)=1.E-10                                                           
      L3=L2-1                                                                   
      DO 106 LL=1,L3                                                            
      L1=L2-LL                                                                  
      FJ(L1)=FLOAT(2*L1+1)*FJ(L1+1)/A+FJ(L1+2)                                  
      IF(FJ(L1)-1.E+30)106,106,111                                              
  111 DO 112 L4=L1,L2                                                           
      FJ(L4)=1.E-10*FJ(L4)      
  112 continue                                                
  106 CONTINUE                                                                  
      ZZ=FJ0/FJ(1)                                                              
      L2=L2-9                                                                   
      DO 109 L1=1,L2                                                            
      FJ(L1)=ZZ*FJ(L1)     
  109 continue                                                     
      GO TO 107                                                                 
  113 U=DELTA*DELTA                                                             
      V=(X+Y)*(X+Y)                                                             
      U=EXP(-U)                                                                 
      V=EXP(-V)                                                                 
      FJ(1)=(U-V)/(2.*A)                                                        
      FJ(2)=((A-1.)*U+(A+1.)*V)/(2.*A*A)                                        
      IF(L-1)107,107,114                                                        
  114 DO 115 L1=2,L                                                             
      FJ(L1+1)=-FLOAT(2*L1-1)*FJ(L1)/A+FJ(L1-1)    
  115 continue                             
  107 VKG=FLOAT(L+L+1)*FJ(L+1)                                                  
      RETURN                                                                    
      END  FUNCTION VKG                                                                     

      FUNCTION HWF(I1,I2,Z1,Z2)                                                 
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      N=I1                                                                      
      L=I2                                                                      
      B=Z1                                                                      
      R=Z2                                                                      
      X=R/B                                                                     
      XSQ=-2.*X*X                                                               
      Y=-X*X/2.                                                                 
      B3=B*B*B                                                                  
      RHO=1.                                                                    
      VNL=1.                                                                    
      IF(N)101,101,102                                                          
  102 XA=1.                                                                     
      K=2*L+1                                                                   
      NS=N+1                                                                    
      NI=0                                                                      
      DO I=1,N                                                              
        K=K+2                                                                     
        NI=NI+1                                                                   
        NS=NS-1                                                                   
        XA=XA*XSQ*NS/(K*NI)                                                       
        VNL=VNL+XA 
        RHO=(K*RHO)/(2.*NI)                                                       
      end do                                                             
  101 DPL=1.                                                                    
      XSL=1.                                                                    
      IF(L)104,104,105                                                          
  105 K=1                                                                       
      DO I=1,L                                                              
       K=K+2                                                                     
       XK=K                                                                      
       DPL=2.*DPL/XK                                                             
       XSL=XSL*X
      end do                                                             
  104 ARG=RHO*DPL/(B3*1.772454)                                                 
      XA=2.*XSL*R*SQRT(ARG)                                                     
      HWF=XA*VNL*EXP(Y)                                                         
      RETURN                                                                    
      END                                                                        


      function yl(la,ja,lb,jb,l)
      implicit double precision (a-h,o-z)
ccc   ja and jb are twice the true values
ccc   other variables have their true values
      qp=4.d0*3.14159265d0
      yl=0.d0
C_GL_2_4_05
      lmina=iabs(ja-jb)/2
      lmaxa=(ja+jb)/2
      if(l.lt.lmina)go to 1
      if(l.gt.lmaxa)go to 1
      ipar=1-2*mod(la+lb+l,2)
      if(ipar.lt.0)go to 1
      l1p=la+lb
      l1m=iabs(la-lb)
      l2p=(ja+jb)/2
      l2m=iabs(ja-jb)/2
      lp=min0(l1p,l2p)
      lm=max0(l1m,l2m)
      l3=(l-lp)*(l-lm)
      l4=la+lb+l-2*((la+lb+l)/2)
      ld=2*l
      ig=l+(jb-1)/2
      sg=1.d0-2.d0*mod(ig,2)
      xja=ja/2.d0
      xjb=jb/2.d0
      xl=l
      jpha=(ja+jb)/2+1
      jpha=1-2*mod(jpha,2)
      c=(ja+1.d0)*(jb+1.d0)/qp
      sg=sg*jpha
      yl=sg*dsqrt(c)*cofcg(xja,xjb,xl,0.5d0,-0.5d0,0.d0)
    1 return
      end  function yl

      function cofcg(a,b,c,x,y,z)
      implicit double precision (a-h,o-z)
      xj1=a
      xj2=b
      xj3=c
      xm1=x
      xm2=y
      xm3=-z
c-rcnp
      xrcnp=dabs(xj1)+dabs(xj2)+dabs(xj3)
      if(xrcnp.lt.200.d0) go to 9999
      write (99,*) 'xj1, xj2, xj3; a, b, c  at point5 =',xj1,xj2,xj3,
     1a,b,c
      stop
 9999 continue
c-rcnp
      call troisj(xj1,xj2,xj3,xm1,xm2,xm3,c3j)
      small=1.d-06
      cof=a-b+z
      cof=dabs(cof)+small
      icof=int(cof)
      cof=1.d0-2.d0*mod(icof,2)
      cofcg=c3j*cof*dsqrt(2.d0*c+1.d0)
      return
      end  function cofcg


      SUBROUTINE TROISJ(XJ1,XJ2,XJ3,XM1,XM2,XM3,C3J)
      implicit double precision (A-H,O-Z)
      DIMENSION FLOG(301)                                               
C                                                                       
C     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS                       
C                                                                       
      DATA(FLOG(I),I=2,31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,
     A4.7874917D0,6.
     15792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0,
     B17.502307D0,19.98721
     24D0,22.552163D0,25.191221D0,27.899271D0,30.671860D0,33.505072D0,
     C36.395445D0,39.3
     339884D0,42.335616D0,45.380139D0,48.471180D0,51.606674D0,54.
     D784729D0,58.003604D0,
     461.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/   
      DATA(FLOG(I),I=32,61)/78.092223D0,81.557959D0,85.054466D0,88.
     A580827D0,92.1
     136175D0,95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.
     B32064D0,114.03421D0,
     2117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,
     C136.80272D0,140.67
     3392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.3311
     D2D0,164.32011D0,16
     48.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0,
     E188.62817D0/        
      DATA(FLOG(I),I=62,91)/192.73904D0,196.86618D0,201.00931D0,
     A205.16820D0,209.
     134258D0,213.53224D0,217.73693D0,221.95644D0,226.19054D0,
     B230.43904D0,234.70172D0,
     2238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,
     C260.56494D0,264.92
     3164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.
     D89313D0,291.32394D0,29
     45.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0,
     E318.15264D0/        
      DATA(FLOG(I),I=92,121)/322.66349D0,327.18529D0,331.71788D0,
     A336.26118D0,340
     1.81505D0,345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.
     B73937D0,368.35449
     2D0,372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,
     C396.24881D0,400.9
     33094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819
     DD0,429.21439D0,4
     433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,
     E457.81238D0/       
      DATA(FLOG(I),I=122,151)/462.60817D0,467.41220D0,472.22438D0,477.
     A04466D0,48
     11.87298D0,486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.
     B13282D0,511.0080
     22D0,515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,
     C540.41692D0,545.
     334417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772
     DD0,575.05753D0,
     4580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,
     E605.02010D0/      
      DATA(FLOG(I),I=152,181)/610.03738D0,615.06126D0,620.09170D0,
     A625.12866D0,63
     10.17208D0,635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.
     B48486D0,660.5662
     26D0,665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,
     C691.18340D0,696.
     330735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.
     D01551D0,727.17456D0,
     4732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,
     E758.24811D0/      
      DATA(FLOG(I),I=182,211)/763.44661D0,768.65061D0,773.86010D0,
     A779.07503D0,78
     14.29539D0,789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.
     B47747D0,815.7297
     23D0,820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,
     C847.35209D0,852.
     364036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.
     D15676D0,884.47488D0,
     4889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,
     E916.48547D0/      
      DATA(FLOG(I),I=212,241)/921.83732D0,927.19391D0,932.55521D0,
     A937.92118D0,94
     13.29181D0,948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.
     B21419D0,975.6123
     25D0,981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,
     C1008.0954D0,1013
     3.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.
     D7367D0,1046.1920D0,
     41051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,
     E1079.0129D0/      
      DATA(FLOG(I),I=242,271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.
     A9768D0,11
     106.4781D0,1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.
     B0452D0,1139.570
     26D0,1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,
     C1172.8063D0,1178
     3.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.
     D1818D0,1211.7577D0,
     41217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,
     E1245.2924D0/      
      DATA(FLOG(I),I=272,301)/1250.8944D0,1256.5003D0,1262.1097D0,
     A1267.7228D0,12
     173.3396D0,1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,
     B1301.4777D0,1307.116
     20D0,1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,
     C1341.0203D0,1346
     3.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,
     D1375.0499D0,1380.7334D0,
     41386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,
     E1414.9058D0/      
      DATA EPS1,EPS2/.1D0,-.2D0/                                       
C                                                                      
C     CALCUL DES COMBINAISONS J,M                                      
C     
c     write(21,*) 'Inside TROISJ'                                       
c     write(21,*) xj1,xj2,xj3,xm1,xm2,xm3,c3j
      XN=XJ2-XM2+EPS1                                                 
      IF(XN.LT.0.D0) XN=XN+EPS2                                       
      N1=XN                                                           
      xiii=dabs(xj1)+dabs(xj2)+dabs(xj3)
      if(xiii.lt.200) go to 9999
      write (99,*) 'xj1 xj2 xj3=', xj1,xj2,xj3
      stop 
 9999 continue
      XN=XJ3+XM3+EPS1                                                 
      IF(XN.LT.0.D0) XN=XN+EPS2                                      
      N2=XN                                                         
      IF(XN.LT.0.D0) XN=XN+EPS2                                    
      XN=XJ3-XM3+EPS1                                             
      N3=XN                                                      
      XN=XJ1+XM1+EPS1                                           
      IF(XN.LT.0.D0) XN=XN+EPS2                                
      N4=XN                                                   
      XN=XJ2+XJ3-XJ1+EPS1                                    
      IF(XN.LT.0.D0) XN=XN+EPS2                             
      N5=XN                                                
      XN=XJ1+XJ3-XJ2+EPS1                                 
      IF(XN.LT.0.D0) XN=XN+EPS2                          
      N6=XN                                             
      XN=XJ1+XJ2-XJ3+EPS1                              
      IF(XN.LT.0.D0) XN=XN+EPS2                       
      N7=XN                                          
      XN=XJ1-XM1+EPS1                               
      IF(XN.LT.0.D0) XN=XN+EPS2                    
      N8=XN                                       
      XN=XJ2+XM2+EPS1                            
      IF(XN.LT.0.D0) XN=XN+EPS2                 
      N9=XN                                    
      XN=XJ3-XJ1-XM2+EPS1                                       
      IF(XN.LT.0.D0) XN=XN+EPS2                                 
      N10=XN                                                   
      XN=XJ3-XJ2+XM1+EPS1                                      
      IF(XN.LT.0.D0) XN=XN+EPS2                               
      N11=XN                                                 
      XN=XJ1+XJ2+XJ3+EPS1                                  
      IF(XN.LT.0.D0) XN=XN+EPS2                           
      N12=XN                                             
      N12=N12+1                                         
C                                                                      
C     TESTS SUR LES J ET M                                            
C                                                                    
      K=N4*N8                                                       
      IF(K.LT.0) GO TO 59                                          
      K=N1*N9                                                     
      IF(K.LT.0) GO TO 59                                        
      K=N2*N3                                                   
      IF(K.LT.0) GO TO 59                                    
      IF(N5.LT.0) GO TO 59                                     
      IF(N6.LT.0) GO TO 59                                    
      IF(N7.LT.0) GO TO 59                                  
      L=N1-N2+N3-N4+N8-N9                                  
      IF(L.NE.0) GO TO 59                                           
      K=N12-1                                                      
      IF(K.GT.0) GO TO 60                                         
      C3J=1.D0                                                   
      RETURN                                                    
   59 C3J=0.D0                                                 
      RETURN                                                  
C                                                            
C     CALCUL DE LA SOMME ALTERNEE                           
C                                                          
   60 K=0                                                 
      L=-N10                                             
      IF(L.GT.K) K=L                                    
      L=-N11                                           
      IF(L.GT.K) K=L                                  
      L=N7                                           
      IF(N8.GT.L) L=N8                              
      IF(N9.GT.L) L=N9                             
      F=1.D0                                      
      S=1.D0                                                         
      I=K+1                                                         
   62 IF(I.GT.L) GO TO 80                                          
      IM1=I-1                                                     
      NN=(N7-IM1)*(N8-IM1)*(N9-IM1)                              
      ND=I*(N10+I)*(N11+I)                                      
      F=-F*FLOAT(NN)/FLOAT(ND)                               
      S=S+F                                                   
      I=I+1                                                  
      GO TO 62                                              
C                                                          
C     CALCUL DE LA RACINE                                 
C                                                        
   80 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N3+1)+FLOG(N4+1)+FLOG(N5+1)+FLOG(N6
     1+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N9+1)                              
      C2N=.5D0*C2N                                                      
      KM1=K-1                                                           
      KP1=K+1                                                           
      C2D=FLOG(KP1)+FLOG(N7-KM1)+FLOG(N8-KM1)+FLOG(N9-KM1)+FLOG(N10+KP1)
     1+FLOG(N11+KP1)+.5D0*FLOG(N12+1)                                  
C                                                                       
C     CALCUL DU C3J SANS PHASE                                          
C                                                                       
      F=C2D-C2N                                                         
      IF(F.GT.80.D0) GO TO 98                                           
      F=C2N/C2D                                                         
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 98                      
      C3J=S*DEXP(C2N-C2D)                                               
      GO TO 106                                                         
   98 IF(S) 100,59,102                                                  
  100 S=DLOG(-S)                                                        
      C3J=-DEXP(S+C2N-C2D)                                              
      GO TO 106                                                         
  102 S=DLOG(S)                                                         
      C3J=DEXP(S+C2N-C2D)                                               
C                                                                       
C     CALCUL DE LA PHASE                                                
C                                                                       
  106 XN=XJ1-XJ2-XM3+EPS1                                               
      IF(XN.LT.0.D0) XN=XN+EPS2                                         
      L=XN                                                              
      L=L+K                                                             
      K=L/2                                                             
      K=2*K                                                             
      IF(L.NE.K) C3J=-C3J                                               
      RETURN                                                            
      END   SUBROUTINE TROISJ



       SUBROUTINE SCHROD(N,L,J,NT,EMIN,E,U,V,HB,hMEN,HMEP,mesmn,mesmp)

       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NNP=10000)
       COMMON/DER/FCT(NNP),DF(NNP),H
       COMMON/BMAX/NMAXT
       COMMON/ENE/EMIN0,EMAX0,ECUT,NITER2
       DIMENSION U(NNP),V(NNP),HMEN(NNP),MESMN(NNP),
     1 MESMP(NNP),HMEP(NNP),DD(NNP)

c      HBAR=20.75D0
c      DD=H**2/HBAR


c      EMIN=EMIN0
       NODE=N-1
       EMAX=ECUT+EMAX0
       WLIM=1E-7

  56   CONTINUE
  

       DO 60 KT=1,NITER2
         
         ETRIAL = (EMIN+EMAX)/2.0
         U(1)=H**(L+1)
         U(2)=(2.*H)**(L+1)
	 ND=0
	 IOUT=0
	 DO I = 2,NMAXT-1
	   R=H*I          
c	   SX=DD*(V(I)-ETRIAL)


c	   SX=DD(i)*(V(I)-ETRIAL)


           IF (NT.EQ.0) THEN
             SX=H**2*(V(I)-ETRIAL/HMEN(I))
           ELSE
             SX=H**2*(V(I)-ETRIAL/HMEP(I))
            END IF 


             U(I+1)=(2+SX)*U(I)-U(I-1)
	   IF(U(I+1)*U(I).LT.0.0d0) ND=ND+1
           IF(ABS(U(I+1)).LT.WLIM.AND.ETRIAL.LT.0) THEN
             IOUT=1
             IK=I
             DO IR= IK,NMAXT-1
               U(IR+1)=0.
             END DO
             GO TO 46
           END IF
         END DO       
  46     CONTINUE
         IF(ND.GT.NODE) GO TO 58
         EMIN=ETRIAL
         IF (IOUT.EQ.1) GO TO 61
         GO TO 60
  58     CONTINUE
         EMAX=ETRIAL
         IF(IOUT.EQ.1) GO TO 61
         
  60   END DO 
  61   CONTINUE

C                                           Cut-Off E-spectrum
       IF (ETRIAL.GT.ECUT) GO TO 91
       E=ETRIAL

C                                           Normalize Radial Wave

       SUM=0.0
       DO 81 I=1,NMAXT
         IF (ABS(U(I)).LT.1E-15) U(I)=0
	 SUM=SUM+U(I)**2*H
  81   CONTINUE    
       SUM=SQRT(SUM)

       ABSWF=0.
       DO 90 I = 1,NMAXT
         U(I)=U(I)/SUM
         IF(ETRIAL.LT.0.AND.ABS(U(I)).GT.ABSWF) THEN
           IABSWF=I
	   ABSWF=ABS(U(I))
         END IF
  90   CONTINUE
 
       IF (IABSWF.EQ.NMAXT.AND.ETRIAL.LT.0) THEN
         WRITE(*,*) ' ERROR IN WAVEFUNCTIONS FOR BOUND STATES'
         WRITE(*,*) ' N,L,J,E',NODE,L,J,ETRIAL
         WLIM=WLIM*10
         WRITE(*,*) ' TRYING WITH WLIM= ',WLIM
         IF (WLIM.GT.0.001) STOP
         GO TO 56
       END IF
 91    CONTINUE
       RETURN
       END   SUBROUTINE SCHROD    





