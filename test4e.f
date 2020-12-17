C**********************************************************************
C                                                                     *
C QRPA Program                                                        *        
C                                                                     *
C**********************************************************************
C                                                                     |
C   Nov. 2001: testing the charge-exchange                            |
C              see notes below                                        |
C----------------------------------------------------------------------
C
C   La subphme da meta` dell' edm SOLO nel caso Isovettoriale ???  
C
C   WRITTEN APRIL 1988 BY BENT LAURITZEN,                          
C   for sd shell nuclei,Wildenthal effective interaction.         
C                                                                  
C   Remaked  since June 1993 by Marco Bianchetti & Mariella Quaglia:  
C   - Any transition with arbitrary J, Parity and Isospin.         
C   - Extension for standard QRPA, isoscalar channel (April 1995).
C     Electric transition probability.
C     Tested against code GIANT for isoscalar 2+ and 3- states in Ca40.
C   - April 1996: added Transition densities for QRPA phonons
C   - March 1997: Phase convention II is assumed everywhere.
C   _ September 1997: added transition densities for isovector charge 
C     exchange pnQRPA phonons
C   - March 1998 (G. Colo`): (1) introduced ishiftrpa (it is an input 
C     parameter) to do the RPA limit; this is actually useless for non charge-
C     exchange modes and it could be removed AS THIS VERSION IS SUPPOSED TO
C     WORK FOR NON CHARGE-EXCHANGE MODES; (2) introduced ECFMIN, related to
C     the problem of two qp configurations with 1=2 (?); (3) try to go 
C     towards a fully consistent QRPA according to Waroquier's paper. 
C     Introduced i_t3 (q.v.). Changed to pn formalism, no isospin ! 
C                                            
C     Parameters are explained in subroutine  HF 
C                                               
C**********************************************************************
C     L'INTERAZIONE DI PAIRING E':
C     (VZERO/2)*(1-Psigma)*(1-RHO(R)/RHOC)**GAMMA      
C     LEGGE DA PROGRAMMA HF-BCS IL POTENZIALE RADIALE:
C     (VZERO/2)*(1-RHO(R)/RHOC)**GAMMA 
C     questo viene ANTISIMMETRIZZATO nella subroutine ppmatrixelements! 

      subroutine qrpa

      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      character msg*20
      dimension J1(100),J2(100),J3(100),J4(100),J0(100),I0(100)
      dimension app(ncf,ncf),aph(ncf,ncf),aph1(ncf,ncf),
     &          bpp(ncf,ncf),bph(ncf,ncf),Bph1(ncf,ncf)
      dimension a(ncf,ncf),b(ncf,ncf)
      dimension freq(ncf),egtp(ncf),egtn(ncf)
      dimension bgtp(ncf),bgtn(ncf),bel_is(ncf),bel_iv(ncf)
      dimension fracts0_is(ncf),fracts1_is(ncf)
      dimension fracts0_iv(ncf),fracts1_iv(ncf)
      dimension up(nsp),un(nsp),vp(nsp),vn(nsp),ep(nsp),en(nsp)
c     dimension Rho(ncf,nnn)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON /PAR/ GPP,GPH,DGPP,DGPH
      COMMON /PAR1/ NGPP,NGPH
      COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/ focc(nsp),delta(nsp),spe(nsp)
      common /hf/nmax,nocc,nunocc,norb
      common /hf1/del
      common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &              phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)
      common /matel1/ppme_so(ncf,ncf)
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1calfv(7),cbetv(7) 
      COMMON /DENS/ dt(nnn),dn(nnn),dp(nnn)
      common /bwf/ wf(nnn,nsp),dwf(nnn,nsp)
      common /biqi/ iqi(ncf)
      common /bcons/ icons,iovert
      common /bres_so/ w2,w2p,c_so,ri_so 
      common/bafa0/afa0
      common/bnosc/nosc
      common/pairing/vpair(nnn),vpairp(nnn)
      common/bmax/nmaxt
      common/b_radii/r2l(12)
      common/b_potential/icoul_d,icoul_e,iso,isj,icm1
      common/btime/time1,sumtime

c --- MAIN CODE : Subroutines are called in the wanted order -------------------

      write(2,*) '   '
      write(2,*) '>>> STARTING PROGRAM QRPA'
      write(2,*) '   '

      Call ReadCom
      Call Read_BCS_Orsay
C     ppp = 0.d0
C     pppx = 0.d0
C     jtspec=0
C     do iii=1,nocc
C      do jjj=1,nocc
C       if(iq(iii).eq.1.and.iq(jjj).eq.1)then
C        call subphme_Coulomb
C    1   (iii,jjj,jjj,iii,jtspec,nmax,del,sphme)
C        call subphme_Coulomb_exchange
C    1   (iii,jjj,jjj,iii,jtspec,nmax,del,sphmex)
C        fac=sqrt(float(lj(iii))+1.0)*sqrt(float(lj(jjj))+1.0)
C        ctr=sphme*fac
C        ctrx=sphmex*fac
C        ppp=ppp+ctr
C        pppx=pppx+ctrx
C        write(23,*) iii,jjj,sphmex,fac,pppx
C       end if
C      end do
C     end do
C     write(23,*) 'ecd=',0.5d0*ppp
C     write(23,*) 'ecx=',0.5d0*pppx
C     stop
      Call label
      If (Ibasis.gt.0) Call Spin_Mat
      MIbasis = Iabs(Ibasis)
      If (Iwild.eq.0) then                         
         Call Skyrme
         msg='Entering_PHME'
c        call timelog(.false.,msg,time1,sumtime,elap)
         Call Phmatrixelements                             
         msg='Exiting_PHME'
c        call timelog(.false.,msg,time1,sumtime,elap)
         if(bcsn.eq.0.and.bcsp.eq.0)go to 50
         if(Gpp.ne.0.or.ngpp.ne.1) then
         if(MIbasis.eq.2) Call Ppmatrixelements
         if(MIbasis.eq.1) Call Pandyainv                 
         if(MIbasis.eq.1.and.Ictr.ge.1) Call Pandyadir   
         endif
50       continue
      else if (Iwild.eq.1.and.Mibasis.eq.1) then           
         stop '>>> OBSOLETE OPTION'
      else if (Iwild.eq.2.and.Mibasis.eq.1) then      
         stop '>>> OBSOLETE OPTION'
      endif
      Call QP(App,Aph,Aph1,Bpp,Bph,Bph1)
c ------------------------------------------------------------------------------
c     Loops over Gph, Gpp have been changed many times: last version has Gph=Gpp
c     and finds by dicotomy the value such that E(lowest state) is close to zero
c     (ISGDR)

      open(unit=97,file='spur_freq.dat',status='unknown')
      Gpp_min = Gpp
      Gpp_max = Gpp+float(ngpp)*dgpp
c     Gph_min = Gph
c     Gph_max = Gph+float(ngph)*dgph
        DO 30 J=1,NGPP                                    
          Gpp = (Gpp_min + Gpp_max)/2.d0 
          Gph = Gpp
          write(97,*) J,Gpp,Gph
          CALL MATRIX(APP,APH,APH1,BPP,BPH,BPH1,A,B)       
          msg='Entering_RPA'
c         call timelog(.false.,msg,time1,sumtime,elap)
          CALL RPA(A,B,FREQ,Iflag)
          msg='Exiting_RPA'
c         call timelog(.false.,msg,time1,sumtime,elap)
          CALL EIGSRT1(FREQ,A,B,NV1,NCF)
          write(97,*) Freq(1),Iflag
          write(97,*) 'Iflag should be 1 if Freq(1) is negative'
          IF (Iflag.EQ.1) goto 400 
            if (Ibel.ge.1) then
              Call Electro(A,B,FREQ,BEL_is,BEL_iv,FRACTS0_is,
     &        FRACTS0_iv,FRACTS1_is,FRACTS1_iv)                
            endif
          if (Freq(1).lt.1e-5.and.Freq(1).gt.0.d0) stop
          Gpp_min = Gpp
          Gpp = (Gpp_min + Gpp_max)/2.d0
          go to 30
400       Gpp_max = Gpp
          Gpp = (Gpp_min + Gpp_max)/2.d0
30      CONTINUE

c     if(ngph.eq.1)goto 666
c     DO 20 I=1,NGPH                
c       Gph = (Gph_min + Gph_max)/2.d0
c         write(97,*) I,Gph
c         CALL MATRIX(APP,APH,APH1,BPP,BPH,BPH1,A,B)       
c         CALL RPA(A,B,FREQ,Iflag)
c         CALL EIGSRT1(FREQ,A,B,NV1,NCF)
c         write(97,*) Freq(1),Iflag
c         write(97,*) 'Iflag should be 1 if Freq(1) is negative'
c         IF (Iflag.EQ.1) goto 440 
c           if (Ibel.ge.1) then
c             Call Electro(A,B,FREQ,BEL_is,BEL_iv,FRACTS0_is,
c    &        FRACTS0_iv,FRACTS1_is,FRACTS1_iv)                
c           endif

c         if (Freq(1).lt.1e-5.and.Freq(1).gt.0.d0) stop
c         Gph_min = Gph
c         Gph = (Gph_min + Gph_max)/2.d0
c         go to 20
c440       Gph_max = Gph
c         Gph = (Gph_min + Gph_max)/2.d0
c20    CONTINUE

      write(2,10)
10    format(//,30('*'),' THE END ',36('*'))

      return
      END

C**********************************************************************

        SUBROUTINE READcom

C	READS PARAMETERS FROM COMMAND FILE 

        implicit double precision (a-h,o-z)
        include 'param.qrpa'
        dimension up(nsp),un(nsp),vp(nsp),vn(nsp)
        dimension spep(nsp),spen(nsp),ep(nsp),en(nsp)
        dimension spe1(nsp),spe2(nsp)
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
        COMMON /PAR/ GPP,GPH,DGPP,DGPH
        COMMON /PAR1/ NGPP,NGPH
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/ focc(nsp),delta(nsp),spe(nsp)
        common/hf/nmax,nocc,nunocc,norb
        common/hf1/del
        common/bafa0/afa0
        common/bnosc/nosc
        COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
        common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1  calfv(7),cbetv(7) 
        common/bshrp/ishiftrpa
        common /bcons/ icons,iovert
        common /bspur/ eta  
        common /bslf/ i_out_slf
        common/bmax/nmaxt
        common/b_radii/r2l(12)

        iprint = 0
        ictr = 0
c       irpa = 2
        write(2,*)
        write(2,*) 'irpa =',irpa
        write(2,*)
        iwild = 0
        iuv = 0
        ispe = 0
        ibcs = 2
        ifile = 0 
c       READ (*,*) Iprint,ICTR,IRPA,IWILD,IUV,ISPE,IBCS,IFILE
        icons = 0
        eta = 5.d0*r2l(2)/3.d0
        if(iovert.eq.2)then
         write(2,*)
         write(2,*) 'eta = ',eta
         write(2,*)
        end if
c       read (*,*) icons,iovert,eta <---
        ishiftrpa = 0
        i_out_slf = 1
c       read *,ishiftrpa,i_so,i_out_slf
        ibasis = -2
c       read *,IBASIS,ISPIN,iorb,isflip,IPAR,irad
c - removed ITOT
        ichex = 0 
        if(ibasis.gt.0)ichex=1
        nmax = nmaxt-2
c       afa0 = 0.d0
c       read *,nmax,afa0 
c       ecfmin = 0.d0
c       ecfmax = 1000.d0
        ephocut = 1000.d0
        s0cut = 0.d0
c       read *,ECFMIN,ECFMAX,ephocut,S0cut
        igraf = 0
        ibeta = 1
        ibel = 1
c       read *,Igraf,IBETA,IBEL
c       READ (*,*),Pgapn,Pgapp
c       READ (*,*),GPP,DGPP,NGPP,GPH,DGPH,NGPH
        write(2,20) Iprint,ICTR,IRPA,IWILD,IUV,ISPE,IBCS,IBEL,IFILE,
     &   IBASIS,ISPIN,IORB,IPAR,ECFMIN,ECFMAX,Ephocut,S0cut,
     &   Igraf,IBETA,
     &   Pgapn,Pgapp,gpp,dgpp,ngpp,gph,dgph,ngph
20      format(//,10('*'),' READING COMMAND FILE ',43('*'),//,
     &   5x,' IPRINT,ICTR,IRPA,IWILD,IUV,ISPE,IBCS,IBEL,IFILE = ',9i3,/,
     &   5X,' IBASIS,ISPIN,IORB,IPAR        = ',4I4/,
     &   5x,' Ecfmin,Ecfmax,Ephocut,S0cut   = ',4(2x,f6.1)/,
     &   5X,' Igraf,Ibeta,                  = ',2I4,/,
     &   5x,' Pgapn, Pgapp    = ',2f9.3,/,
     &   5x,' Gpp, DGpp, NGpp; Gph, DGph, NGph = ',2(2e12.5,i4))

ccc-------------EXPLANATION OF PARAMETERS --------------------------------------

c Iprint = 0  -> prints something
c        = 1  -> print matrix elements
c        = 2  -> print   "       "     complete evaluation   
c Ictr   = 0  -> controls nothing
c        = 1  -> Control of matrix elements recoupling
c        = 2  -> ALSO control of rpa matrix diagonalization
c Irpa   = 0  -> QRPA (don't use this value, it's confusing)
c        < 0  -> particle-hole basis for RPA & pnRPA
c        = -1 -> TDA  (no RPA B matrix)
c        = -2 -> RPA 
c        > 0  -> 2 quasi-particle basis for QRPA & pnQRPA 
c        = 1  -> QTDA (no QRPA B matrix)
c        = 2  -> QRPA 
c        = 3  -> HF-BCS
c Iwild  = 0  -> Skyrme interaction
c        = 1  -> Wildenthal sd-shell int. ALL J configurations are needed
c        = 2  -> Wildenthal+McGrory+Millener-Kurath int. for sd+fp shell
c                for Ar40. ALL J config. are needed 
c Iuv    = 0  -> BCS U,V,factors 
c        = 1  -> Input U, V factors for some nuclei from sub. UVdata
c Ispe   = 0  -> Hartree-Fock Single Particle Energies 
c        = 1  -> Input SPE from sub. SPEdata (before BCS) 
c        = 2  -> Input QPE from sub. SPEdata (after BCS)
c Ibcs   = 0  -> Normal  BCS procedure (see Bent)
c        = 1  -> Refined BCS procedure (see Civitarese et al.,J.Phys.G,17,1991)
c        = 2  -> BCS made by Orsay code (E. Khan) 
c Ifile  = 0  -> Doesn't write anything in files
c        = 1  -> Writes in files uv factors, phonon energies, transition density
c                matrix (different files according to Ispin parameter)
c Ibasis > 0  -> Proton-neutron states are built (for pn-QRPA)
c        < 0  -> Proton-proton and neutron-neutron states (for QRPA) are built
c        = 1  -> all  |p,n;J> states are built
c        = 2  -> only |p,n;J=Ispin> states are built
c        = 0  -> Nothing happens
c        =-1  -> all  |p,p';J> & |n,n';J> states are built
c        =-2  -> only |p,p';J=Ispin> & |n,n';J=Ispin> states are built
c Ispin  = J  -> select J=Ispin 2-body states 
c Iorb   = L  -> select L multipolarity for spin-flip multipole transitions
c        = 0  -> for GT or Fermi transitions 
c                [see Drozdz, Osterfeld et al., Pyhs. Lett. B189, 271 (1987)].
c Ipar   =+1  -> positive parity states are selected
c        = 0  -> any parity states are built
c        =-1  -> negative parity states are selected
c Ecf*    =    -> cutoff energies for the 2-body basis
c Ephocut =    -> phonon energy cutoff
c S0cut   =    -> S0 fraction lower limit
c Ibel   = 0  -> No Electromagnetic transitions
c        = 1  -> Calculate electromagnetic transition probabilities with
c                multipolarity J = Ispin (Standard QRPA, Ibasis < 0).
c Igraf  = 0  -> Nessun grafico!!
c        = -1 -> B(E,l) graph (energy and fraction of EWSR)
c        = 1  -> Qrpa strength histogram graph
c        = 2  ->  "   summed strength graph
c        = 3  -> Both graphs
c Ibeta  = -1 -> graphs for beta - strength
c        = +1 -> graphs for beta + strength
c i_t3   = 0  -> When calculating a phme, the t3 contribution is excluded
c                and the part t0...t2 is treated as before
c                (needed for the pp: Pandya for t0...t2, t3 is treated 
c                separately)
c        = 1     t3 is included, therefore we have the whole interaction 
c                (needed for the ph: all interaction is treated on the 
c                same footing)

        RETURN
        END

C**********************************************************************

      SUBROUTINE Read_BCS_Orsay

      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common/hf/nmax,nocc,nunocc,norb
      common/hf1/del
      COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/bwg/wg(nnn,nmx),ewg(nmx),xor(nmx),dwg(nnn,nmx)
      common/bpot/xmn(nnn),xmp(nnn),vvn(nnn),vvp(nnn),yn(nnn),yp(nnn),
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)
      common/bafa0/afa0
      common/bnosc/nosc
      common/bnuc/anucl,znucl 
      common/bibcs/i_bcs_p,i_bcs_n
      common /bres_so/ w2,w2p,c_so,ri_so 
      common /bslf/ i_out_slf
      COMMON/DER/FCT(NNN),DF(NNN),H 
      common/pairing/vpair(nnn),vpairp(nnn) 
      dimension up(nsp),un(nsp),vp(nsp),vn(nsp),ep(nsp),en(nsp)
      dimension spep(nsp),spen(nsp),u(nsp),v(nsp) 
      dimension en1(nsp),ep1(nsp),pgp(nsp),pgn(nsp)
      DIMENSION DAL(3),MASH(3)
      dimension nqqx(20,2),nqqy(20,2)
      dimension taut(nnn),dt1(nnn),dt2(nnn)
      data hbc,mc2 /197.3d0,938.73d0/ 
103   FORMAT(2X,'NUMBER OF S.P. STATES I=',I3,' IS LARGER THAN',
     &' RESERVED DIMENSION NSP=',I3)
208   FORMAT(6E12.5) 
209   FORMAT(15I3)
210   format(6(1x,e12.5))

      pi=4.d0*datan(1.d0)

      write(2,10)
10    FORMAT(/,10('*'),' READING OUTPUT OF BCS-ORSAY (E. KHAN) ',
     &26('*'),/)
      write(2,11)
11    FORMAT(/,10('*'),
     &' >>> 20,21 must be linked to facoc.in, lecpot.dat '
     &,15('*'),/)

      open(unit=27,file='facoc.in',status='old')
      open(unit=10,file='fort.10',status='old')
      open(unit=50,file='occupied.dat',status='old')
      open(unit=51,file='unoccupied.dat',status='old')

      read(27,*) fp,fn,i_bcs_p,i_bcs_n
      if(i_bcs_p.eq.0.and.i_bcs_n.eq.0) go to 104 
      open(20,file='delta.dat',status='old')
      open(12,file='vpair.dat',status='old')
  104 read(10,*) nnmax
      READ(10,208) T0,T1,T2,T3,T13,X0
      READ(10,208) ALFE,X3,X1,X2,W2,W2P
      READ(10,208) ANUCL,ZNUCL,H,(DAL(I),I=1,3)
      iprot = int(znucl)
      ineut = int(anucl-znucl)
      del = h 
      READ(10,209)(MASH(I),I=1,3)
      READ(10,209) I1,I2,I3,I4,I5
      READ(10,208) (XMN(J),J=1,NNMAX),(XMP(J),J=1,NNMAX)
      READ(10,208) ( VVN(J),J=1,NNMAX),( VVP(J),J=1,NNMAX)
      READ(10,208) (YN(J),J=1,NNMAX),(YP(J),J=1,NNMAX)
      READ(10,208) (VSN(J),J=1,NNMAX),(VSP(J),J=1,NNMAX)
      READ(10,208) (VC(J),J=1,NNMAX)
      READ(10,208) (DT(J),J=1,NNMAX),(TAUT(J),J=1,NNMAX)
      READ(10,208) (DT1(J),J=1,NNMAX),(DT2(J),J=1,NNMAX)
      READ(10,208) (DP(J),J=1,NNMAX),(DN(J),J=1,NNMAX)
ccc   from routine lecpot; wpd, wnd not defined 

      if(i_bcs_p.eq.0.and.i_bcs_n.eq.0) go to 105 
      do i=1,nnmax
      read(12,*) ir,vpair(i)
cc potenziale per protone-protone
      vpairp(i)=vpair(i)
      enddo
  105 continue

c     read *,Nocc,nosc     

C     hbom = 41 * (anucl**(-0.33333333d0)) 
C     afa0=HBOM*MC2/(HBC**2)

c     nosc=0

      ehmax_p=-1000.0
      ihmax_p=-10
      epmin_p=1000.0
      ipmin_p=-10
      ehmax_n=-1000.0
      ihmax_n=-10
      epmin_n=1000.0
      ipmin_n=-10

      ibox=0
      if(dabs(afa0).lt.0.001)then
      ibox=1 
      write(2,*) '>>> BOX calculation '
      go to 2000
      end if 

      write(2,*) '>> H.O. calculation'
      write(2,*) '>>> afa0 = ',afa0
      write(2,*) '>>> nosc = ',nosc
      BOSC=1.0/afa0  
      BOSC=dsqrt(BOSC)   
 2000 open(unit=99,status='old',file='lines')
      rewind(99)
      DO 110 I=1,Nocc   
c     READ *, NN(I),LL(I),LJ(I),iq(I)
      read(99,*) NN(I),LL(I),LJ(I),iq(I)
      lev(i)=i
      noc=nn(i)
      lp=ll(i)
      jp=lj(i)
      lt1=iq(i)+1   
      delta(i)=0.d0
      if(i_bcs_p.ne.0.or.i_bcs_n.ne.0) read(20,*) delta(i)
      read(27,*) focc(i)  
      v(i)=dsqrt(focc(i))
      u(i)=dsqrt(1.0-focc(i))
      if(lt1.eq.2)vp(i)=v(i)
      if(lt1.eq.1)vn(i)=v(i)
      if(lt1.eq.2)up(i)=u(i)
      if(lt1.eq.1)un(i)=u(i)
      if(ibox.eq.0)CALL QWG_NONAG(bosc,nosc,lp,jp,lt1)
      if(ibox.eq.1)then
       read(50,*) spe(i),nn(i),ll(i),lj(i),iq(i)
       if(iq(i).eq.1)spep(i)=spe(i)
       if(iq(i).eq.0)spen(i)=spe(i)
       read(50,*) (wf(j,i),j=1,nmax)
       read(50,*)(dwf(j,i),j=1,nmax)
       go to 301
      end if
 
      spe(i) = ewg(noc)
      if(iq(i).eq.1)spep(i)=ewg(noc)
      if(iq(i).eq.0)spen(i)=ewg(noc)

      do 300 j=1,nmax
      dwf(j,i)=dwg(j,noc)
      wf(j,i)=wg(j,noc)

  300 continue
 
  301 if(i_bcs_p.eq.0.and.iq(i).eq.1)then
       if(spe(i).gt.ehmax_p)then 
        ehmax_p=spe(i)
        ihmax_p=i
       end if
      end if
      if(i_bcs_n.eq.0.and.iq(i).eq.0)then
       if(spe(i).gt.ehmax_n)then 
        ehmax_n=spe(i)
        ihmax_n=i
       end if
      end if

  110 continue  
      nocc1=nocc+1
      I=nocc 
  200 continue
c     READ *, LP,JP,LTP,N1,N2
      read(99,*,end=201) LP,JP,LTP,N1,N2    
      LTP1=LTP+1
      IBB=(JP-2*LP+3)/2
      NQQY(LP+1,IBB)=N1
      NQQX(LP+1,IBB)=N2 
      if(ibox.eq.0)CALL QWG_NONAG(BOSC,NOSC,LP,JP,LTP1)
      DO 202 ii= N1,N2
      I=I+1 
      IF(I.LE.NSP) GO TO 203
      WRITE (*,103) I,NSP
      STOP 
  203 continue  
      lev(i)=i 
      v(i)=0.d0
      u(i)=1.d0
      if(ltp1.eq.2)vp(i)=0.d0
      if(ltp1.eq.1)vn(i)=0.d0
      if(ltp1.eq.2)up(i)=1.d0
      if(ltp1.eq.1)un(i)=1.d0
      NN(I)=II  
      LL(I)=LP  
      LJ(I)=JP  
      iq(I)=LTP 
       
      write(2,*) nn(i),ll(i),lj(i),iq(i)

      if(ibox.eq.1)then
       read(51,*) spe(i),nn(i),ll(i),lj(i),iq(i)
       if(iq(i).eq.1)spep(i)=spe(i)
       if(iq(i).eq.0)spen(i)=spe(i)
       read(51,*) (wf(j,i),j=1,nmax)
       read(51,*)(dwf(j,i),j=1,nmax)
       go to 2040
      end if

      if(iq(i).eq.1)spep(i)=ewg(ii)
      if(iq(i).eq.0)spen(i)=ewg(ii)

      do 204 j=1,nmax
      dwf(j,i)=dwg(j,ii)
      wf(j,i)=wg(j,ii)
  204 continue

 2040 if(i_bcs_p.eq.0.and.iq(i).eq.1)then
c      write(2,*) i,epmin_p,spe(i)
       if(spe(i).lt.epmin_p)then 
        epmin_p=spe(i)
        ipmin_p=i
c       write(2,*) 'Changed',ipmin_p,epmin_p
       end if
      end if
      if(i_bcs_n.eq.0.and.iq(i).eq.0)then
       if(spe(i).lt.epmin_n)then 
        epmin_n=spe(i)
        ipmin_n=i
       end if
      end if
  202 CONTINUE  

      GO TO 200 
  201 continue

      if(i_bcs_p.eq.0.and.
     1abs(focc(ihmax_p)-1.0).gt.0.01)
     2epmin_p=ehmax_p
      if(i_bcs_n.eq.0.and.
     1abs(focc(ihmax_n)-1.0).gt.0.01)
     2epmin_n=ehmax_n
      write(2,*)
c     write(2,*) ihmax_p,ehmax_p,focc(ihmax_p)
c     write(2,*) epmin_p
c     write(2,*) ihmax_n,ehmax_n,focc(ihmax_n)
c     write(2,*) epmin_n

      Norb=i
      nunocc = norb - nocc
c     write(2,*) 'Exiting from implicit loops'
c     write(2,*) norb

      i_harmonic=0
      if(i_harmonic.eq.1)then
c ---- Test with H.O. wave functions -------------------------------------------
       b_hwf=1.d0
       do i=1,norb
       n_hwf=nn(i)-1
       l_hwf=ll(i)
       if(i.eq.1)then 
       write(94,*) i,nn(i),ll(i),lj(i),iq(i)
       write(94,*) n_hwf,l_hwf,b_hwf
       end if
       do j=1,nmax
       rdist=float(j)*del
       wf(j,i)=hwf(n_hwf,l_hwf,b_hwf,rdist)
       if(dabs(wf(j,i)).lt.1e-10)wf(j,i)=0.d0
       if(i.eq.1)write(94,*) rdist,wf(j,i)
       end do
       end do
c ------------------------------------------------------------------------------
      end if

      if(i_bcs_p.eq.0)fp=(epmin_p+ehmax_p)/2.0
      if(i_bcs_n.eq.0)fn=(epmin_n+ehmax_n)/2.0
      write(2,*) '>>> Proton and neutron Fermi energies:'
      write(2,*) '>>> fp,fn = ',fp,fn

      fp1 = fp
      fn1 = fn
   
      do i=1,norb
      pgapp = delta(i)*float(i_bcs_p)
      pgapn = delta(i)*float(i_bcs_n)
      if(iq(i).eq.1)Ep(i)=dsqrt((spep(i)-fp)**2+pgapp**2)
      if(iq(i).eq.0)En(i)=dsqrt((spen(i)-fn)**2+pgapn**2)
      end do

      box = float(nmax)*del           ! box dimension
      IF(NORB.GT.nsp) then            ! parameter nsp check
         write(2,20)         
20       format(/,10x,' WARNING : PARAMETER nsp 
     & ( = max No of single particle states) IS TOO SMALL ',/)
         STOP
      ELSE IF(NNMAX.GT.NNN)THEN         ! parameter nnn check
         write(2,30)
30       format(/,10x,' WARNING : PARAMETER nnn 
     & ( = max No.of points for radial quantities) IS TOO SMALL',/)
         STOP
      ENDIF

      if(i_out_slf.eq.1)then
      open(unit=98,status='unknown',file='slf.dat')
      write(98,*) nocc,norb,nmax,del
      end if

      write(2,52) Iuv,Ispe
52    format(10x,' QUASI PARTICLE STATES (IUV = ',i1,' ISPE = ',i1,
     &') : ',//,
     &  5x,' n  Iq ll lj  Focc    Spe     Qpe     Spe`    Qpe`
     &    Gap     V^2     U^2',/,
     &  5x,'--- -- -- --  ----   -----   -----   -----   -----
     &   -----   -----   -----',/,5X,'Neutron States ',/)
        
      iflag=0 
      do i=1,norb 
          if (Iq(i).eq.0) then ! neutron states
             if (spen(i).gt.fn1.and.iflag.eq.0) then
                write(2,111) fn1,fn
                iflag=1
             endif
             write(2,80) i,Iq(i),ll(I),lJ(I),focc(i),
     &       spen(I),En1(i),dabs(SPEn(i)-fn1),En(I),delta(I),
     &       Vn(I)**2,Un(i)**2
             if(i_out_slf.eq.1)then
             write(98,*) i,nn(i),ll(I),lJ(I),iq(i)
             write(98,*) spen(i),En(i)
             write(98,210) (wf(j,i),j=1,nmax)
             write(98,210) (dwf(j,i),j=1,nmax)
             end if 
          endif
          if (iflag.eq.0.and.i.eq.norb) write(2,111) fn1,fn
      enddo
      write(2,100)
      iflag=0
      do i=1,norb 
          if (Iq(i).eq.1) then ! proton states
             if (spep(i).gt.fp1.and.iflag.eq.0) then
                write(2,111) fp1,fp
                iflag=1
             endif
             write(2,80) i,Iq(i),ll(I),lJ(I),focc(i),
     &       spep(I),Ep1(i),dabs(SPEp(i)-fp1),Ep(I),delta(I),
     &       Vp(I)**2,Up(i)**2
             if(i_out_slf.eq.1)then
             write(98,*) i,nn(i),ll(I),lJ(I),iq(i)
             write(98,*) spep(i),Ep(i)
             write(98,210) (wf(j,i),j=1,nmax)
             write(98,210) (dwf(j,i),j=1,nmax)
             end if 
          endif
          if (iflag.eq.0.and.i.eq.norb) write(2,111) fp1,fp
      enddo
80    format(5x,i3,3i2,'/2',2x,f5.2,7f8.3)
100   format(/,5x,'Proton States ',/)
111   format(5x,'--Fermi Surface--',2f8.3,1x,38('-'))

      RETURN
      END

C**********************************************************************

      SUBROUTINE LABEL

c     Set up basis for angular momentum coupled two-body states 
c     for (Q)RPA calculations.
c     if Ibasis>0 Isovector Skyrme channel is needed and
c        proton-neutron |p,n;J,py> states are built;
c     if Ibasis<0 Isoscalar Skyrme channel is needed and
c        proton-proton |p,p';J,py> & neutron-neutron |n,n';J,py> states
c        are built.
 
      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      dimension econf(ncf)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common/hf/nmax,nocc,nunocc,norb
      common/hf1/del
      common /biqi/ iqi(ncf)
      common/bibcs/i_bcs_p,i_bcs_n
      common /bslf/ i_out_slf
      COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
      common/bshrp/ishiftrpa
      common /bwf/ wf(nnn,nsp),dwf(nnn,nsp)
      dimension up(nsp),un(nsp),vp(nsp),vn(nsp),ep(nsp),en(nsp)

      open(unit=390,status='unknown',file='reduced_EM_sp.dat')
      write(390,*) iorb,isflip,ispin
c     IBasis>0 : |(np,lp,Jp ; np',lp',Jp')JJ,prty>,
c                |(nn,ln,Jn ; nn',ln',Jn')JJ,prty> ; nv states
c     IBasis<0 : |(np,lp,Jp ; nn,ln,Jn)JJ,prty> ; nv states
c        for sd shell     : NV = 28 (any J)
c                         : J=1+  -> Nv=7,   J=0+ -> Nv=3
c        for sd+pfg shell : Any J -> Nv = ?, J=1+ -> Nv=18, J=0+ -> Nv=8
c        for sp+sd  shell : Any J -> Nv = ?, J=1+ -> Nv=? , J=0+ -> Nv=?
  
      if (ibasis.lt.0) write(2,40) IBASIS,ISPIN,IORB,IPAR
40    format(/,10('*'),' TWO BODY PROTON-PROTON & NEUTRON-NEUTRON
     & BASIS STATE ',11('*'),//,
     &10x,'Ibasis =',i2,' Ispin =',i2,' Iorb =',i2,' Ipar =',i2,//
     &10X,'State lev1,iq1,l1,J1   lev2,iq2,l2,J2     J  P    Econf',/,
     &10x,'----- --------------   --------------     -  -   ------'/)
      if (ibasis.gt.0) write(2,41) IBASIS,ISPIN,IORB,IPAR
41    format(/,10('*'),' TWO BODY PROTON-NEUTRON BASIS STATE ',
     &12('*'),//,
     &10x,'Ibasis =',i2,' Ispin =',i2,' Iorb =',i2,' Ipar =',i2,//
     &10X,'State levp,iqp,lp,Jp   levn,iqn,ln,Jn     J  P   Econf',/,
     &10x,'----- --------------   --------------     -  -  -------'/)

      write(2,*)
     1'espmin,espmax,ecfmin,ecfmax=',espmin,espmax,ecfmin,ecfmax
     
      gsp = 5.585694713d0
      gsn = -3.82608545d0

c-------------------------------------------------------------------------------
      nv=0                        
      if(ichex.eq.1)then
       nv1=0
       nv2=0
      end if
      ecf_f_min=1000.d0
      ecf_f_max=0.d0
c---------- First state --------------------------------------------------------
      do 10 i1=1,norb                
         lev1  = lev(i1)
         iq1   = iq(i1)
         focc1 = focc(i1)
         nn1   = nn(i1)
         l1    = ll(i1)         
         j1    = lj(i1)
         iq1   = iq(i1)
         spe1  = spe(i1)
c---------- Second state -------------------------------------------------------
         do 20 i2=1,norb             
            lev2  = lev(i2)
            iq2   = iq(i2)
            focc2 = focc(i2)
            nn2   = nn(i2)
            l2    = ll(i2)         
            j2    = lj(i2)
            iq2   = iq(i2)
            spe2  = spe(i2)
c---------- Before any selection rule, we calculate s.p. e.m. matrix el. -------
            
            ix2=lev1
            iy2=lev2
            Jx2=lj(ix2)
            Jy2=lj(iy2)
            lx2=ll(ix2)
            ly2=ll(iy2)

            cccc = 0.d0

            if(iq1.ne.iq2) go to 20307            

            lltry = lx2 + ly2 + ispin
            ipartry = 1-2*mod(lltry,2) 

            if (iq2.eq.1) then                              
              if(isflip.eq.0.and.ipartry.eq.1)then
               idummy=0
c the quantity in the following line is < ix2 || i^L O(EL) || iy2 > already in phase convention II
               Iexp = (-lx2+ly2+iorb)/2
               Isgn = Ipot(Iexp)
               rdint = 0.d0

                do ipoint=1,nmax
                rdist=float(ipoint)*del
                rdint=rdint+wf(ipoint,ix2)*wf(ipoint,iy2)
     &                *rdist**iorb*del                        !!! non ci va ispin??
c               write(390,*) rdist,wf(ipoint,ix2),wf(ipoint,iy2)
c               write(390,*) del,iorb,rdint
                enddo

               cccc = rdint*yl(lx2,jx2,ly2,jy2,iorb)*Isgn     !!! non ci va ispin??
c              write(390,*) 'prova',rdint,Isgn,yl(lx2,jx2,ly2,jy2,iorb)
               write(390,9701) 
     &         Iq(ix2),nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2)
     &         ,lj(iy2),
     &         idummy,Ispin,cccc
              end if
 9701         format(10(1x,i3),2x,e13.6) 
              if(isflip.eq.1.and.ipartry.eq.-1) then
               idummy=1
c the quantity is now < ix2 || i^L-1 O(ML) || iy2 > already in phase convention II
               iphase1 = (lj(iy2)-lj(ix2))/2 + Ispin  - 1
               iphase1 = 1-2*mod(abs(iphase1),2)
               iphase2 = (ll(iy2) - ll(ix2) + Ispin- 1) / 2
               iphase2 = 1-2*mod(abs(iphase2),2)
               fact1 = dsqrt(dfloat(lj(iy2)+1)*dfloat(2*Ispin+1))
     &         /dsqrt(16.d0*datan(1.d0))
               Rdint=0.d0
               do ipoints=1,nmax
                rdist=dfloat(ipoints)*del
                Rdint=Rdint+wf(ipoints,ix2)*wf(ipoints,iy2)
     &          *rdist**(Ispin-1)*del                           
               enddo
!               iphase3 = 1
!               if(lj(iy2).lt.2*ll(iy2))iphase3 = -1 
               iphase3 = 1-2*mod(abs(ll(iy2)+(1-lj(iy2))/2),2)
               iphase4 = 1-2*mod(abs((lj(ix2)+lj(iy2))/2-Ispin),2)
               squarebra = 1.d0 + 1.d0/2.d0/dfloat(Ispin)*iphase3
     &         *((lj(iy2)+1)+iphase4*(lj(ix2)+1))
               curl1 = (gsp - 2.d0/dfloat(Ispin+1))*dfloat(Ispin)/2.d0*
     &         cofcg(dfloat(lj(iy2))/2.d0,dfloat(Ispin),dfloat(lj(ix2))
     &         /2.d0,
     &         0.5d0,0.d0,0.5d0)*squarebra
               iphase5 = 1-2*mod(abs((lj(ix2)+lj(iy2))/2+Ispin),2)
               cg5 = cofcg(dfloat(lj(iy2))/2.d0,dfloat(Ispin-1),
     &         dfloat(lj(ix2))/2.d0,0.5d0,0.d0,0.5d0)
               call SIXJ(dfloat(lj(iy2))/2.d0,1.d0,dfloat(lj(iy2))/2.d0,
     &               dfloat(Ispin-1),dfloat(lj(ix2))/2.d0,dfloat(Ispin),
     &               csj2)
               curl2 = iphase5*2.d0/(dfloat(Ispin+1))
     &           *sqrt(
     &           dfloat(Ispin)
     &           *(2.d0*dfloat(Ispin)-1.d0)
     &           *(2.d0*dfloat(Ispin)+1.d0)
     &           /8.d0
     &           *dfloat(lj(iy2))
     &           *(dfloat(lj(iy2))+2.d0)
     &           *(2.d0*dfloat(lj(iy2))+2.d0))
     &           *cg5*csj2
!               write(*,*) curl2
               curl = curl1 + curl2 
               cccc = iphase1*iphase2*fact1*Rdint*curl
               write(390,9701) Iq(ix2),
     &         nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2),lj(iy2),
     &         idummy,Ispin,cccc
              end if
            else if (iq2.eq.0) then         ! messo i minuscola                 ! neutrons
              cccc = 0.d0
              if(isflip.eq.1.and.ipartry.eq.-1) then
               idummy=1
c the quantity is now < ix2 || i^L-1 O(ML) || iy2 > already in phase convention II
               iphase1 = (lj(iy2)-lj(ix2))/2 + Ispin  -1
               iphase1 = 1-2*mod(abs(iphase1),2)
               iphase2 = (ll(iy2) - ll(ix2) + Ispin - 1) / 2       
               iphase2 = 1-2*mod(abs(iphase2),2)
               fact1 = dsqrt(dfloat(lj(iy2)+1)*dfloat(2*Ispin+1))
     &         /dsqrt(16.d0*datan(1.d0))
               Rdint=0.d0
               do ipoints=1,nmax
                rdist=dfloat(ipoints)*del
                Rdint=Rdint+wf(ipoints,ix2)*wf(ipoints,iy2)
     &          *rdist**(Ispin-1)*del                        ! perchè ci metto Ispin-1????
               enddo
!               squarebra = 1.d0 + 1.d0/2.d0/dfloat(Ispin)*
!     &         (-1)**(dfloat(ll(iy2))+1.d0/2.d0-dfloat(lj(iy2)))
!     &         *((dfloat(lj(iy2))+1.d0)+
!     &         (-1)**(dfloat(lj(ix2))+dfloat(lj(iy2))-dfloat(Ispin))
!     &         *(dfloat(lj(ix2))+1.d0))
               iphase3 = 1-2*mod(abs(ll(iy2)+(1-lj(iy2))/2),2)
               iphase4 = 1-2*mod(abs((lj(ix2)+lj(iy2))/2-Ispin),2)     !  ----------> nonn ci va Ispin?????????
               squarebra = 1.d0 + 1.d0/2.d0/dfloat(Ispin)*iphase3
     &         *((lj(iy2)+1)+iphase4*(lj(ix2)+1))
               curl1 = gsn*dfloat(Ispin)/2.d0*
     &         cofcg(dfloat(lj(iy2))/2.d0,dfloat(Ispin),dfloat(lj(ix2))
     &         /2.d0,
     &         0.5d0,0.d0,0.5d0)*squarebra
!               iphase3 = 1                                         !  ----------> non capisco questa logica
!               if(lj(iy2).lt.2*ll(iy2))iphase3 = -1
!               iphase4 = 1-2*mod((lj(ix2)+lj(iy2))/2-Iorb,2)
               curl2 = 0.d0
               curl = curl1 + curl2
               cccc = iphase1*iphase2*fact1*Rdint*curl
               write(390,9701) Iq(ix2),
     &         nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2),lj(iy2),
     &         idummy,Ispin,cccc
              end if
            endif

20307       continue

c---------- Total charge -------------------------------------------------------
            iqtot = iq1 + iq2
c---------- Total angular momentum ---------------------------------------------
            jmin=iabs(j1-j2)/2       
            jmax=(j1+j2)/2           
            do 30 Jt=jmin,jmax    
c---------- Selection rules ----------------------------------------------------
            if (ibasis.lt.0) then               ! QRPA case
             if (iqtot.eq.2) ecf = Ep(i1) + Ep(i2) 
             if (iqtot.eq.0) ecf = En(i1) + En(i2) 
             if (iqtot.eq.2.and.Ep(i1).gt.(espmax-Fp).
     &       or.Ep(i2).gt.(espmax-Fp))
     &       go to 30 
             if (iqtot.eq.0.and.En(i1).gt.(espmax-Fn).
     &       or.En(i2).gt.(espmax-Fn))
     &       go to 30 
             if (ecf.lt.ecfmin)   goto 30
             if (ecf.gt.ecfmax)   goto 30      ! cutoff
             if (i_bcs_p.eq.0.and.iq1.eq.1.and.iq2.eq.1)then
              if(focc1.eq.1.d0.and.focc2.eq.1.d0)  goto 30    
              if(focc1.eq.0.d0.and.focc2.eq.0.d0)  goto 30   
             endif
             if (i_bcs_n.eq.0.and.iq1.eq.0.and.iq2.eq.0)then
              if(focc1.eq.1.d0.and.focc2.eq.1.d0)  goto 30 
              if(focc1.eq.0.d0.and.focc2.eq.0.d0)  goto 30 
             endif
             if (i_bcs_p.eq.1.and.iq1.eq.1.and.iq2.eq.1)then
              if (i1.gt.nocc.and.i2.gt.nocc) go to 30
             end if
             if (i_bcs_n.eq.1.and.iq1.eq.0.and.iq2.eq.0)then
              if (i1.gt.nocc.and.i2.gt.nocc) go to 30
             end if
            endif
            Iqtot = Iq1 + Iq2                 
            if (ibasis.lt.0) then            
              if (Iqtot.eq.1)    goto 30     
              if (lev2.lt.lev1)  goto 30     
            else if (ibasis.gt.0) then      
              if (iq1.ne.0)      goto 30   
              if (iq2.ne.1)      goto 30  
              isign = 1                          
              if(spe(i2).le.Fp)isign = -1 
c             ecf=En(i1)+Ep(i2)
              ecf=En(i1)+Ep(i2)+(Fp-Fn)*isign*ishiftrpa 
              if (ecf.lt.ecfmin) goto 30
              if (ecf.gt.ecfmax) goto 30      ! cutoff
              if (i_bcs_p.eq.0.and.i_bcs_n.eq.0) then
               if(focc1.eq.1.d0.and.focc2.eq.1.d0)  goto 30      
               if(focc1.eq.0.d0.and.focc2.eq.0.d0)  goto 30     
              endif
            endif
            modibasis=iabs(ibasis)
            if (modIbasis.eq.2) then       
               if (Jt.ne.Ispin) goto 30    
            endif 
            if (Ipar.eq.0) then           
               Iprty=0                    
            else if(Ipar.ne.0) then     
               Iprty=(-1)**(l1+l2)           
               if (Iprty.ne.Ipar) goto 30    
            endif
c---------- Building the configurations-----------------------------------------
               nv=nv+1                       
               if (nv.gt.ncf) then           
                  write(2,80) nv
                  stop
               endif
               ipp(nv)=i2
               inn(nv)=i1
               JJ(nv)=Jt                    
               iqi(nv)=-2*iq1+1
c              write(2,*) En(i1),Ep(i2)
c              write(2,*) (Fp-Fn),isign,ishiftrpa
c              write(2,*) ecf
               write(2,60) nv,ipp(nv),iq(i2),ll(i2),lJ(i2),
     &                        inn(nv),iq(i1),ll(i1),lJ(i1),
     &                        JJ(nv),Iprty,ecf
               if(ecf.lt.ecf_f_min)ecf_f_min=ecf
               if(ecf.gt.ecf_f_max)ecf_f_max=ecf
               if(i_out_slf.eq.1)
     &         write(98,*) nv,ipp(nv),inn(nv)
30          continue
20       continue
10    continue
60    format(10x,i5,3X,2(4i3,'/2',3x),2i3,f9.3)
80    format(//,10x,' WARNING : PARAMETER ncf 
     &( = max no. of two body states ) IS TOO SMALL ',/i7)

      write(2,70) Ispin,nv,nv*nv,nv*(nv+1)/2
70    format(/,10x,
     &     ' No.of two body J = ',i1,' basis states : nv   = ',
     &    i6,/,10x,' No.of matrix elements             : nv^2 = ',i6,/,
     &         10x,' No.of independent matrix elements :      = ',i6)

      if(i_out_slf.eq.1)then
      nterm = -1
      write(98,*) nterm,nterm,nterm
      end if

      return
      END

C**********************************************************************

      SUBROUTINE SPIN_MAT

C     REDUCED MATRIX ELEMENTS OF THE MULTIPOLE SPIN OPERATOR : 
C     <P||F||N>, F = r^L * [Y_L * sigma]_J=ispin. 
C     Function TJL is used.

      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      dimension spe1(nsp),spe2(nsp),spen(nsp),spep(nsp)
      COMMON /SPIN/ sigma(nsp,nsp)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common/hf/nmax,nocc,nunocc,norb
      common/hf1/del
      common /bwf/ wf(nnn,nsp),dwf(nnn,nsp)

      PI = 3.14159265358d0          ! Pi for Pitagora...
      write(2,60) ISPIN,IORB
60    format(//,10('*'),' REDUCED MATRIX ELEMENTS OF THE MULTIPOLE 
     &SPIN OPERATOR : ',11('*'),//,10x,'Ispin = ',i1,' Iorb = ',i1,/,
     &10x,' i, ip,iqp,Lp,Jp    in,iqn,Ln,Jn     Rad.Int.  Geom.part
     &   Mat. El.',/,
     &10x,'--- --------------  --------------   ---------  ---------
     &  ----------',/)
c-------------------------------------------------------------------------------
      do 10 i=1,nv
         ip=ipp(i)
         in=inn(i)
         iqp=iq(ip)
         iqn=iq(in)
         llp=ll(ip)
         lln=ll(in)
         lJp=lj(ip)
         lJn=lj(in)
c----- Selection rules ---------------------------------------------------------
         sigma(ip,in)=0.d0
c----- Calculation of reduced matrix element : Geometrical part ----------------
         if (Ispin.eq.0) then                       ! Fermi case
            if (lJp.ne.lJn) goto 10 
            xlJp = float(lJp)*0.5d0                ! true value
            Geom = dsqrt(2*xlJp+1)                  
         else if (Ispin.ne.0) then                  ! Multipole case
            TJLx = TJL(llp,lJp,lln,lJn,Ispin,Iorb)      
            if (Iorb.eq.0) TJLx=TJLx*dsqrt(4.d0*pi) 
            Iexp = (-llp+lln-Iorb)/2  
c  questa divisione di interi e' corretta perche' nei casi in cui tjl e' 
c  non nullo (che ovviamente sono gli unici interessanti) -llp+lln-Iorb 
c  e' pari (conservazione della parita' proveniente da el. di matrice 
c  ridotti di YL)
            TJLx = Ipot(Iexp)*TJLx ! Phase Convention II 
            Geom = TJLx                             ! Geometrical part
         endif
c ---- Radial integral ---------------------------------------------------------
         Rint = 0.d0    
         DO 20 k=1,Nmax
            rdist = float(k)*del                  ! radial distance
            xOp = rdist**Iorb                      ! radial operator 
            Rint = Rint + WF(k,ip)*Xop*WF(k,in)
20       enddo
         Rint = Del*Rint                          ! Radial Integral
         sigma(ip,in) = Rint*Geom                 ! Final m.e.
c -------- Printing the results -----------------------------------------------
         if (Iprint.ge.2) write(2,50) 
c         write(2,50) 
     &   i,ip,iqp,llp,lJp,in,iqn,lln,lJn,Rint,Geom,sigma(ip,in)
10    continue
50    format(10x,i3,1x,2(4i3,'/2',2x),3(f9.4,2x))
c---------- Explicit hand-made values for sd+pfg-shell -------------------------
ccc   sigma(1,1) = dsqrt(6.)            Fermi transitions 
ccc   sigma(2,2) = dsqrt(2.)            (others are zero)
ccc   sigma(3,3) = dsqrt(4.)            (Iorb=0,Ispin=0)
ccc   sigma(4,4) = dsqrt(8.)
ccc   sigma(5,5) = dsqrt(4.)
ccc   sigma(6,6) = dsqrt(6.)
ccc   sigma(7,7) = dsqrt(2.)
ccc   sigma(8,8) = dsqrt(10.)
ccc                                            Gamow-Teller transitions
ccc   sigma(1,1) = dsqrt(42./5.)        (others are zero)
ccc   sigma(1,3) =-dsqrt(48./5.)        (Iorb=0,Ispin=1)
ccc   sigma(2,2) = dsqrt(6.)
ccc   sigma(3,1) = dsqrt(48./5.)
ccc   sigma(3,3) =-dsqrt(12./5.)
ccc   sigma(4,4) = dsqrt(72./7.)
ccc   sigma(4,6) =-dsqrt(96./7.)
ccc   sigma(5,5) =-dsqrt(20./3.)
ccc   sigma(5,7) =-dsqrt(16./3.)
ccc   sigma(6,4) = dsqrt(96./7.)
ccc   sigma(6,6) =-dsqrt(30./7.)
ccc   sigma(7,5) = dsqrt(16./3.)
ccc   sigma(7,7) =-dsqrt(2./3.)
ccc   sigma(8,8) = dsqrt(110./9.)
c---- Independent particle Beta Strength ---------------------------------------
      SpBgtp = 0.d0
      SpBgtn = 0.d0
      DO 80 i=1,nv
         ip=ipp(i)
         in=inn(i)
         iqp=iq(ip)
         iqn=iq(in)
         llp=ll(ip)
         lln=ll(in)
         lJp=lj(ip)
         lJn=lj(in)
         Jt=JJ(i)           
         if (Jt.eq.ispin) then
         SpBgtp = SpBgtp + SIGMA(ip,in)**2*(1.d0-focc(ip))*focc(in)
         SpBgtn = SpBgtn + SIGMA(ip,in)**2*focc(ip)*(1.d0-focc(in))
         if (Iprint.ge.2) write(2,81) ip,in,sigma(ip,in)**2,
     &   SpBgtp,SpBgtn
         endif
80    continue
81    format(10x,2i4,3f10.4)
c --------- Gamow-Teller or Fermi Sum Rule -------------------------------------
      SUM = SpBgtp-SpBgtn
      if (Ispin.eq.0) RULE = float(Ineut-Iprot)        
      if (Ispin.eq.1) RULE = 3.d0*float(Ineut-Iprot)   
      if (Ispin.gt.1) RULE = 0.d0                       
      If (Rule.ne.0.) diffperc = dabs(Sum/Rule)*100.d0  

      RETURN
      END

C***********************************************************************

      SUBROUTINE SKYRME

c     Set up parameters to calculate Isoscalar or Isovector 
c     Skyrme particle-hole matrix elements in subphme

      implicit double precision (a-h,o-z)      
      include 'param.qrpa'
      common /hf/nmax,nocc,nunocc,norb
      common /hf1/del
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1calfv(7),cbetv(7) 
      COMMON/Dens/dt(nnn),dn(nnn),dp(nnn)
      common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)
      Common/Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &              iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common/Basis1/focc(nsp),delta(nsp),spe(nsp)
      common /bslf/ i_out_slf
      common /bres_so/ w2,w2p,c_so,ri_so 

      write(98,1) t0,t1,t2,t3
      write(98,1) x0,x1,x2,x3
      write(98,1) w2,w2p,alfe
      do j=1,nmax                               
      write(98,1) dp(j),dn(j),dt(j)
	end do
c     k_F
c     we choose k_F corresponding to half saturation density 
c     for all nuclei
      akf_test = 1.33d0
ccx   IF (ITOT.EQ.0) THEN           ! Isoscalar case (ITOT = 0)     
         f1(1)=0.75d0*t0       
         g1(1)=-0.25d0*t0*(1.d0-2.d0*x0)  
         do j=1,nmax                               
c           write(98,1) dp(j),dn(j),dt(j)
            f0(j,1)=f1(1)+(t3/48.D0)*(dt(j)**alfe)*(3.D0*(alfe+1.D0)*
     &      (alfe+2.D0)+alfe*(1.D0-alfe)*(1.D0+2.D0*x3)*
     &      (((dn(j)-dp(j))/dt(j))**2))  
            ctr_extra_pvc = 
     &      (akf_test**2/8.d0)*(3.d0*t1+t2*(5.d0+4.d0*x2))
c           if(i_out_slf.eq.1)write(98,1)f0(j,1),f0(j,1)+ctr_extra_pvc
C_TEST
c           f0(j,1)=f0(j,1)+ctr_extra_pvc   
            if(i_out_slf.eq.1)write(98,*) f0(j,1)
            g0(j,1)=g1(1)-(t3/24.D0)*(1.D0-2.D0*x3)*(dt(j)**alfe)
            corr_res(j)=-(t3/24.d0)*(2.d0*x3+1.d0)*alfe*(dt(j)**
     &      (alfe-1.d0))*(dn(j)-dp(j)) 
c           corr_res(j)=0.d0
         enddo
         avv1s=-3.d0*t1/32.d0
         avv2s=-2.d0*avv1s
         avv4s=-0.25d0*t2*(x2+5.d0/4.d0)
         bvv1s=-0.125d0*t1*(x1/2.d0-0.25d0)
         bvv2s=-2.d0*bvv1s
         bvv4s=-0.25d0*t2*(x2/2.d0+0.25d0)
ccx   ELSE if (Itot.eq.1) then                ! Isovector case (ITOT = 1)
         f1(2)=-0.25D0*t0*(1.D0+2.D0*x0)               
         g1(2)=-0.25D0*t0
         do j=1,nmax                                
            f0(j,2)=f1(2)-(t3/24.D0)*(1.D0+2.D0*x3)*(dt(j)**alfe)
            ctr_extra_pvc = 
     &      (akf_test**2/8.d0)*(-t1*(1.d0+2.d0*x1)+t2*(1.d0+2.d0*x2))
c           if(i_out_slf.eq.1)write(98,1)f0(j,2),f0(j,2)+ctr_extra_pvc
C_TEST
c           f0(j,2)=f0(j,2)+ctr_extra_pvc   
            if(i_out_slf.eq.1)write(98,*) f0(j,2)
            g0(j,2)=g1(2)-(t3/24.D0)*(dt(j)**alfe)
         enddo
         avv1v=(-t1/8.d0)*(-x1/2.d0-0.25d0)
         avv2v=-2.d0*avv1v
         avv4v=-0.25d0*t2*(x2/2.d0+0.25d0)
         bvv1v=-0.125d0*t1*(-0.25d0)        
         bvv2v=-2.d0*bvv1v
         bvv4v=-0.25d0*t2*0.25d0
ccx   ENDIF
      calfs(1)=avv1s
      calfs(2)=avv2s
      calfs(3)=avv2s
      cbets(1)=bvv1s
      cbets(2)=bvv2s
      cbets(3)=bvv2s
C_TEST
c     cbets(1)=0.d0
c     cbets(2)=0.d0
c     cbets(3)=0.d0
      calfv(1)=avv1v
      calfv(2)=avv2v
      calfv(3)=avv2v
      cbetv(1)=bvv1v
      cbetv(2)=bvv2v
      cbetv(3)=bvv2v
C_TEST
c     cbetv(1)=0.d0
c     cbetv(2)=0.d0
c     cbetv(3)=0.d0
      do iy=4,7
         iyy=1-2*(iy/6)
         calfs(iy)=iyy*avv4s
         calfv(iy)=iyy*avv4v
         cbets(iy)=iyy*bvv4s
         cbetv(iy)=iyy*bvv4v
C_TEST
c        cbets(iy)=0.d0
c        cbetv(iy)=0.d0
      enddo
      
    1 format(5(1x,e12.5))

      return
      end

C**********************************************************************

        SUBROUTINE PHMATRIXELEMENTS

C    PROBLEMA NORMALIZZAZIONE ANCORA APERTO !!!!!
c    Calculate proton-neutron particle-hole matrix elements
c    <p1,n1;J|V|p2,n2;J> of Skyrme effective nuclear interaction , 
c    coupled to total J with the scheme : p1-n1, p2-n2,
c    (1 = final states, 2 = initial states),
c    starting from single particle Hartree-Fock wave functions and densities

        include 'param.qrpa'
        implicit double precision (a-h,o-z)
        integer counter
c       real*4 sla,slb,slc,sld,sja,sjb,sjc,sjd,sjt,sj,c1,c2,c3,c4,
c    &         sl,slk1,slk2,sll,slka,slkb,slkc,slkd
        character msg*20
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
        common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &              phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)
        COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)           
        COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
        common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1  calfv(7),cbetv(7) 
        common/brint/ri0,ris,rx1,rx2(2,2)
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
        common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)
        common /biqi/ iqi(ncf)
        common /bres_so/ w2,w2p,c_so,ri_so 
        common/b_potential/icoul_d,icoul_e,iso,isj,icm1
        common/btime/time1,sumtime

C -------- MATRIX ELEMENTS CALCULATION------------------------------------------
c    Remember : from subphme we have : H(p1,n2,p2,n1;J),
c               J coupled with the scheme : p1-n1, p2-n2;
c               we need : <p1,n1|V|p2,n2> = - H(p1,n2,p2,n1;J) 
c               (1 = final states; 2 = initial states)
c    CAUTION : THAT MINUS SIGN WILL NOT BE USED

      write(2,10)
10    format(/,10('*'),' PARTICLE-HOLE MATRIX ELEMENTS ',34('*'),//,
     &'      LEV.   Iq  l   J       SPE    FOCC    JJ It
     &  MATRIX ELEMENT  COUNTER',/,
     &'    -------- -- --  ---    ------   ----    -- --
     & ---------------  -------',/)
  
      counter     = 0
      Icountminus = 0
      Icountplus  = 0
      Icountzero  = 0
      Sumt        = 0.d0
      Sumt2       = 0.d0
      Summinus    = 0.d0
      Sumplus     = 0.d0
      sum_time_skyrme = 0.d0
      sum_time_so = 0.d0
      DO 20 i1=1,nv                          ! final states
         ip1=ipp(i1)
         in1=inn(i1)
         JJ1=JJ(i1)
         It1=It(i1)
         Jp1=lj(ip1)
         Jn1=lj(in1)
         lp1=ll(ip1)
         ln1=ll(in1)
         Iprty1=(-1)**(lp1+ln1)
         Iq1=Iq(ip1)+Iq(in1)
         DO 30 i2=i1,nv                       ! initial states
            ip2=ipp(i2)
            in2=inn(i2)
            JJ2=JJ(i2)
            It2=It(i2)
            Jp2=lj(ip2)
            Jn2=lj(in2)
            lp2=ll(ip2)
            ln2=ll(in2)
            Iprty2=(-1)**(lp2+ln2)
            Iq2=Iq(ip2)+Iq(in2)
c ---------- Selection rules ---------------------------------------------------
            if (JJ1.ne.JJ2)                 goto 30 
            if (Iprty1.ne.iprty2)           goto 30 
            if (Ibasis.gt.0.and.Iq1.ne.Iq2) goto 30 
c ---------- Matrix element calculation ----------------------------------------
            Jt=JJ1                           ! Total angular momentum
            phme(i1,i2)=0.d0
            phme(i2,i1)=0.d0
            phmeex(i1,i2)=0.d0
            phmeex(i2,i1)=0.d0
            result = 0.d0
            resultex = 0.d0
            if(IRPA.eq.3)go to 1041
            iphase=iqi(i1)*iqi(i2)
            if(ichex.eq.1)iphase=1
            f1(3)=float(1-ichex)*f1(1) 
     1       + iphase*f1(2)
            g1(3)=float(1-ichex)*g1(1)
     1       + iphase*g1(2)
            do 1039 j=1,nmax
            f0(j,3)=float(1-ichex)*f0(j,1) 
     1             + iphase*f0(j,2) 
     1             - float(1-ichex)*float((1+iphase)/2)
     1              *float(iqi(i1))*2.d0*corr_res(j)
            g0(j,3)=float(1-ichex)*g0(j,1) 
     1       + iphase*g0(j,2)
 1039       continue
            do 1040 k=1,7
            calf(k)=float(1-ichex)*calfs(k) 
     1       + iphase*calfv(k)
            cbet(k)=float(1-ichex)*cbets(k) 
     1       + iphase*cbetv(k)
 1040       continue
C_phme of the Skyrme interaction (terms t0, t1, t2, t3)
            msg='After_initial'
c           call timelog(.false.,msg,time1,sumtime,elap)
            CALL subphme(1,Ip1,In2,Ip2,In1,Jt,nmax,del,result)  ! phme
            result0 = result
            msg='Skyrme_phme'
c           call timelog(.false.,msg,time1,sumtime,elap)
            sum_time_skyrme = sum_time_skyrme + elap 
C_Coulomb 
            result_Coulomb=0.d0
            result_Coulomb_ex=0.d0
            result_Coulomb_2=0.d0
            result_Coulomb_ex_2=0.d0
            if(ichex.eq.1)go to 1042
            if(iqi(i1).ne.-1.or.iqi(i2).ne.-1)go to 1042
            if(icoul_d.eq.1)then
             call subphme_Coulomb
     1       (Ip1,In2,Ip2,In1,Jt,nmax,del,result_Coulomb)
             msg='Coulomb_direct'
c            call timelog(.false.,msg,time1,sumtime,elap)
             if(icoul_e.eq.1) then 
              call subphme_Coulomb_exchange
     1        (Ip1,In2,Ip2,In1,Jt,nmax,del,result_Coulomb_ex)
              msg='Coulomb_exchange'
c             call timelog(.false.,msg,time1,sumtime,elap)
             end if
            end if
 1042       continue
            result_so=0.d0
            result_so_2=0.d0
C_Spin-orbit
            if(iso.eq.1)
     &      call subphme_so(JT,ip1,in2,ip2,in1,result_so)
            msg='Spin-orbit'
c           call timelog(.false.,msg,time1,sumtime,elap)
            sum_time_so = sum_time_so + elap 
c           call subphme_so(JT,ip2,in1,ip1,in2,result_symm)
c           diff_frac=dabs(result_so-result_symm)/dabs(result_so)
C_TEST
            result = result
     &      +result_Coulomb
     &      +result_Coulomb_ex
     &      +result_so
c           if(dabs(result_so).gt.0.05d0.and.diff_frac.gt.0.02)then
c           write(2,*)
c           write(2,*) i1,i2,ip1,in2,ip2,in1
c           write(2,*) 'Skyrme=',result0
c           write(2,*) 'spin-o=',result_so
c           write(2,*) 'Coulomb=',result_Coulomb
c           write(2,*) 'Coulomb-e=',result_Coulomb_ex
c           if(i1.eq.1.and.i2.eq.121)stop
c           stop
c           end if
            if (ibasis.lt.0) then                               ! QRPA
            CALL subphme(1,Ip1,Ip2,In2,In1,Jt,nmax,del,resultex)! PHMEex
             if(iso.eq.1)then
             call subphme_so(JT,ip1,ip2,in2,in1,result_so_2)
             end if
             if(icoul_d.eq.1)then
             call subphme_Coulomb
     1       (Ip1,Ip2,In2,In1,Jt,nmax,del,result_Coulomb_2)
             if(icoul_e.eq.1) then
              call subphme_Coulomb_exchange
     1        (Ip1,Ip2,In2,In1,Jt,nmax,del,result_Coulomb_ex_2)
             end if
             end if
            end if
C_TEST
            resultex = resultex
     &      +result_Coulomb_2
     &      +result_Coulomb_ex_2
     &      +result_so_2
c --------- Normalization ------------------------------------------------------
            Rnorm1 = 1.d0
            Rnorm2 = 1.d0
            if (ip1.eq.in1) Rnorm1=1.d0/dsqrt(2.d0)
            if (ip2.eq.in2) Rnorm2=1.d0/dsqrt(2.d0)
            result   = Rnorm1*Rnorm2*result
            resultex = Rnorm1*Rnorm2*resultex
c --------- Building particle-hole matrix --------------------------------------
 1041       continue
            phme(i1,i2) = result
            phme(i2,i1) = result
            phmeex(i1,i2) = resultex
            phmeex(i2,i1) = resultex           
c --------- Counting -----------------------------------------------------------
            counter = counter + 1                 ! counts m.e.
            Sumt = Sumt + phme(i1,i2)
            Sumt2 = Sumt2 + phme(i1,i2)**2
            if (phme(i1,i2).lt.0) then 
               Icountminus = Icountminus + 1      ! counts negative m.e.
               Summinus = Summinus + phme(i1,i2)  ! sums negative m.e.
            else if (phme(i1,i2).gt.0) then 
               Icountplus  = Icountplus  + 1      ! counts positive m.e.
               Sumplus  = Sumplus  + phme(i1,i2)  ! sums positive m.e.
            else if (phme(i1,i2).eq.0) then 
               Icountzero  = Icountzero  + 1      ! counts zero m.e.
            endif
C --------- Print the result----------------------------------------------------
      if (Iprint.ge.2) then
ccc   if (ip1.ne.ip2.or.in1.ne.in2) goto 30  ! print diagonal m.e. only
      write(2,50) Ip1,IQ(Ip1),Lp1,Jp1,spe(Ip1),FOCC(Ip1),
     &            In1,IQ(In1),Ln1,Jn1,spe(In1),FOCC(In1),
     &            Ip2,IQ(Ip2),Lp2,Jp2,spe(Ip2),FOCC(Ip2),
     &            In2,IQ(In2),Ln2,Jn2,spe(In2),FOCC(In2),
     &            JT,RESULT,resultex,COUNTER
50    FORMAT(4X,'p1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'p2 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n2 : ',4i3,'/2',f10.3,f7.2,3x,1i3,3x,f10.5,3x,
     &       f10.5,3x,i6,/)
      endif
c ---------- End of cycles -----------------------------------------------------
30       continue                               ! loop over Ip2,In2
20    continue                                  ! loop over Ip1,In1
      write(2,60) Icountminus,Icountplus,Icountzero,Counter,
     &            Summinus,Sumplus,Sumt
60    format(10x,' No. of negative m.e. = ',i11,/,
     &       10x,' No. of positive m.e. = ',i11,/,
     &       10x,' No. of zero     m.e. = ',i11,/,
     &       10x,' Total No. of    m.e. = ',i11,/
     &       10x,' Sum of negative m.e. = ',f11.5,/,
     &       10x,' Sum of positive m.e. = ',f11.5,/,
     &       10x,' Total Sum of    m.e. = ',f11.5)
      write(2,*) 
      write(2,*) 'Mean of phme: ',Sumt / counter
      write(2,*) 'R.m.s.: ',dsqrt(Sumt2/counter)
      ave_skyrme = sum_time_skyrme / counter
      ave_so = sum_time_so / counter
      write(2,*)
      write(2,*) ' Average time for t0...t3 and for spin-orbit : '
      write(2,*) ave_skyrme,ave_so

      return
      end

C**********************************************************************

        SUBROUTINE PPMATRIXELEMENTS

        include 'param.qrpa'
        implicit double precision (a-h,o-z)
        integer counter
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
        common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &                phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)

        common/bafa0/afa0
        COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)           
        COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
        common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1  calfv(7),cbetv(7) 
        common/brint/ri0,ris,rx1,rx2(2,2)
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
        common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)
        common /biqi/ iqi(ncf)
        common /bcons/ icons,iovert
        common /bres_so/ w2,w2p,c_so,ri_so 
        common/pairing/vpair(nnn),vpairp(nnn)
      common/b_potential/icoul_d,icoul_e,iso,isj,icm1
      write(2,10)
10    format(/,10('*'),' PARTICLE-PARTICLE MATRIX ELEMENTS ',30('*'),//,
     &'      LEV.   Iq  l   J       spe    FOCC    JJ It
     &  MATRIX ELEMENT  COUNTER',/,
     &'    -------- -- --  ---    ------   ----    -- --
     & ---------------  -------',/)

      if(icons.eq.0)i_t3_used = 1
      if(icons.eq.1)i_t3_used = 0

      counter     = 0
      Icountminus = 0
      Icountplus  = 0
      Icountzero  = 0
      Sum         = 0.d0
      Summinus    = 0.d0
      Sumplus     = 0.d0
      do 20 i1=1,nv
         ip1=ipp(i1)
         in1=inn(i1)
         Jp1=lj(ip1)
         Jn1=lj(in1)
         lp1=ll(ip1)
         ln1=ll(in1)
         Iprty1=(-1)**(lp1+ln1)   
         JJ1=JJ(i1)
         It1=It(i1)
         Iq1=Iq(ip1)+Iq(in1)
         sjp1=float(jp1)*0.5d0
         sjn1=float(jn1)*0.5d0
         do 30 i2=i1,nv         
            ip2=ipp(i2)
            in2=inn(i2)
            Jp2=lj(ip2)
            Jn2=lj(in2)
            lp2=ll(ip2)
            ln2=ll(in2)
            Iprty2=(-1)**(lp2+ln2)
            JJ2=JJ(i2)
            It2=It(i2)
            Iq2=Iq(ip2)+Iq(in2)
            sjp2=float(jp2)*0.5d0
            sjn2=float(jn2)*0.5d0
            if(IRPA.eq.3)go to 1042
c------Selection rules ---------------------------------------------------------
            if (JJ1.ne.JJ2)       goto 30       
            if (Iprty1.ne.Iprty2) goto 30       
            if (Ibasis.gt.0.and.It1.ne.It2) goto 30 
            if (Ibasis.gt.0.and.Iq1.ne.Iq2) goto 30 
            Jt1=JJ1
            sjt1=float(jt1)                             
c------p1-n2,p2,n1 pairs tot.ang.mom.coupling ----------------------------------
            Jp1n2max=(jp1+jn2)/2                
            Jp1n2min=iabs(jp1-jn2)/2            
            Iprtyp1n2=(-1)**(lp1+ln2)           
            Jp2n1max=(jp2+jn1)/2                
            Jp2n1min=iabs(jp2-jn1)/2            
c_new_Pandya
            Jp1p2max=(jp2+jp1)/2                
            Jp1p2min=iabs(jp2-jp1)/2            
            Jn1n2max=(jn2+jn1)/2                
            Jn1n2min=iabs(jn2-jn1)/2            
            Iprtyp2n1=(-1)**(lp2+ln1)           
            if (Iprtyp1n2.ne.Iprtyp2n1) goto 30 ! parity conservation
            Jt2min=max(Jp1p2min,Jn1n2min)       ! extrema of Jt
            Jt2max=min(Jp1p2max,Jn1n2max)       !    "
1042        continue
            ppme(i1,i2)=0.d0
            ppme(i2,i1)=0.d0
            if(IRPA.eq.3)go to 1041
            result=0.d0

c           dstep=0.1d0
C_GL_28_6_05
            dstep=del
 
c           goto 1043
c********* se gli stati sono neutronici il calcolo e' fatto con
c********* pairing density dependent interaction neutronica
            if(Iq(ip1).eq.0.and.Iq(in1).eq.0.and.
     &      Iq(ip2).eq.0.and.Iq(in2).eq.0) then
            
C************************************************
CCC  Parte radiale: wave function j wf(ir,j)
C************************************************ 
            sumrad=0.0d0
            do ir=1,nmax
            radp=dstep*ir
            sumrad=sumrad+dstep*VPAIR(ir)*wf(ir,ip1)*
     &      wf(ir,in1)*wf(ir,ip2)*wf(ir,in2)/radp/radp 
            enddo       
c           if(ip1.eq.12.and.in1.eq.12.and.ip2.eq.14.and.in2.eq.14)
c    &      write(98,*) 'R=',sumrad
C_GL_12_5_05
            lcount_min_1=iabs(jp1-jp2)/2
            lcount_max_1=(jp1+jp2)/2
            lcount_min_2=iabs(jn1-jn2)/2
            lcount_max_2=(jn1+jn2)/2
            lcount_min=max(lcount_min_1,lcount_min_2)
            lcount_max=min(lcount_max_1,lcount_max_2)
            ang=0.d0
            iphased=1-2*mod((jn1+jp2)/2+jt1,2)
c           if(ip1.eq.12.and.in1.eq.12.and.ip2.eq.14.and.in2.eq.14)
c    &      write(98,*) 'l1,l2=',lcount_min,lcount_max
            do lcount=lcount_min,lcount_max
            slam=float(lcount)
            call sixj(sjp1,sjp2,slam,sjn2,sjn1,sjt1,c1)  
            ridmat1=yl(lp1,jp1,lp2,jp2,lcount)
            ridmat2=yl(ln1,jn1,ln2,jn2,lcount)
            ang=ang+iphased*c1*ridmat1*ridmat2
c           if(ip1.eq.12.and.in1.eq.12.and.ip2.eq.14.and.in2.eq.14)
c    &      write(98,*) 'lcount,r1,r2,6j,iphase,ang',lcount,
c    &      ridmat1,ridmat2,c1,iphased,ang
            end do

            lcount_min_1=iabs(jp1-jn2)/2
            lcount_max_1=(jp1+jn2)/2
            lcount_min_2=iabs(jn1-jp2)/2
            lcount_max_2=(jn1+jp2)/2
            lcount_min=max(lcount_min_1,lcount_min_2)
            lcount_max=min(lcount_max_1,lcount_max_2)
            ang_ex=0.d0
            iphasex=1-2*mod((jp2+jn1)/2,2)
c           if(ip1.eq.12.and.in1.eq.12.and.ip2.eq.14.and.in2.eq.14)
c    &      write(98,*) 'l1,l2=',lcount_min,lcount_max
            do lcount=lcount_min,lcount_max
            slam=float(lcount)
            call sixj(sjp1,sjn2,slam,sjp2,sjn1,sjt1,c1)  
            ridmat1=yl(lp1,jp1,ln2,jn2,lcount)
            ridmat2=yl(ln1,jn1,lp2,jp2,lcount)
            ang_ex=ang_ex+iphasex*c1*ridmat1*ridmat2
c           if(ip1.eq.12.and.in1.eq.12.and.ip2.eq.14.and.in2.eq.14)
c    &      write(98,*) 'lcount,r1,r2,6j,iphase,ang_ex',lcount,
c    &      ridmat1,ridmat2,c1,iphasex,ang_ex
            end do

C********** normalizzazione *********************
            fac=1.0d0     
            if(ip1.eq.in1) fac=fac/dsqrt(2.0d0)
            if(ip2.eq.in2) fac=fac/dsqrt(2.0d0)

C************************************************
C           elemento di matrice            
C************************************************
            PPME(i1,i2)=fac*sumrad*(ang+ang_ex)
c           phase convention
            iphase_conv=iabs((-lp1-ln1+lp2+ln2)/2)
            iphase_conv=1-2*mod(iphase_conv,2)
            PPME(i1,i2)=PPME(i1,i2)*iphase_conv
c           write(98,98) ip1,lp1,jp1,spe(ip1),in1,ln1,jn1,spe(in1)
c           write(98,98) ip2,lp2,jp2,spe(ip2),in2,ln2,jn2,spe(in2)
c           write(98,*) 'ppme=',ppme(i1,i2)
c           write(98,*) 'factors:',fac,sumrad,ang,ang_ex,iphase_conv
98          format(1x,3(1x,i4),1x,f7.3,1x,3(1x,i4),1x,f7.3)

            goto 1041

            endif

            if(Iq(ip1).eq.1.and.Iq(in1).eq.1.and.
     &      Iq(ip2).eq.1.and.Iq(in2).eq.1) then

C************************************************
CCC  Parte radiale: wave function j wf(ir,j)
C************************************************ 
            sumrad=0.0d0
            do ir=1,nmax
            radp=dstep*ir
            sumrad=sumrad+dstep*VPAIRp(ir)*wf(ir,ip1)*
     &      wf(ir,in1)*wf(ir,ip2)*wf(ir,in2)/radp/radp 
            enddo       
C_GL_12_5_05
            lcount_min_1=iabs(jp1-jp2)/2
            lcount_max_1=(jp1+jp2)/2
            lcount_min_2=iabs(jn1-jn2)/2
            lcount_max_2=(jn1+jn2)/2
            lcount_min=max(lcount_min_1,lcount_min_2)
            lcount_max=min(lcount_max_1,lcount_max_2)
            ang=0.d0
            iphased=1-2*mod((jn1+jp2)/2+jt1,2)
c           if(ip1.eq.39.and.in1.eq.5.and.ip2.eq.39.and.in2.eq.5)
c    &      write(98,*) 'Limits of l=',lcount_min,lcount_max
            do lcount=lcount_min,lcount_max
            slam=float(lcount)
            call sixj(sjp1,sjp2,slam,sjn2,sjn1,sjt1,c1)  
            ridmat1=yl(lp1,jp1,lp2,jp2,lcount)
            ridmat2=yl(ln1,jn1,ln2,jn2,lcount)
            ang=ang+iphased*c1*ridmat1*ridmat2
c           if(ip1.eq.39.and.in1.eq.5.and.ip2.eq.39.and.in2.eq.5)
c    &      write(98,*) lcount,iphased,c1,ridmat1,ridmat2,ang
            end do

            lcount_min_1=iabs(jp1-jn2)/2
            lcount_max_1=(jp1+jn2)/2
            lcount_min_2=iabs(jn1-jp2)/2
            lcount_max_2=(jn1+jp2)/2
            lcount_min=max(lcount_min_1,lcount_min_2)
            lcount_max=min(lcount_max_1,lcount_max_2)
            ang_ex=0.d0
            iphasex=1-2*mod((jp2+jn1)/2,2)
c           if(ip1.eq.39.and.in1.eq.5.and.ip2.eq.39.and.in2.eq.5)
c    &      write(98,*) 'Limits of l=',lcount_min,lcount_max
            do lcount=lcount_min,lcount_max
            slam=float(lcount)
            call sixj(sjp1,sjn2,slam,sjp2,sjn1,sjt1,c1)  
            ridmat1=yl(lp1,jp1,ln2,jn2,lcount)
            ridmat2=yl(ln1,jn1,lp2,jp2,lcount)
            ang_ex=ang_ex+iphasex*c1*ridmat1*ridmat2
c           if(ip1.eq.39.and.in1.eq.5.and.ip2.eq.39.and.in2.eq.5)
c    &      write(98,*) lcount,iphasex,c1,ridmat1,ridmat2,ang_ex
            end do

C********** normalizzazione *********************
            fac=1.0d0     
            if(ip1.eq.in1) fac=fac/dsqrt(2.0d0)
            if(ip2.eq.in2) fac=fac/dsqrt(2.0d0)

C************************************************
C           elemento di matrice            
C************************************************
            PPME(i1,i2)=fac*sumrad*(ang+ang_ex)
c           phase convention
            iphase_conv=1-2*mod((-lp1-ln1+lp2+ln2)/2,2)
            PPME(i1,i2)=PPME(i1,i2)*iphase_conv
c           write(98,98) ip1,lp1,jp1,spe(ip1),in1,ln1,jn1,spe(in1)
c           write(98,98) ip2,lp2,jp2,spe(ip2),in2,ln2,jn2,spe(in2)
c           write(98,*) 'factors:',fac,sumrad,ang,ang_ex,iphase_conv
c           write(98,*) 'ppme=',ppme(i1,i2)

            goto 1041

            endif

c           iphase=iqi(i1)*iqi(i2)
ccx         if(ichex.eq.1)iphase=1
c           f1(3)=f1(1)
c    1       + iphase*f1(2)
c           g1(3)=g1(1)  
c    1       + iphase*g1(2)
c           do 1039 j=1,nmax
c           f0(j,3)=f0(j,1) 
c    1             + iphase*f0(j,2) -
c    1             float((1+iphase)/2)*float(iqi(i1))*
c    1             2.d0*corr_res(j)
c           g0(j,3)=g0(j,1) 
c    1       + iphase*g0(j,2)
c1039       continue
c           do 1040 k=1,7
c           calf(k)=calfs(k) 
c    1       + iphase*calfv(k)
c           cbet(k)=cbets(k)  
c    1       + iphase*cbetv(k)
c1040       continue
            do 40 kj=Jt2min,Jt2max       
               Jt2=kj                           
               sJt2=float(Jt2)                  
c              call subphme(0,Ip1,In1,In2,Ip2,Jt2,nmax,del,result)
c              If icons=0, i_t3_used=1 and all pp matrix element is treated 
c              with the Pandya relation; if icons=1, i_t3_used=0 and only
c              the t0...t2 part comes out here, the t3 part is added below
               call subphme(i_t3_used,Ip1,In1,In2,Ip2,Jt2,nmax,del
     1         ,result)
c              if(i_so.eq.1)then
c               call subphme_so(ip1,in1,in2,ip2,jt2,nmax,del,resultso)
c               jtprime_min_1=iabs(jp1-jn2)/2
c               jtprime_min_2=iabs(jp2-jn1)/2
c               jtprime_min=max(jtprime_min_1,jtprime_min_2)
c               jtprime_max_1=(jp1+jn2)/2
c               jtprime_max_2=(jp2+jn1)/2
c               jtprime_max=min(jtprime_max_1,jtprime_max_2)
c               resultso_ex=0.d0
c               do jtprime=jtprime_min,jtprime_max
c                call subphme_so(Ip1,In1,Ip2,In2,Jtprime,nmax,del,
c    1           resultso_ex_term)
c                iphase=(jp2+jn2)/2+Jt2+Jtprime
c                isgn=ipot(iphase)
c                sjtprime=float(jtprime)
c                call sixj(sjp2,sjp1,sjt2,sjn2,sjn1,sjtprime,c2)
c                resultso_ex=resultso_ex + isgn*float(2*Jtprime+1)*
c    1            c2*resultso_ex_term
c               end do
c               result=result+resultso-resultso_ex           
c              end if
c ---------- 6-J symbol calculation --------------------------------------------
               sJt1=float(Jt1)                  
               call sixj(sjp1,sjn1,sJt1,sjn2,sjp2,sJt2,c1)  ! 6j symbol 
               iphase_Pan = ipot ((jp2+jn2)/2 - Jt1)
c --------- Inverse Pandya relation --------------------------------------------
               ppme(i1,i2) = ppme(i1,i2) + iphase_Pan*     
     &                      (2.d0*sJt2+1.d0)*c1*result  
               if (iprint.ge.2) write(2,62) ip1,in2,ip2,in1,
     &                          Jt2,result,c1,ppme(i1,i2)
62             format(4x,5i3,3f8.3)
40          continue
            res_t3 = 0.d0
            if(icons.eq.1)call vppt3(ip1,in1,ip2,in2,iphase,jt1,
     &      nmax,del,res_t3)
c ADD THE PP MATRIX ELEMENT OF THE T3 PART CALCULATED DIRECTLY !             
c           call vppt3(ip1,in1,ip2,in2,iphase,jt1,nmax,del,res_t3) 
c           write(49,*) i1,i2,ppme(i1,i2),res_t3
            ppme(i1,i2) = ppme(i1,i2) + res_t3 
c --------- Normalization ------------------------------------------------------
            Rnorm1 = 1.d0
            Rnorm2 = 1.d0
            if (ip1.eq.in1) Rnorm1=1.d0/dsqrt(2.d0)
            if (ip2.eq.in2) Rnorm2=1.d0/dsqrt(2.d0)
            ppme(i1,i2)=ppme(i1,i2)*rnorm1*rnorm2 
c --------- Symmetry -----------------------------------------------------------
1041        continue
           
c            write(86,*)'********************************'
c            write(86,*)'config',i1,i2
c            write(86,*)'in1',ln1,jn1,lev(in1)
c            write(86,*)'ip1',lp1,jp1,lev(ip1)
c            write(86,*)'in2',ln2,jn2,lev(in2)
c            write(86,*)'ip2',lp2,jp2,lev(ip2)
c            write(86,*)'J,Q,PPME',JJ1,iq1,iq2,ppme(i1,i2)
c            write(86,*)sumrad,ridmat1,ridmat2
c            write(86,*)'*******************************'
         
            ppme(i2,i1)=ppme(i1,i2)
        

       
c --------- Counting -----------------------------------------------------------
            counter = counter + 1                 ! counts m.e.
            Sum = Sum + ppme(i1,i2)
            if (ppme(i1,i2).lt.0) then 
               Icountminus = Icountminus + 1      ! counts negative m.e.
               Summinus = Summinus + ppme(i1,i2)  ! sums negative m.e.
            else if (ppme(i1,i2).gt.0) then 
               Icountplus  = Icountplus  + 1      ! counts positive m.e.
               Sumplus  = Sumplus  + ppme(i1,i2)  ! sums positive m.e.
            else if (ppme(i1,i2).eq.0) then 
               Icountzero  = Icountzero  + 1      ! counts zero m.e.
            endif
C --------- Print the result----------------------------------------------------
      if (iprint.ge.2) then
ccc   if (ip1.ne.ip2.or.in1.ne.in2) goto 30  ! print diagonal m.e. only
      write(2,50) Ip1,Iq(Ip1),Lp1,jp1,Spe(Ip1),Focc(Ip1),
     &            In1,Iq(In1),Ln1,jn1,Spe(In1),Focc(In1),
     &            Ip2,Iq(Ip2),Lp2,jp2,Spe(Ip2),Focc(Ip2),
     &            In2,Iq(In2),Ln2,jn2,Spe(In2),Focc(In2),
     &            Jt1,PPME(I1,I2),Counter
50    FORMAT(4X,'p1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'p2 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n2 : ',4i3,'/2',f10.3,f7.2,3x,1i3,3x,f10.5,3x,i6,/)
      endif
c ---------- End of cycles -----------------------------------------------------
30       continue
20    continue
      write(2,60) Icountminus,Icountplus,Icountzero,Counter,
     &            Summinus,Sumplus,Sum
60    format(10x,' No. of negative m.e. = ',i11,/,
     &       10x,' No. of positive m.e. = ',i11,/,
     &       10x,' No. of zero     m.e. = ',i11,/,
     &       10x,' Total No. of    m.e. = ',i11,/
     &       10x,' Sum of negative m.e. = ',f11.5,/,
     &       10x,' Sum of positive m.e. = ',f11.5,/,
     &       10x,' Total Sum of    m.e. = ',f11.5)

      return
      end   

C**********************************************************************

        SUBROUTINE Pandyainv
 
c       Calculate proton-neutron particle-particle matrix elements
c       starting from the complete set (all J) of particle-hole matrix elements,
c       via inverse Pandya relation :
c       Ppme1(abcd)J = - sum{J'} {(2J'+1)sixj(adJ',cbJ)Phme(adcb)J'}
c       (be careful : subroutine sixj wants single precision variables)
c       Remeber the following scheme :
c       1) calculate ALL J configurations with any parity
c       2) calculate ALL J phme
c       3) call Pandyainv ---> you have Ppme1 (and Ppme for the following)
c       4) now call Pandyadir ---> you have Phme1 : they must coincide with Phme

        include 'param.qrpa'
        implicit double precision (a-h,o-z)
        integer counter
c       real*4 sjp3,sjn3,sjp4,sjn4,sJJ3,sJJ2,c1
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
        common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &              phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)
        COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)           
        COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
        common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1  calfv(7),cbetv(7) 
        common/brint/ri0,ris,rx1,rx2(2,2)
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
        common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)

      write(2,70)
70    format(//,10('*'),' PARTICLE-PARTICLE MATRIX ELEMENTS ',43('*')
     &,//,
     &'      LEV.   Iq  l   J       spe    FOCC    JJ It
     &  MATRIX ELEMENT  COUNTER',/,
     &'    -------- -- --  ---    ------   ----    -- --
     & ---------------  -------',/)

        counter     = 0
        Icountminus = 0
        Icountplus  = 0
        Icountzero  = 0
        Sum         = 0.d0
        Summinus    = 0.d0
        Sumplus     = 0.d0
        DO 10 i1=1,NV     
           ip1=ipp(i1)
           in1=inn(i1)
           Jp1=lJ(ip1)
           Jn1=lJ(in1)
           lp1=ll(ip1)
           ln1=ll(in1)
           Iprty1=(-1)**(lp1+ln1)
           JJ1=JJ(i1)
           It1=It(i1)
           Iq1=Iq(ip1)+Iq(in1)
           DO 20 i2=1,i1       
           ip2=ipp(i2)
           in2=inn(i2)
           Jp2=lJ(ip2)
           Jn2=lJ(in2)
           lp2=ll(ip2)
           ln2=ll(in2)
           Iprty2=(-1)**(lp2+ln2)
           JJ2=JJ(i2)
           It2=It(i2)
           Iq2=Iq(ip2)+Iq(in2)
           ppme1(i1,i2) = 0.d0
           ppme1(i2,i1) = 0.d0
           ppme (i1,i2) = 0.d0
           ppme (i2,i1) = 0.d0
c--------- Selectiion rules ----------------------------------------------------
           IF (JJ1.NE.JJ2)       GOTO 20          
           If (Iprty1.ne.Iprty2) goto 20          
           if (Ibasis.gt.0.and.It1.ne.It2) goto 20 
           if (Ibasis.gt.0.and.Iq1.ne.Iq2) goto 20 
           Jt=JJ1
c-------------------------------------------------------------------------------
           DO 30 i3=1,NV    
              ip3=ipp(i3)
              in3=inn(i3)
              Jp3=lJ(ip3)
              Jn3=lJ(in3)
              lp3=ll(ip3)
              ln3=ll(in3)
              Iprty3=(-1)**(lp3+ln3)
              JJ3=JJ(i3)
              It3=It(i3)
              Iq3=Iq(ip3)+Iq(in3)
              IF (ip3.NE.ip1) GOTO 30           ! p3=p1
              IF (in3.NE.in2) GOTO 30           ! n3=n2
              DO 40 i4=1,NV  
                 ip4=ipp(i4)
                 in4=inn(i4)
                 Jp4=lJ(ip4)
                 Jn4=lJ(in4)
                 lp4=ll(ip4)
                 ln4=ll(in4)
                 Iprty4=(-1)**(lp4+ln4)
                 JJ4=JJ(i4)
                 It4=It(i4)
                 Iq4=Iq(ip4)+Iq(in4)
                 IF (ip4.NE.ip2)       GOTO 40  
                 IF (in4.NE.in1)       GOTO 40 
c-------- Selection rules ------------------------------------------------------
                 IF (JJ3.NE.JJ4)       GOTO 40  
                 if (Iprty3.ne.Iprty4) goto 40  
                 if (Ibasis.gt.0.and.It3.ne.It4) goto 40 
                 if (Ibasis.gt.0.and.Iq3.ne.Iq4) goto 40 
c-------------------------------------------------------------------------------
                 sjp3=float(jp3)*0.5            
                 sjn3=float(jn3)*0.5            
                 sjp4=float(jp4)*0.5            
                 sjn4=float(jn4)*0.5            
                 sJJ2=float(JJ2)                
                 sJJ3=float(JJ3)                
                 call sixj(sjp3,sjn3,sJJ3,sjp4,sjn4,sJJ2,c1)
                 ppme1(i1,i2) = ppme1(i1,i2) - 
     &                          (2.d0*sJJ3+1.d0)*c1*PHME(i3,i4)
                 ppme(i1,i2)=ppme1(i1,i2)
                 If (Iprint.ge.2) write(2,50) jp1,jn1,jp2,jn2,JJ2,
     &                                        jp3,jn3,jp4,jn4,JJ3,
     &                                 c1,phme(i3,i4),ppme1(i1,i2)
50               format(4x,2(4(i2,'/2'),i3),3f8.3)
               GOTO 30
40            continue
30         continue
           ppme1(i2,i1) = ppme1(i1,i2)
           ppme (i2,i1) = ppme (i1,i2)
c --------- Counting -----------------------------------------------------------
            counter = counter + 1                 ! counts m.e.
            Sum = Sum + ppme1(i1,i2)
            if (ppme1(i1,i2).lt.0) then 
               Icountminus = Icountminus + 1      ! counts negative m.e.
               Summinus = Summinus + ppme1(i1,i2) ! sums negative m.e.
            else if (ppme1(i1,i2).gt.0) then 
               Icountplus  = Icountplus  + 1      ! counts positive m.e.
               Sumplus  = Sumplus  + ppme1(i1,i2) ! sums positive m.e.
            else if (ppme1(i1,i2).eq.0) then 
               Icountzero  = Icountzero  + 1      ! counts zero m.e.
            endif
C --------- Print the result----------------------------------------------------
      if (iprint.ge.1) then
ccc   if (ip1.ne.ip2.or.in1.ne.in2) goto 30  ! print diagonal m.e. only
      write(2,60) Ip1,Iq(Ip1),Ll(Ip1),Lj(Ip1),Spe(Ip1),Focc(Ip1),
     &            In1,Iq(In1),Ll(In1),Lj(In1),Spe(In1),Focc(In1),
     &            Ip2,Iq(Ip2),Ll(Ip2),Lj(Ip2),Spe(Ip2),Focc(Ip2),
     &            In2,Iq(In2),Ll(In2),Lj(In2),Spe(In2),Focc(In2),
     &            Jt,PPME1(I1,I2),Counter
60    FORMAT(4X,'p1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'p2 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n2 : ',4i3,'/2',f10.3,f7.2,3x,1i3,3x,f10.5,3x,i6,/)
      endif
c ---------- End of cycles -----------------------------------------------------
20      continue
10      continue
      write(2,61) Icountminus,Icountplus,Icountzero,Counter,
     &            Summinus,Sumplus,Sum
61    format(10x,' No. of negative m.e. = ',i11,/,
     &       10x,' No. of positive m.e. = ',i11,/,
     &       10x,' No. of zero     m.e. = ',i11,/,
     &       10x,' Total No. of    m.e. = ',i11,/
     &       10x,' Sum of negative m.e. = ',f11.5,/,
     &       10x,' Sum of positive m.e. = ',f11.5,/,
     &       10x,' Total Sum of    m.e. = ',f11.5)
      
        return
        END

C**********************************************************************

        SUBROUTINE Pandyadir
 
c       Calculate proton-neutron particle-hole matrix elements PHME1
c       starting from the whole set of particle-particle matrix elements PPME1,
c       via direct Pandya relation :
c       Phme1(abcd)J = - sum{J'} (2J'+1)sixj(adJ',cbJ)Ppme1(adcb)J'
c       (be careful : subroutine sixj wants single precision variables)
c       Remeber the following scheme :
c       1) calculate ALL J configurations
c       2) calculate ALL J phme
c       3) call Pandyainv ---> you have Ppme1 (and Ppme for the following)
c       4) now call Pandyadir ---> you have Phme1 : they must coincide with Phme

        include 'param.qrpa' 
        implicit double precision (a-h,o-z)
        integer counter
c       real*4 sjp3,sjn3,sjp4,sjn4,sJJ3,sJJ2,c1
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
        common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &              phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)
      write(2,70)
70    format(//,10('*'),' REEVALUATED PARTICLE-HOLE MATRIX ELEMENTS '
     &         ,25('*'),//,
     &'      LEV.   Iq  l   J       spe    FOCC    JJ It
     &  MATRIX ELEMENT  COUNTER',/,
     &'    -------- -- --  ---    ------   ----    -- --
     & ---------------  -------',/)

        counter     = 0
        Icountminus = 0
        Icountplus  = 0
        Icountzero  = 0
        Sum         = 0.d0
        Summinus    = 0.d0
        Sumplus     = 0.d0
        DO 10 i1=1,NV     
           ip1=ipp(i1)
           in1=inn(i1)
           Jp1=lJ(ip1)
           Jn1=lJ(in1)
           lp1=ll(ip1)
           ln1=ll(in1)
           Iprty1=(-1)**(lp1+ln1)
           JJ1=JJ(i1)
           It1=It(i1)
           Iq1=Iq(ip1)+Iq(in1)
        DO 20 i2=1,i1       
           ip2=ipp(i2)
           in2=inn(i2)
           Jp2=lJ(ip2)
           Jn2=lJ(in2)
           lp2=ll(ip2)
           ln2=ll(in2)
           Iprty2=(-1)**(lp2+ln2)
           JJ2=JJ(i2)
           It2=It(i2)
           Iq2=Iq(ip2)+Iq(in2)
           phme1(i1,i2) = 0.d0
           phme1(i2,i1) = 0.d0
           if (Iwild.eq.1.or.Iwild.eq.2) then
              phme(i1,i2) = 0.d0
              phme(i2,i1) = 0.d0
           endif
c--------- Selectiion rules ----------------------------------------------------
           IF (JJ1.NE.JJ2)       GOTO 20          
           If (Iprty1.ne.Iprty2) goto 20          
           if (Ibasis.gt.0.and.It1.ne.It2) goto 20 
           if (Ibasis.gt.0.and.Iq1.ne.Iq2) goto 20 
           Jt=JJ1
c-------------------------------------------------------------------------------
           DO 30 i3=1,NV    
              ip3=ipp(i3)
              in3=inn(i3)
              Jp3=lJ(ip3)
              Jn3=lJ(in3)
              lp3=ll(ip3)
              ln3=ll(in3)
              Iprty3=(-1)**(lp3+ln3)
              JJ3=JJ(i3)
              It3=It(i3)
              Iq3=Iq(ip3)+Iq(in3)
              IF (ip3.NE.ip1) GOTO 30       ! p3=p1
              IF (in3.NE.in2) GOTO 30       ! n3=n2
              DO 40 i4=1,NV  
                 ip4=ipp(i4)
                 in4=inn(i4)
                 Jp4=lJ(ip4)
                 Jn4=lJ(in4)
                 lp4=ll(ip4)
                 ln4=ll(in4)
                 Iprty4=(-1)**(lp4+ln4)
                 JJ4=JJ(i4)
                 It4=It(i4)
                 Iq4=Iq(ip4)+Iq(in4)
                 IF (ip4.NE.ip2) GOTO 40    
                 IF (in4.NE.in1) GOTO 40   
c-------- Selection rules ------------------------------------------------------
                 IF (JJ3.NE.JJ4)       GOTO 40  
                 if (Iprty3.ne.Iprty4) goto 40  
                 if (Ibasis.gt.0.and.It3.ne.It4) goto 40 
                 if (Ibasis.gt.0.and.Iq3.ne.Iq4) goto 40 
c-------------------------------------------------------------------------------
                 sjp3=float(jp3)*0.5        
                 sjn3=float(jn3)*0.5        
                 sjp4=float(jp4)*0.5        
                 sjn4=float(jn4)*0.5        
                 sJJ2=float(JJ2)            
                 sJJ3=float(JJ3)            
                 call sixj(sjp3,sjn3,sJJ3,sjp4,sjn4,sJJ2,c1)
                 phme1(i1,i2) = phme1(i1,i2) - 
     &                          (2.d0*sJJ3+1.d0)*c1*PPME1(i3,i4)
                 If (Iprint.ge.2) write(2,50) jp1,jn1,jp2,jn2,JJ2,
     &                                        jp3,jn3,jp4,jn4,JJ3,
     &                                 c1,ppme1(i3,i4),phme1(i1,i2)
50               format(4x,2(4(i2,'/2'),i3),3f8.3)
                GOTO 30
40            continue
30         continue
           phme1(i2,i1)=phme1(i1,i2)
           if (Iwild.eq.1.or.Iwild.eq.2) then
              phme(i1,i2)=phme1(i1,i2)
              phme(i2,i1)=phme1(i2,i1)
           endif
c --------- Counting -----------------------------------------------------------
            counter = counter + 1                 ! counts m.e.
            Sum = Sum + phme1(i1,i2)
            if (phme1(i1,i2).lt.0) then 
               Icountminus = Icountminus + 1      ! counts negative m.e.
               Summinus = Summinus + phme1(i1,i2) ! sums negative m.e.
            else if (phme1(i1,i2).gt.0) then 
               Icountplus  = Icountplus  + 1      ! counts positive m.e.
               Sumplus  = Sumplus  + phme1(i1,i2) ! sums positive m.e.
            else if (phme1(i1,i2).eq.0) then 
               Icountzero  = Icountzero  + 1      ! counts zero m.e.
            endif
C --------- Print the result----------------------------------------------------
      if (Iprint.ge.1) then
ccc   if (ip1.ne.ip2.or.in1.ne.in2) goto 30  ! print diagonal m.e. only
      write(2,60) Ip1,Iq(Ip1),Lp1,jp1,Spe(Ip1),Focc(Ip1),
     &            In1,Iq(In1),Ln1,jn1,Spe(In1),Focc(In1),
     &            Ip2,Iq(Ip2),Lp2,jp2,Spe(Ip2),Focc(Ip2),
     &            In2,Iq(In2),Ln2,jn2,Spe(In2),Focc(In2),
     &            Jt,PhME1(I1,I2),Counter
60    FORMAT(4X,'p1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n1 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'p2 : ',4i3,'/2',f10.3,f7.2,/,
     &       4X,'n2 : ',4i3,'/2',f10.3,f7.2,3x,1i3,3x,f10.5,3x,i6,/)
      endif
c ---------- End of cycles -----------------------------------------------------
20      continue
10      continue
      write(2,61) Icountminus,Icountplus,Icountzero,Counter,
     &            Summinus,Sumplus,Sum
61    format(10x,' No. of negative m.e. = ',i11,/,
     &       10x,' No. of positive m.e. = ',i11,/,
     &       10x,' No. of zero     m.e. = ',i11,/,
     &       10x,' Total No. of    m.e. = ',i11,/
     &       10x,' Sum of negative m.e. = ',f11.5,/,
     &       10x,' Sum of positive m.e. = ',f11.5,/,
     &       10x,' Total Sum of    m.e. = ',f11.5)

        return
        END

C***********************************************************************

      SUBROUTINE QP(APP,APH,APH1,BPP,BPH,BPH1)

C     Set up BCS factors for Qrpa A and B matrices
C     ready for pnQRPA, pnQTDA, QTDA, not for QRPA.

      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      integer counter
      dimension up(nsp),un(nsp),vp(nsp),vn(nsp),
     &          vn2(nsp),vp2(nsp),un2(nsp),up2(nsp),
     &          ep(nsp),en(nsp),en1(nsp),ep1(nsp)
      dimension app(ncf,ncf),aph(ncf,ncf),aph1(ncf,ncf),
     &          bpp(ncf,ncf),bph(ncf,ncf),Bph1(ncf,ncf)
      dimension u(nsp),v(nsp)
      COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
      COMMON /PAR/ GPP,GPH,DGPP,DGPH
      COMMON /PAR1/ NGPP,NGPH
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common /hf/nmax,nocc,nunocc,norb
      common /hf1/del

      if (iprint.ge.2) write(2,100)
100   format(/,10('*'),' U, V, factors for QRPA A and B Matrices ',
     &37('*'),//,2x,'  p  n  p` n` J  T          U,V factors      
     &        result  count',/,
     &            2x,' ------------------    --------------------- 
     &        ------  ----  ',//)
c---- ReSet U,V factors for isoscalar calculation ------------------------------
      if (ibasis.lt.0) then
        do k=1,norb
          if (iq(k).eq.1) then
            U(k)=Up(k)
            V(k)=Vp(k)
          else if(iq(k).eq.0) then
            U(k)=Un(k)
            V(k)=Vn(k)
          endif
        enddo
      endif
c---- Calculation of BCS factors -----------------------------------------------
      counter=0
      do 10 i1=1,nv
         ix1=ipp(i1)
         iy1=inn(i1)
         lx1=ll(ix1)
         ly1=ll(iy1)
         Iprty1=(-1)**(lx1+ly1)
         Jt1=JJ(i1)
         It1=It(i1)
         do 20 i2=1,i1
            ix2=ipp(i2)
            iy2=inn(i2)
            lx2=ll(ix2)
            ly2=ll(iy2)
            Iprty2=(-1)**(lx2+ly2)
            Jt2=JJ(i2)
            It2=It(i2)
c---------- Selectiion rules ---------------------------------------------------
            IF (Jt1.NE.Jt2)       GOTO 20          ! J conservation
            if (It1.ne.It2)       goto 20          ! T conservation
            If (Iprty1.ne.Iprty2) goto 20          ! parity conservation
c---- Construction of BCS factors ----------------------------------------------
      if(Ibasis.gt.0) then                               ! pnQRPA
         if(i1.eq.1.and.i2.eq.1)then 
          write(2,*) 'u,v factors: '
          write(2,*) ix1,up(ix1),vp(ix1)
         end if
         Aph(i1,i2) = Up(ix1)*Vn(iy1)*Up(ix2)*Vn(iy2)+   ! A, p-h
     &                Vp(ix1)*Un(iy1)*Vp(ix2)*Un(iy2)    
         App(i1,i2) = Up(ix1)*Un(iy1)*Up(ix2)*Un(iy2)+   ! A, p-p
     &                Vp(ix1)*Vn(iy1)*Vp(ix2)*Vn(iy2)    
         Bph(i1,i2) = Up(ix1)*Vn(iy1)*Vp(ix2)*Un(iy2)+   ! B, p-h
     &                Vp(ix1)*Un(iy1)*Up(ix2)*Vn(iy2)
         Bpp(i1,i2) = Up(ix1)*Un(iy1)*Vp(ix2)*Vn(iy2)+   ! B, p-p
     &                Vp(ix1)*Vn(iy1)*Up(ix2)*Un(iy2)
         Aph(i2,i1) = Aph(i1,i2)
         App(i2,i1) = App(i1,i2)
         Bph(i2,i1) = Bph(i1,i2)
         Bpp(i2,i1) = Bpp(i1,i2)
         counter = counter+1
         if (iprint.ge.2) write(2,30) Ix1,Iy1,Ix2,Iy2,Jt1,It1,
     &      Aph(i1,i2),App(i1,i2),Bph(i1,i2),Bpp(i1,i2),counter
      else if(ibasis.lt.0) then                           ! QRPA
         Aph(i1,i2) =  U(ix1)*V(iy1)*U(ix2)*V(iy2)+       ! A, p-h
     &                 V(ix1)*U(iy1)*V(ix2)*U(iy2)    
         Aph1(i1,i2)=  U(ix1)*V(iy1)*V(ix2)*U(iy2)+       ! A, p-h
     &                 V(ix1)*U(iy1)*U(ix2)*V(iy2)    
         App(i1,i2) =  U(ix1)*U(iy1)*U(ix2)*U(iy2)+       ! A, p-p
     &                 V(ix1)*V(iy1)*V(ix2)*V(iy2)    
         Bph(i1,i2) =  U(ix1)*V(iy1)*V(ix2)*U(iy2)+       ! B, p-h
     &                 V(ix1)*U(iy1)*U(ix2)*V(iy2)
         Bph1(i1,i2) = U(ix1)*V(iy1)*U(ix2)*V(iy2)+       ! B, p-h
     &                 V(ix1)*U(iy1)*V(ix2)*U(iy2)
         Bpp(i1,i2) =  U(ix1)*U(iy1)*V(ix2)*V(iy2)+       ! B, p-p
     &                 V(ix1)*V(iy1)*U(ix2)*U(iy2)
         Aph(i2,i1)  = Aph(i1,i2)
         Aph1(i2,i1) = Aph1(i1,i2)
         App(i2,i1)  = App(i1,i2)
         Bph(i2,i1)  = Bph(i1,i2)
         Bph1(i2,i1) = Bph1(i1,i2)
         Bpp(i2,i1)  = Bpp(i1,i2)
         counter = counter+1
         if (iprint.ge.2) write(2,31) Ix1,Iy1,Ix2,Iy2,Jt1,It1,
     &      Aph(i1,i2),Aph1(i1,i2),App(i1,i2),Bph(i1,i2),Bph1(i1,i2),
     &      Bpp(i1,i2),counter
      endif
c---- End of cycles ------------------------------------------------------------
20       continue
10    continue

30    format(2x,6i3,'  Aph = ',f6.3,14x,'  App = ',f6.3,/,
     &          20x,'  Bph = ',f6.3,14x,'  Bpp = ',f6.3,i6,/)
31    format(2x,6i3,'  Aph = ',f6.3,'  Aph1 = ',f6.3,'  App = ',f6.3,/,
     &          20x,'  Bph = ',f6.3,'  Bph1 = ',f6.3,'  Bpp = ',f6.3,i6
     &,/)


      RETURN
      END

C**********************************************************************


        SUBROUTINE MATRIX(APP,APH,APH1,BPP,BPH,BPH1,A1,B1)

C	SET UP A AND B RPA MATRICES 
c       and eventually reduce their dimension to the subspace of interest 

        include 'param.qrpa'
        implicit double precision (a-h,o-z)
        integer   counter
        dimension A(ncf,ncf),B(ncf,ncf),A1(ncf,ncf),B1(ncf,ncf)
        dimension up(nsp),un(nsp),vp(nsp),vn(nsp)
        dimension spep(nsp),spen(nsp),ep(nsp),en(nsp)
        dimension app(ncf,ncf),aph(ncf,ncf),aph1(ncf,ncf),
     &            bpp(ncf,ncf),bph(ncf,ncf),Bph1(ncf,ncf)
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
        COMMON /AB/ NUM(ncf)
        COMMON /PAR/ GPP,GPH,DGPP,DGPH
        COMMON /PAR1/ NGPP,NGPH
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /matel/phme(ncf,ncf),phmeex(ncf,ncf),ppme(ncf,ncf),
     &              phme1(ncf,ncf),ppme1(ncf,ncf),tbint(1000)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
        common/bshrp/ishiftrpa
        common /bcons/ icons,iovert
        common /bgfb/ ecf(ncf),akeep(ncf)
      
        if(ngph.eq.1.and.ngpp.eq.1) then
        write(2,*)
        write(2,*) 'Entering MATRIX: Fp - Fn = ',Fp-Fn
        write(2,*)
        endif

        if (Iprint.ge.1) write(2,100) Gph,Gpp
100     format(//,10('*'),' QRPA A and B MATRICES ',55('*'),//,
     &  10x,' Gph and Gpp multiplicative parameters : ',2e12.5,//,2x, 
     &  '  p  n  p` n` J  T            E    PHME   UVph   PPME   UVpp 
     & RPAmat count',/,2x,
     &  ' ------------------        ------ ------ ------ ------ ------
     & ------ -----',/)
c101    format(2x,i4,2x,i4,2x,i4,3x,3e12.5)
        counter=0
        modIRPA=abs(IRPA)
        DO 10 i1=1,NV
           ip1=ipp(i1)
           in1=inn(i1)
           lp1=ll(ip1)
           ln1=ll(in1)
           Iprty1=(-1)**(lp1+ln1)
           Jt1=JJ(i1)
           It1=It(i1)
           iq1=iq(Ip1)+iq(In1)
           DO 20 i2=i1,nv
              ip2=ipp(i2)
              in2=inn(i2)
              lp2=ll(ip2)
              ln2=ll(in2)
              jp2=lj(ip2)
              jn2=lj(in2)
              Iprty2=(-1)**(lp2+ln2)
              Jt2=JJ(i2)
              It2=It(i2)
              A(i1,i2) = 0.d0
              B(i1,i2) = 0.d0
c-------- Selection rules -----------------------------------------------------
              IF (Jt1.NE.Jt2)       GOTO 20          
              if (It1.ne.It2)       goto 20          
              If (Iprty1.ne.Iprty2) goto 20          
c-------- Diagonal Energy term ------------------------------------------------
              Econf = 0.d0
              IF (i1.EQ.i2) then                     
                 if (ibasis.lt.0) then               
                    if (iq1.eq.2) Econf = Ep(Ip1) + Ep(In1) 
                    if (iq1.eq.0) Econf = En(Ip1) + En(In1) 
                    write(2,*) 'conf,ih,ip,Econf: ',i1,in1,ip1,Econf
                 endif
                 if (ibasis.gt.0) then  
                  isign = 1                          
                  if(spe(Ip1).le.Fp)isign = -1 
                  Econf = Ep(Ip1) + En(In1) + (Fp - Fn)*isign*ishiftrpa 
cx                write(2,*) 'i1,E,E+shift=',i1,Ep(Ip1) + En(In1),
cx   &            Econf 
                 Endif
              Endif
c ------- A & B Matrices -------------------------------------------------------
      if (ibasis.gt.0) then                                    
c        write(2,*) i1,i2
c        write(2,*) Gph,Gpp
c        write(2,*) phme(i1,i2),aph(i1,i2)
c        write(2,*) bph(i1,i2)
         A(i1,i2) = Econf +     Gph * PHME(i1,i2)  * APH(i1,i2) +
     &                          Gpp * PPME(i1,i2)  * APP(i1,i2)
         B(i1,i2) =             Gph * PHME(i1,i2)  * BPH(i1,i2) -
     &                          Gpp * PPME(i1,i2)  * BPP(i1,i2)
         if(IRPA.eq.3)then
          A(i1,i2) = Econf
          B(i1,i2) = 0.d0
         end if
c        write(2,*) 'i1,i2,A,B=',i1,i2
c        write(2,*) A(i1,i2)-Econf,B(i1,i2)
c        write(2,*) PHME(i1,i2),BPH(i1,i2)
      else if (ibasis.lt.0) then                               ! QRPA
         sjp2=float(jp2)*0.5d0   ! true jp2
         sjn2=float(jn2)*0.5d0   ! true jn2
         jppjn=int(sjp2+sjn2)
         jpmjn=int(sjp2-sjn2)
          A(i1,i2) = Econf +     Gph * PHME(i1,i2)  * APH(i1,i2)
     &    -(-1)**(Jt1+jppjn)  * Gph * PHMEex(i1,i2)* APH1(i1,i2)
     &                        + Gpp * PPME(i1,i2)  * APP(i1,i2)
          B(i1,i2) =             Gph * PHME(i1,i2)  * BPH(i1,i2)
     &    +(-1)**(Jt1+jpmjn)  * Gph * PHMEex(i1,i2)* BPH1(i1,i2)
     &                        - Gpp * PPME(i1,i2)  * BPP(i1,i2)
      endif
      IF (modIRPA.EQ.1) B(i1,i2)=0.d0      ! (Q)TDA or pn(Q)TDA limit
      if(IRPA.eq.3)then
       A(i1,i2) = Econf
       B(i1,i2) = 0.d0
      end if 
      A(i2,i1)=A(i1,i2)
      B(i2,i1)=B(i1,i2)      
      if(i1.eq.i2) then 
       ecf(i1)=Econf
       akeep(i1)=A(i1,i1)-Econf
c      if(ngph.eq.1.and.ngpp.eq.1) then
       write(2,1010) ecf(i1),akeep(i1),B(i1,i2)
c      endif
 1010  format(1x,'Energy,A,B = ',3(e12.5,2x))
      end if
      counter=counter+1
c-------- Print Qrpa A and B matrices ------------------------------------------
      if (Iprint.ge.2) then
      write(2,30) ip1,in1,ip2,in2,Jt1,It1,Econf,
     &gph*Phme(i1,i2),Aph(i1,i2),Aph1(i1,i2),
     &gpp*ppme(i1,i2),App(i1,i2),A(i1,i2),
     &gph*Phme(i1,i2),Bph(i1,i2),Bph1(i2,i2),
     &gpp*ppme(i1,i2),Bpp(i1,i2),B(i1,i2),counter
      endif
c-------- End of cycles --------------------------------------------------------
20      END DO
10      END DO
30      format(2x,6i3,'   A :  ',7f7.3,/,2x,18x,'   B :  ',7x,6f7.3,i4)

c-------- Dimension reduction --------------------------------------------------
c       If the code runs with ALL spin configurations,
c       HERE A and B matrices can be reduced in dimension;
C       new dimension = NV1.

      modibasis=abs(ibasis)
      If (Modibasis.eq.1.and.Iprint.ge.1) write(2,60) Ispin
60    format(//,10x,' REDUCED DIMENSION OF QRPA MATRICES TO NV1
     & SPIN = ',i2,' TRANSITIONS ',//,
     &10x,'old counter   new counter    J    P    T ',/,
     &10x,'------------  ------------  ---  ---  ---',/)

        i3 = 1    
        i4 = 1  
        DO 40 i1=1,NV
           ip1=ipp(i1)
           in1=inn(i1)
           lp1=ll(ip1)
           ln1=ll(in1)
           Iprty1=(-1)**(lp1+ln1)
           Jt1=JJ(i1)
           It1=It(i1)
           IF (Jt1.NE.ISPIN)   GOTO 40     
cx         if (Iprty1.ne.Ipar) goto 40     
           NUM(i3) = i1                    
           if (modibasis.eq.1.and.Iprint.ge.1) 
     &        write(2,80) i3,i1,Jt1,Iprty1,It1
80         format(10x,2(3x,i3,6x,2x),3(i2,1x,2x))
           DO 50 i2=1,NV
              ip2=ipp(i2)
              in2=inn(i2)
              lp2=ll(ip2)
              ln2=ll(in2)
              Iprty2=(-1)**(lp2+ln2)
              Jt2=JJ(i2)
              It2=It(i2)
c ------- Selection ------------------------------------------------------------
              IF (Jt2.NE.ISPIN)   GOTO 50 
cx            if (Iprty2.ne.Ipar) goto 50 
c ------- Dimension reduced matrices -------------------------------------------
              A1(i3,i4) = A(i1,i2)
              B1(i3,i4) = B(i1,i2)
c ------------------------------------------------------------------------------
           i4=i4+1
50         continue
           i4=1
           i3=i3+1
40      continue
        NV1=i3-1                          
c number of transition of interest.
        If (Iprint.ge.1) write(2,70) ISPIN,NV1
70      format(//,10X,'No.of delta(j) = ',I1,' transitions: NV1 = '
     &         ,I3,//)


        write(50,*) 'Printing the RPA matrix'
c       do i1=1,7
         write(50,*) a(1,1),b(1,1)
c       end do
c700    format(1x,7(f8.5,1x))

        RETURN
        END

C**********************************************************************

        SUBROUTINE RPA(A,B,D,Iflag)

C	DIAGONALIZATION OF RPA MATRICES, SEE RING & SCHUCK, pag. 306.

        include 'param.qrpa'
        implicit double precision (a-h,o-z)
        dimension a(ncf,ncf),b(ncf,ncf),d(ncf)
        dimension temp(ncf,ncf)
        dimension e(ncf),c(ncf,ncf),xpy(ncf,ncf),xmy(ncf,ncf)
        dimension temp2(ncf,ncf) 
        COMMON /PAR/ GPP,GPH,DGPP,DGPH
        COMMON /PAR1/ NGPP,NGPH
        COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
        COMMON /DIM1/ PINT
        COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
        common/b_to_do/ichf,ihf,irpa
        COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
        COMMON /CTR3/ IGRAF,IBETA,IBEL
        COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del

c       The stability matrix S must be positive definite.
c       If rpa matrices A and B are real, this condition is satisfied 
c       if A+B and A-B are positive definite (see Ring Shuk, Ullah Rowe)
 
        Iflag=0
        DO 2 I=1,NV1                                              
        DO 2000 J=1,NV1
         TEMP(I,J)=A(I,J)+B(I,J)
2000    CONTINUE
2       CONTINUE

        CALL TRED2(TEMP,NV1,NCF,D,E)
        CALL TQLI(D,E,NV1,NCF,TEMP)
        IF (Iprint.ge.1) write(2,22)
22      format(/,10('*'),' A + B ROOTS ',65('*'),//,
     &       10x,'    ROOTS    ',/,10x,'-------------')
        DO 6 I=1,NV1
           if (Iprint.ge.2) write(2,23) d(i)
23      Format(10x,f9.4)
        IF(D(I).LT.0) THEN
           write(2,*) ' WARNING : NO SQUARE ROOT OF A+B !!',D(i)
           Iflag=1
           write(2,501) GPH,GPP
           RETURN
        ENDIF
6       CONTINUE

        Iflag = 0                                                    
        DO 1 I=1,NV1
        DO 1000 J=1,NV1
          C(I,J)=0.d0
          XPY(I,J)=0.d0
          XMY(I,J)=0.d0
          TEMP(I,J)=A(I,J)-B(I,J)
1000    CONTINUE
1       CONTINUE

        CALL TRED2(TEMP,NV1,NCF,D,E)
        CALL TQLI(D,E,NV1,NCF,TEMP)
        IF (Iprint.ge.1) write(2,20)
20      format(/,10('*'),' A - B ROOTS ',65('*'),//,
     &       10x,'    ROOTS    ',/,10x,'-------------')
        DO 7 I=1,NV1
           if (Iprint.ge.2) write(2,21) d(i)
21      Format(10x,f9.4)
        IF(D(I).LT.0) THEN
           write(2,*) ' WARNING : NO SQUARE ROOT OF A-B !!'
           Iflag=1 
           write(2,501) GPH,GPP
501        format(//,10x,' RPA INSTABILITY for Gph, Gpp = ',2e25.18)
         RETURN 
        ENDIF

        D(I)=dsqrt(D(I))
7       CONTINUE
        DO 8 I=1,NV1
        DO 80 J=1,NV1
        DO 800 K=1,NV1
        C(I,J)=C(I,J)+TEMP(I,K)*D(K)*TEMP(J,K)
800     CONTINUE
80      CONTINUE
8       CONTINUE
        DO 9 I=1,NV1
        DO 90 J=1,NV1
        TEMP2(I,J)=0.d0
        DO 900 K=1,NV1
        TEMP2(I,J)=TEMP2(I,J)+C(I,K)*(A(K,J)+B(K,J))
900     CONTINUE
90      CONTINUE
9       CONTINUE
        DO 92 I=1,NV1
        DO 920 J=1,NV1
        TEMP(I,J)=0.d0
        DO 9200 K=1,NV1
        TEMP(I,J)=TEMP(I,J)+TEMP2(I,K)*C(K,J)
9200    CONTINUE
920     CONTINUE
92      CONTINUE
        CALL TRED2(TEMP,NV1,NCF,D,E)
        CALL TQLI(D,E,NV1,NCF,TEMP)
        DO 10 I=1,NV1
        SGN=D(I)/DABS(D(I))
        D(I)=SGN*dsqrt(DABS(D(I)))
10      continue
        DO 11 I=1,NV1
        DO 110 J=1,NV1
        XPY(I,J)=0.d0
        DO 1100 K=1,NV1
        XPY(I,J)=XPY(I,J)+C(I,K)*TEMP(K,J)/dsqrt(DABS(D(J)))
1100    continue
110     continue
11      continue
        DO 12 I=1,NV1
        DO 120 J=1,NV1
        XMY(I,J)=0.d0
        DO 1200 K=1,NV1
        XMY(I,J)=XMY(I,J)+(A(I,K)+B(I,K))*XPY(K,J)/D(J)
1200    continue
120     continue
12      continue
        DO 13 I=1,NV1
        DO 1300 J=1,NV1
        temp(I,J)=0.5d0*(XPY(I,J)+XMY(I,J))
        temp2(I,J)=0.5d0*(XPY(I,J)-XMY(I,J))
1300    continue
13      continue
        do 130 i=1,nv1
        do 131 j=1,nv1
        a(i,j)=temp(i,j)
        b(i,j)=temp2(i,j)
131     continue
130     continue

        RETURN
        END

C**********************************************************************

      SUBROUTINE EIGSRT1(D,V,W,N,NP)

C     REARRANGE ORDER AFTER INCREASING FREQUENCY D(I)

      implicit double precision (a-h,o-z)
      dimension d(np),v(np,np),w(np,np)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /PAR/ GPP,GPH,DGPP,DGPH
      COMMON /PAR1/ NGPP,NGPH

      Call Eigsrt (D,V,W,N,NP)

c---------print the results----------------------------------------------------
      if (Iprint.ge.1) write(2,40)
40    format(/,10('*'),' REORDERED FREQUENCIES ',55('*'),//,
     &       10x,' Frequencies ',/
     &       10x,'-------------')
      if (iprint.ge.2) then
        do i=1,nv1
          write(2,51) d(i)
51      format(10x,f9.4)
        enddo
      endif
c-------------------------------------------------------------------------------
      RETURN
      END

C**********************************************************************

      SUBROUTINE ELECTRO(A,B,FREQ,BEL_is,BEL_iv,FRACTS0_is,
     &FRACTS0_iv,FRACTS1_is,FRACTS1_iv)                

C     EVALUATE TRANSITION DENSITIES AND ELECTRIC TRANSITION 
C     PROBABILITIES B(EL) 

      include 'param.qrpa'
      implicit double precision (a-h,o-z)
      dimension up(nsp),un(nsp),vp(nsp),vn(nsp)
      dimension spep(nsp),spen(nsp),ep(nsp),en(nsp)
      dimension freq(ncf),a(ncf,ncf),b(ncf,ncf)
      dimension Econf(ncf),BcsEconf(ncf)
      dimension BEL_is(ncf),BEL_iv(ncf),sq_BEL_em(ncf),  
     &    hfBEL(ncf),BcsBEL(ncf),
     &    FractS0_is(ncf),fracts0_iv(ncf),fracthfS0(ncf),
     &    fractbcsS0(ncf),
     &    FractS1_is(ncf),fracts1_iv(ncf),fracthfS1(ncf),
     &    fractbcsS1(ncf)
      dimension rmat(2,ncf,ncf)
      dimension Rhon(ncf,nnn),Rhop(ncf,nnn),rmat1(2,ncf,ncf)
      COMMON /DIM/ NINT,NV,NV1,IPROT,INEUT
      COMMON /DIM1/ PINT
      COMMON /CTR1/ IPRINT,ICTR,IWILD,IUV,ISPE,IBCS,IFILE
      common/b_to_do/ichf,ihf,irpa
      COMMON /CTR2/ IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      COMMON /CTR3/ IGRAF,IBETA,IBEL
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON /BCS2/ UP,UN,VP,VN,EP,EN,bcsp,bcsn,Fp,Fn,Pgapp,Pgapn
      COMMON /AB/ NUM(ncf)
      COMMON /PAR/ GPP,GPH,DGPP,DGPH
      COMMON /PAR1/ NGPP,NGPH
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common /hf/nmax,nocc,nunocc,norb
      common /hf1/del
      Common /dens/dt(nnn),dn(nnn),dp(nnn)
      Common /bwf/ wf(nnn,nsp),dwf(nnn,nsp)
      common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)
      common /biqi/ iqi(ncf)
      common /bcons/ icons,iovert
      common /bspur/ eta  
      common/bafa0/afa0
      common/bnosc/nosc
      common/bnuc/anucl,znucl 
      common /bgfb/ ecf(ncf),akeep(ncf)
      common/b_radii/r2l(12)
      common/brc/e_centr_min,e_centr_max
      common/b_dipo/DipoKap
      dimension sqbel_is(ncf),sqbel_iv(ncf)
      dimension belelec(ncf)
      dimension nord(ncf)
  101 format(6(1x,e12.5))

      i_check_m1=0
      open(unit=39,status='unknown',file='reduced_EM_sp_old.dat')
      write(39,*) iorb,isflip,ispin
      open(unit=40,status='unknown',file='reduced_EM_phon.dat')
      write(40,*) iorb,isflip,ispin
c---- Constants ----------------------------------------------------------------
      PI    = 3.14159265358d0       ! Pi for Pitagora
      Hbarc = 197.327053            ! hbar*c in MeV*Fermi
      protonmass  = 938.27231d0     ! proton mass in MeV/c^2
      aneutronmass = 939.56563d0    ! neutron  "   "  "
      a_med_mass = (protonmass+aneutronmass)/2.d0
      hbdm = hbarc**2/2.d0/a_med_mass
      gsp = 5.585694713d0
      gsn = -3.82608545d0
c     write(*,*) hbdm
c ---- Yj operator reduced matrix elements ------------------------------------ 
      write(2,115) 
115   format(/,10('*'),'TRANSITION DENSITY CALCULATION',35('*'),/)
       IF (Iprint.ge.2) write(2,120)
120    format(10x,' i1,iq1,L1,J1    i2,iq2,L2,J2    Geom.part   
     &   Mat.el. ',/,10x,'-----------   -------------- 
     &   ------------     -----------',/)
      do 130 k=1,nv1
         ix=ipp(num(k))
         iy=inn(num(k))
         rmat1(1,ix,iy)=0.d0
         rmat1(2,ix,iy)=0.d0
         Jt=JJ(num(k))
         lJx=lj(ix)
         lJy=lj(iy)
         llx=ll(ix)
         lly=ll(iy)
         if(isflip.eq.0)rmat1(1,ix,iy) = yl(llx,lJx,lly,lJy,Jt)         
         if(isflip.eq.1)then
c         if(ispin.eq.1.and.iorb.eq.0.and.irad.eq.0)then
c          iq1=iq(ix)
c          rmat1(1,ix,iy) = rmat_m1(llx,lJx,lly,lJy,iq1)
c         else
           rmat1(1,ix,iy) = TJL(llx,lJx,lly,lJy,Ispin,Iorb)
          end if
         rmat1(2,ix,iy) = rmat1(1,ix,iy)*iqi(num(k))
c        if (Iprint.ge.2)   
         write(2,125) 
     &      ix,iq(ix),llx,lJx,iy,iq(iy),lly,lJy,rmat1(1,ix,iy),
     &      rmat1(2,ix,iy)
130   enddo
125   format(10x,2(4i3,'/2',2x),2(f9.4,2x))

c --- Transition Density calculation -----------------------------------------
      DO 160 i1=1,NV1
        ix1=ipp(num(i1))
        iy1=inn(num(i1))
        Jt1=JJ(num(i1))
        It1=It(num(i1))
        Ex1=Spe(ix1)
        Ey1=Spe(iy1)
        freq1=freq(i1)
        rhointp=0.d0
        DO 150 i=1,nmax
c         if (i1.eq.1.and.i.le.10)then 
c         write(2,*) 'iorb,iovert=',iorb,iovert 
c         end if
          rhon(i1,i)=0.d0
          rhop(i1,i)=0.d0
          rdist=float(i)*del
          if (iorb.ne.0) then
          xop=rdist**(iorb+2+2*iovert)
            if(iovert.eq.2)then
             xop=rdist**(iorb+4)-eta*rdist**(iorb+2)
             if (i1.eq.1.and.i.le.10)then 
c             write(2,*) ' '
c             write(2,*) 'eta,rdist,xop: '
c             write(2,*) eta,rdist,xop
c             write(2,*) ' '
             end if 
            end if 
          else
          xop=rdist**(4+2*iovert)
          endif
          if(irad.eq.0)xop=dsqrt(16.d0*datan(1.d0))*rdist**2
c         if (i1.eq.1.and.i.le.10)then 
c         write(2,*) 'rdist,xop,eta=',rdist,xop,eta 
c         end if
          DO 140 i2=1,NV1  
            ix2=ipp(num(i2))
            iy2=inn(num(i2))
            Jt2=JJ(num(i2))
            It2=It(num(i2))
            Jx2=lj(ix2)
            Jy2=lj(iy2)
            lx2=ll(ix2)
            ly2=ll(iy2)
            Iprty2=(-1)**(lx2+ly2)
            Iq2=Iq(ix2)+Iq(iy2)
            Ex2=Spe(ix2)
            Ey2=Spe(iy2)
c ---- Test with H.O. wave functions -------------------------------------------
c           b=2.423d0
c           nrx2 = nn(ix2)-1
c           nry2 = nn(iy2)-1
c           wf(i,ix2)=hwf(nrx2,lx2,b,rdist)
c           wf(i,iy2)=hwf(nry2,ly2,b,rdist)
c ---- Selection rule ----------------------------------------------------------
c           if (Ex2.lt.Ey2) goto 140            ! Already in the basis
            if (Jt1.ne.Jt2) goto 140
c ---- Warnings ----------------------------------------------------------------
            if (Jt1.ne.Ispin) write(2,60) Jt1,Ispin        
60          format(/,10x,' WARNING : SPIN ERROR ! ',
     &                2x,'Jt1=',i1,' Ispin=',i1)
c ---- QRPA BEL calculation ----------------------------------------------------
            i_rpa_fac = 1
            if(isflip.eq.1)i_rpa_fac = ipot(ispin+iorb) 
            XY = a(i2,i1) + i_rpa_fac*b(i2,i1)               
            Iexp = (-lx2+ly2+iorb)/2 
            Isgn = Ipot(Iexp)                                
            itau = ipot(isflip)
            if(ibasis.lt.0)then 
             if (iq2.eq.2) then                              
              UVp = Up(ix2)*Vp(iy2) + itau*Vp(ix2)*Up(iy2)   
              RHOp(i1,i)= RHOp(i1,i) + Isgn*rmat1(1,ix2,iy2)*XY*UVp*  
     &         WF(i,ix2)*WF(i,iy2)/(rdist**2)
             else if (Iq2.eq.0) then                        
              UVn = Un(ix2)*Vn(iy2) + itau*Vn(ix2)*Un(iy2)   
              RHOn(i1,i)= RHOn(i1,i) + Isgn*rmat1(1,ix2,iy2)*XY*UVn*
     &         WF(i,ix2)*WF(i,iy2)/(rdist**2)
             endif
            else
c            if(i.eq.1)then 
c             write(2,*) ix2,iy2
c             write(2,*) Un(ix2),Vn(ix2),Up(ix2),Vp(ix2)
c             write(2,*) Un(iy2),Vn(iy2),Up(iy2),Vp(iy2)
c            end if
             UV = Up(ix2)*Vn(iy2) + itau*Vp(ix2)*Un(iy2)
             if(i.eq.1) write(2,*) UV
             RHOp(i1,i)= RHOp(i1,i) + Isgn*rmat1(1,ix2,iy2)*XY*UV*
     &         WF(i,ix2)*WF(i,iy2)/(rdist**2)
             RHOn(i1,i)= 0.d0
            end if

140        END DO
           rhop(i1,i)=-rhop(i1,i)/sqrt(float(2*Ispin+1))
           rhon(i1,i)=-rhon(i1,i)/sqrt(float(2*Ispin+1))
150   ENDDO
160   ENDDO
c     if(ngph.eq.1.and.ngpp.eq.1) then
c     write(87,*)freq(1)
c     do i=1,nmax
c     write(87,*)i,rhop(1,i),rhon(1,i)
c     enddo
c     endif
    
c --- BEL calculation ----------------------------------------------------------
      HFS0    = 0.d0
      BCSS0   = 0.d0
      S0_is   = 0.d0
      S0_iv   = 0.d0
      HFS1    = 0.d0
      BCSS1   = 0.d0
      S1_is   = 0.d0
      S1_iv   = 0.d0
      Sm1_is  = 0.d0
      Sm1_iv  = 0.d0
      S3_is  = 0.d0
      S3_iv  = 0.d0

      s0_is_part = 0.d0
      s1_is_part = 0.d0
      s0_iv_part = 0.d0
      s1_iv_part = 0.d0

      do 200 i1=1,nv1
        ix1=ipp(num(i1))
        iy1=inn(num(i1))
        Ex1=Spe(ix1)
        Ey1=Spe(iy1)
        rhoelec=0.0d0
        rhoint_is=0.d0                      
        rhoint_iv=0.d0 
        Rdint=0.d0
        fac1=2.d0*znucl/anucl
        fac2=2.d0*(anucl-znucl)/anucl

         do i=1,nmax 
          rdist=float(i)*del
          if (iorb.ne.0) then
          xop=rdist**(iorb+2+2*iovert)
            if(iovert.eq.2)then
             xop=rdist**(iorb+4)-eta*rdist**(iorb+2)
            end if 
          else
            xop=rdist**(4+2*iovert)
          endif
          if(irad.eq.0)xop=dsqrt(16.d0*datan(1.d0))*rdist**2
          rhoelec=rhoelec+(rhop(i1,i))*xop  
          rhoint_is=rhoint_is + (rhop(i1,i)+rhon(i1,i))*xop  
          if(iorb.eq.1.and.isflip.eq.0)then     
c IV DIPOLE 
           rhoint_iv=rhoint_iv + (fac1*rhon(i1,i)-fac2*rhop(i1,i))*xop 
          else
           rhoint_iv=rhoint_iv + (rhop(i1,i)-rhon(i1,i))*xop
          end if
          Rdint=Rdint+wf(i,ix1)*wf(i,iy1)*xop/(rdist**2)
        enddo
        
        rhoelec=rhoelec*del
        rhoint_is=rhoint_is*del
        rhoint_iv=rhoint_iv*del
        write(99,*) i1,rhoelec,rhoint_is,rhoint_iv
        rdint=rdint*del
        rmat(1,ix1,iy1)=rmat1(1,ix1,iy1)*rdint       
        rmat(2,ix1,iy1)=rmat1(2,ix1,iy1)*rdint       
        qrhoelec=rhoelec**2
        qrhoint_is=rhoint_is**2
        qrhoint_iv=rhoint_iv**2

        sqbel_is(i1)=rhoint_is*sqrt(float(2*ispin+1))
        sqbel_iv(i1)=rhoint_iv*sqrt(float(2*ispin+1))
        bel_is(i1)=qrhoint_is*(2*ispin+1)
        bel_iv(i1)=qrhoint_iv*(2*ispin+1)

        belelec(i1)=qrhoelec*(2*ispin+1)
c       Remember that the transition density is normalized 
c       in order to reproduce the electric strength B(EL,J->0),  
c       while B(EL,0->J)=(2J+1)*B(EL,J->0)

c ----- Single Particle limits : HF and BCS ------------------------------------
        HfBEL(I1) = 0.d0                                   ! HF
        If (ix1.eq.iy1) goto 11
        Econf(i1)=dabs(Ex1-Ey1)                            ! Conf.energy
        hfocc = dsqrt((1.d0-focc(ix1))*(focc(iy1)))+       
     &          dsqrt((focc(ix1))*(1.d0-focc(iy1)))
        hfBEL(i1) = (rmat(1,ix1,iy1)*hfocc)**2             ! Hf B(E,L)
        Iq1=Iq(ix1)+Iq(iy1)
        itau = ipot(isflip)
        if (Iq1.eq.2) then                                 ! Protons
           bcsEconf(I1) = Ep(Ix1) + Ep(Iy1)
           UVp = Up(Ix1)*Vp(Iy1) + itau*Vp(Ix1)*Up(Iy1)    
           BCSBEL(I1) = (rmat(1,Ix1,Iy1)*UVp)**2           ! Bcs BEL
        else if (Iq1.eq.0) then                            ! neutrons
           bcsEconf(I1) = En(Ix1) + En(Iy1)
           UVn = Un(Ix1)*Vn(Iy1) + itau*Vn(Ix1)*Un(Iy1)    ! Bcs factor
           BCSBEL(I1) = (rmat(1,Ix1,Iy1)*UVn)**2           ! Bcs BEL
        endif
11      continue

c ---- Strength Function Moments -------------------------------------------
        if(ispin.ne.1.or.i1.ne.1)  then
        hfS0  = hfS0  + hfBEL(i1)                    ! Hf limit
        BcsS0 = BcsS0 + bcsBEL(i1)                   ! bcs limit
        hfS1  = hfS1  + Econf(i1)*hfBEL(i1)          ! Hf limit
        BcsS1 = BcsS1 + BcsEconf(i1)*bcsBEL(i1)      ! bcs limit
        S0_is = S0_is + BEL_is(I1)                   ! NEWSR 
        S1_is = S1_is + Freq(i1)*BEL_is(i1)          ! EWSR  
        Sm1_is = Sm1_is + BEL_is(i1)/Freq(i1)        ! m_{-1} 
        S3_is = S3_is + BEL_is(i1)*Freq(i1)**3       ! m_3 
        S0_iv = S0_iv + Bel_iv(i1)                   ! NEWSR 
        S1_iv = S1_iv + Freq(i1)*BEL_iv(i1)          ! EWSR 
        Sm1_iv = Sm1_iv + BEL_iv(i1)/Freq(i1)        ! m_{-1}
        S3_iv = S3_iv + BEL_iv(i1)*Freq(i1)**3       ! m_3 
        if(freq(i1).gt.e_centr_min.and.
     &  freq(i1).lt.e_centr_max)then
         s0_is_part = s0_is_part + BEL_is(I1)
         s1_is_part = s1_is_part + BEL_is(I1)*Freq(i1)
         s0_iv_part = s0_iv_part + BEL_iv(I1)
         s1_iv_part = s1_iv_part + BEL_iv(I1)*Freq(i1)
        endif
        endif
200   enddo

c --- Fraction of the sum rules ------------------------------------------------
      do k=1,nv1
        if (S0_is.ne.0.d0) then
          FractS0_is(k) = BEL_is(k)/S0_is *100.d0
        else 
          FractS0_is(k) = 0.d0
        endif
        if (S1_is.ne.0.d0) then
          FractS1_is(k) = freq(k)*BEL_is(k)/S1_is *100.d0
        else 
          FractS1_is(k)=0.d0
        endif
        if (S0_iv.ne.0.d0) then
          FractS0_iv(k) = BEL_iv(k)/S0_iv *100.d0
        else 
          FractS0_iv(k) = 0.d0
        endif
        if (S1_iv.ne.0.d0) then
          FractS1_iv(k) = freq(k)*BEL_iv(k)/S1_iv *100.d0
        else 
          FractS1_iv(k)=0.d0
        endif
      enddo

      do k=1,nv1
        if (HfS0.ne.0.d0) then
          FracthfS0(k)  =  hfBEL(k) /HfS0 *100.d0  ! Hf limit
        else 
          FracthfS0(k) = 0.d0
        endif
        if (HfS1.ne.0.d0) then
          FracthfS1(k)  = Econf(k)*hfBEL(k) /HfS1 *100.d0  !  "
        else 
          FracthfS1(k)=0.d0
        endif
      enddo

      do k=1,nv1
        if (bcsS0.ne.0.d0) then
          FractbcsS0(k) = bcsBEL(k)/bcsS0*100.d0  ! BCS limit
        else 
          FractbcsS0(k) = 0.d0
        endif
        if (bcsS1.ne.0.d0) then
          FractbcsS1(k) = BcsEconf(k)*bcsBEL(k)/bcsS1*100.d0  !  "
        else 
          FractbcsS1(k)=0.d0
        endif
      enddo
  
c --- Theoretical EWSR (see B&M2,p.399->; Rowe,p.29->, R&S,p.332_>) ------------
c --- (Caution : funziona per multipolarita` Jt2 fissate a priori ! )-----------
      Xintp  = 0.d0
      Xintn  = 0.d0
      Xint   = 0.d0
      Xnormp = 0.d0
      Xnormn = 0.d0
      Xnorm  = 0.d0
      iiexp = (2*Jt2-2)
      if (Jt2.eq.0) iiexp = 2                       
c      exp = float(2*Jt2-2)
c      if (Jt2.eq.0) exp = 2.D0                     
      do 40 k=1,nmax                               
         rdist  = k*del 
         rdist2 = rdist**2 
c         op     = rdist**(exp)
         op     = rdist**(iiexp)
         Xintp  = Xintp  + rdist2 * dp(k) * op
         Xintn  = Xintn  + rdist2 * dn(k) * op
         Xnormp = Xnormp + rdist2 * dp(k)
         Xnormn = Xnormn + rdist2 * dn(k)
40    enddo
      Xmeanp = Xintp/Xnormp
      Xmeann = Xintn/Xnormn
      fact   = Jt2*(2*Jt2+1)**2/(4.d0*pi)*hbarc**2/2.d0
      Xewsrp = fact/protonmass *Iprot*Xmeanp
      Xewsrn = fact/aneutronmass*Ineut*Xmeann
      Xewsr  = Xewsrp + Xewsrn
c     write(2,*) '>>> ONLY THE IS THEORETICAL SUM RULE IS CALCULATED'
c --- Percentual of Theoretical S1 ---------------------------------------------
      if (xewsr.ne.0d0) then
      diffpercS1_is = dabs(S1_is /Xewsr)*100.d0 
      hfdiffpercS1  = dabs(hfS1 /Xewsr)*100.d0 
      BcsdiffpercS1 = dabs(BcsS1/Xewsr)*100.d0 
    
c     write(2,*)'rel. diff. with theor. is.ewsr',(S1_is-xewsr)/xewsr
c     write(2,*)fact,protonmass,xmeanp,iprot,ineut
c     write(2,*)'fact,prmass,rmedio2,zprot,zneu'
      endif
c --- Try to get the B in a different way 
      do i1=1,nv1

        prova_is = 0.d0
        prova_iv = 0.d0 
        prova_em = 0.d0
c       prova_em is the electromagnetic reduced matrix element associated
c       with the L, S and J quantum numbers that are in the (Q)RPA input
        do i2=1,NV1
            ix2=ipp(num(i2))
            iy2=inn(num(i2))
            Jt2=JJ(num(i2))
            Jx2=lj(ix2)
            Jy2=lj(iy2)
            lx2=ll(ix2)
            ly2=ll(iy2)
            Iprty2=(-1)**(lx2+ly2)
            Iq2=Iq(ix2)+Iq(iy2)
            Ex2=Spe(ix2)
            Ey2=Spe(iy2)
            i_rpa_fac = 1 
            if(isflip.eq.1)i_rpa_fac = ipot(ispin+iorb) 
            XY = a(i2,i1) + i_rpa_fac*b(i2,i1)               
            Iexp = (-lx2+ly2+Iorb)/2
            Isgn = Ipot(Iexp)                               
            itau = ipot(isflip)
            if (iq2.eq.2) then                              
              cccc = 0.d0
              UVp = Up(ix2)*Vp(iy2) + itau*Vp(ix2)*Up(iy2)  
              prova_is = prova_is + Isgn*rmat(1,ix2,iy2)*XY*UVp
              prova_iv = prova_iv + Isgn*rmat(2,ix2,iy2)*XY*UVp
              if(isflip.eq.0)then
               idummy=0
c the quantity in the following line is < ix2 || i^L O(EL) || iy2 > already in phase convention II
               cccc = rmat(1,ix2,iy2)*Isgn
               write(39,9701) 
     &         Iq(ix2),nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2)
     &         ,lj(iy2),
     &         idummy,Ispin,rmat(1,ix2,iy2)*Isgn
              end if
 9701         format(10(1x,i3),2x,e13.6) 
              if(isflip.eq.1) then
               idummy=1
c the quantity is now < ix2 || i^L-1 O(ML) || iy2 > already in phase convention II
               iphase1 = (lj(iy2)-lj(ix2))/2 + Ispin  - 1
               iphase1 = 1-2*mod(iphase1,2)
               iphase2 = (ll(iy2) - ll(ix2) + Ispin- 1) / 2
               iphase2 = 1-2*mod(iphase2,2)
               fact1 = dsqrt(dfloat(lj(iy2)+1)*dfloat(2*Ispin+1))
     &         /dsqrt(16.d0*datan(1.d0))
               Rdint=0.d0
               do ipoints=1,nmax
                rdist=dfloat(ipoints)*del
                Rdint=Rdint+wf(ipoints,ix2)*wf(ipoints,iy2)
     &          *rdist**(Ispin-1)
               enddo
               curl = curl1 + curl2 
               curl1 = (gsp - 2.d0/dfloat(Ispin+1))*dfloat(Ispin)/2.d0*
     &         cofcg(dfloat(lj(iy2))/2.d0,dfloat(Ispin),dfloat(lj(ix2))
     &         /2.d0,
     &         0.5d0,0.d0,0.5d0)*squarebra
               iphase3 = 1
               if(lj(iy2).lt.2*ll(iy2))iphase3 = -1 
               iphase4 = 1-2*mod((lj(ix2)+lj(iy2))/2-Iorb,2)
               squarebra = 1.d0 + 1.d0/2.d0/dfloat(Ispin)*iphase3
     &         *((lj(iy2)+1)+iphase4*(lj(ix2)+1))
               curl2 = 0.d0
               cccc = iphase1*iphase2*fact1*Rdint*curl
               write(39,9701) Iq(ix2),
     &         nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2),lj(iy2),
     &         idummy,Ispin,cccc
              end if
              if(I1.eq.1)then
              write(2,9700) i2,ix2,iy2,a(i2,i1),UVp,rmat(1,ix2,iy2),Isgn
 9700         format(3i3,3(1x,e12.5),1x,i2)
              end if
              prova_em = prova_em + cccc*XY*UVp
            else if (Iq2.eq.0) then                          ! neutrons
              cccc = 0.d0
              UVn = Un(ix2)*Vn(iy2) + itau*Vn(ix2)*Un(iy2)      
              prova_is = prova_is + Isgn*rmat(1,ix2,iy2)*XY*UVn
              prova_iv = prova_iv + Isgn*rmat(2,ix2,iy2)*XY*UVn
              if(isflip.eq.1) then
               idummy=1
c the quantity is now < ix2 || i^L-1 O(ML) || iy2 > already in phase convention II
               iphase1 = (lj(iy2)-lj(ix2))/2 + Ispin  -1
               iphase1 = 1-2*mod(abs(iphase1),2)
               iphase2 = (ll(iy2) - ll(ix2) + Ispin - 1) / 2
               iphase2 = 1-2*mod(abs(iphase2),2)
               fact1 = dsqrt(dfloat(lj(iy2)+1)*dfloat(2*Ispin+1))
     &         /dsqrt(16.d0*datan(1.d0))
               Rdint=0.d0
               do ipoints=1,nmax
                rdist=dfloat(ipoints)*del
                Rdint=Rdint+wf(ipoints,ix2)*wf(ipoints,iy2)
     &          *rdist**(Ispin-1)
               enddo
               curl = curl1 + curl2
               curl1 = gsn*dfloat(Ispin)/2.d0*
     &         cofcg(dfloat(lj(iy2))/2.d0,dfloat(Ispin),dfloat(lj(ix2))
     &         /2.d0,
     &         0.5d0,0.d0,0.5d0)*squarebra
               iphase3 = 1
               if(lj(iy2).lt.2*ll(iy2))iphase3 = -1
               iphase4 = 1-2*mod((lj(ix2)+lj(iy2))/2-Iorb,2)
               squarebra = 1.d0 + 1.d0/2.d0/dfloat(Ispin)*iphase3
     &         *((lj(iy2)+1)+iphase4*(lj(ix2)+1))
               curl2 = 0.d0
               cccc = iphase1*iphase2*fact1*Rdint*curl
               write(39,9701) Iq(ix2),
     &         nn(ix2),ll(ix2),lj(ix2),Iq(iy2),nn(iy2),ll(iy2),lj(iy2),
     &         idummy,Ispin,cccc
!               write(*,*) iphase3, ' ',iphase4
              end if
              prova_em = prova_em + cccc*XY*UVn
              if(I1.eq.1)then
              write(2,9700) i2,ix2,iy2,a(i2,i1),UVn,rmat(1,ix2,iy2),Isgn
              end if
            endif

        end do
c_TEST
        write(2,*) 'b_is and prova_is**2',bel_is(i1),prova_is**2
        write(2,*) 'b_iv and prova_iv**2',bel_iv(i1),prova_iv**2

        sq_BEL_em(i1) = prova_em 
        write(40,9702) iorb,isflip,Ispin,i1,prova_em,prova_em**2
 9702   format(3(1x,i3),2x,i5,2(1x,e13.6))

        end do

c --- Write the results --------------------------------------------------------
      write(2,70)
70    format(/,10('*'),' QRPA ELECTROMAGNETIC TRANSITIONS ',31('*'),//,
     &'  N   J       E          BEL_is       %S0_is      BEL_em       
     &BEL_iv      %S0_iv  ',/,
     &' --- ---   -------    ---------    ---------    --------- 
     &   -------       --------',/)
        do k=1,nv1
        write(2,80) k,JJ(num(k)),Freq(k),BEL_is(k),
     &  fractS0_is(k),belelec(k),bel_iv(k),fractS0_iv(k)
        enddo
80    Format(i4,i4,6(e12.5,1x))

      write(2,102) S0_is,s0_iv,s1_is,s1_iv,Sm1_is,S3_is,Gpp,Gph
      write(2,*) '  m(-1)(IS) =',sm1_is
      write(2,*) '  m(1) (IS) =',S1_is        
      write(2,*) '  m(-1)(IV) =',sm1_iv
      write(2,*) '  m(1) (IV) =',S1_iv        
      icho=iorb
c     write(*,*) hbdm
c     write(*,*) anucl,iorb
c     write(*,*) icho,r2l(icho)
      write(2,*) 'quantum numbers: ',iorb,ispin,isflip,ipar,irad

c     Monopole
      if(iorb.eq.0.and.ispin.eq.0.and.isflip.eq.0.and.ipar.eq.1.
     &and.irad.eq.1) then
       am1_dc = hbdm*anucl*r2l(2)/pi
       go to 81
      end if 

c     IV Dipole
      if(iorb.eq.1.and.ispin.eq.1.and.isflip.eq.0.and.ipar.eq.-1.
     &and.irad.eq.1.and.iovert.eq.0) then
       write(2,*) 'm(1) D.C. for the IV dipole'
       akappa = DipoKap
c      write(2,*) 60.d0*(anucl-znucl)*znucl/anucl
       write(2,*) 'k = ',akappa
c      write(2,*) 60.d0*(anucl-znucl)*znucl/anucl*(1+akappa)
       am1_dc = 60.d0*(anucl-znucl)*znucl/anucl*(1+akappa)
c      write(2,*) am1_dc
       go to 81
      end if 

c     IS Dipole
      if(iorb.eq.1.and.ispin.eq.1.and.isflip.eq.0.and.ipar.eq.-1.
     &and.irad.eq.1.and.iovert.eq.2) then
       icho=iorb+iovert
       write(2,*) 'm(1) D.C. for the IS dipole'
       am1_dc = hbdm*anucl*(33.d0*r2l(icho)-25.d0*r2l(2)**2)
     & /4.d0/pi
       write(2,*) am1_dc
       go to 81
      end if 

c     L equal to or greater than 2
      if(iorb.eq.ispin.and.isflip.eq.0.and.iorb.ne.0.and.irad.eq.1)then
       am1_dc = hbdm*anucl*iorb*(2*iorb+1)**2
     & /4.d0/pi*r2l(icho)
       go to 81
      end if
   
      go to 82 

81    write(2,*) '  m(1) D.C. =',am1_dc
      if(iorb.ne.1.or.ispin.ne.1.or.isflip.ne.0.or.ipar.ne.-1.
     1or.irad.ne.1)then
       write(2,*) '  % exhau.  =',S1_is/am1_dc*100.d0,' %'
      else
       write(2,*) '  % exhau.  =',S1_iv/am1_dc*100.d0,' %'
      end if
82    write(2,*) '  m(3)      =',S3_is
      write(2,*) '  sqrt[m(1)/m(-1)] (IS) =',dsqrt(s1_is/sm1_is)
      write(2,*) '  m(1)/m(0) (IS) =',s1_is/s0_is
      write(2,*) '  sqrt[m(3)/m(1)] (IS)  =',dsqrt(s3_is/s1_is)
      write(2,*) '  sqrt[m(1)/m(-1)] (IV) =',dsqrt(s1_iv/sm1_iv)
      write(2,*) '  m(1)/m(0) (IV) =',s1_iv/s0_iv
      write(2,*) '  sqrt[m(3)/m(1)] (IV)  =',dsqrt(s3_iv/s1_iv)
102   FORMAT(/,10x,'QRPA STRENGTH FUNCTION MOMENTS :',//,
     &10x,'     S0_is     S0_iv      S1_is    S1_iv      S_{-1}_is
     &     S3_is',/,
     &10x,'----------- ---------- ---------- ---------- ----------
     & --------',/,
     &3x,'QRPA :',6(e12.5,1x),//,
     &3x,'Gpp, Gph = ',2e25.18)

      if(e_centr_min.gt.0.1d0.and.e_centr_max.gt.0.1d0)then
      write(2,*)
      write(2,*) 'Between ',e_centr_min,' and ',e_centr_max,' :'
      write(2,*) 'm(0),m(1),m(1)/m(0) IS = '
      write(2,*) s0_is_part,s1_is_part,s1_is_part/s0_is_part
      write(2,*) 'm(0),m(1),m(1)/m(0) IV = '
      write(2,*) s0_iv_part,s1_iv_part,s1_iv_part/s0_iv_part
      write(2,*)
      end if

c --- Print out all t.d.'s -----------------------------------------------------
c --- (and more data about phonons) --------------------------------------------
c --- AND CORRELATION ENERGY (Sep. 2001) ---------------------------------------
      icount=0
c     write(2,*) ephocut,S0cut
      open(unit=95,status='unknown',file='TD.dat')
      open(unit=96,status='unknown',file='corr.dat')

      do i1=1,nv1
      nord(i)=i1
      end do

      do 226 i1=1,nv1
      k=i1
      x=ecf(i1)
      aaa=akeep(i1)
       do 227 i2=i1,nv1
       if(x-ecf(i2)) 227,227,228     
  228  k=i2
       x=ecf(i2)
       aaa=akeep(i2)    
  227  continue
      kk=nord(k)
      nord(k)=nord(i1)   
      nord(i1)=kk
      ecf(k)=ecf(i1)
      akeep(k)=akeep(i1)
      ecf(i1)=x
      akeep(i1)=aaa
  226 continue

      ecorr=0.d0
      ecorr2=0.d0
      do 103 i1=1,nv1
      write(2,*) ecf(i1),akeep(i1),freq(i1)
      ecorr=ecorr+0.5d0*(freq(i1)-ecf(i1)-akeep(i1))
      sumy2=0.d0
      do 104 i2=1,nv1
       sumy2=sumy2+b(i2,i1)**2
  104 continue
      ecorr2=ecorr2-freq(i1)*sumy2
      write(96,*) freq(i1),float(2*ispin+1)*ecorr,
     1float(2*ispin+1)*ecorr2
      if(freq(i1).gt.ephocut)go to 103
      if(fractS0_is(i1).lt.S0cut.and.fractS0_iv(i1).lt.S0cut)go to 103
      icount=icount+1
c     write(95,106) icount,iorb,ispin,isflip,freq(i1)
C_TEST
      write(95,107) freq(i1),fractS0_is(i1),fractS1_is(i1)
c     write(95,*) fractS0_is(i1),fractS0_iv(i1)
      do 105 i=1,nmax
       rdist = float(i)*del
       write(95,*) rdist,rhop(i1,i),rhon(i1,i)
  105 continue
  103 continue
  106 format(1x,4(i4,1x),e12.5)    
  107 format(2x,e12.5,4x,e12.5,4x,e12.5)

      close(96)
      close(95) 

c--write in file 100 phonon quantum number, energies and transition densities---
c     WARNING! IF IFILE=1 ANGULAR MOMENTUM MUST BE FIXED --> JT=ISPIN
c-------------------------------------------------------------------------------

        if (nv1.eq.0) go to 112
        if (ifile.eq.0) go to 112
        write(2,*) '>>> DO NOT USE IFILE=1 YET ! '
        stop

c  in phonon file 100 are written different quantities, according with 
c  itot value. In order to allow reading phonon file in the correct way, it is 
c  necessary to know if itot=0 or =1
 
c       do i1=1,nv1
c         if (freq(i1).lt.ephocut.and.fractS1(i1).gt.S1cut) then ! energy and 
c           write (100,*) ispin,ipar,freq(i1),fractS1(i1)  ! fractS1 cutoff
c           write (100,*) (rho(i1,k),k=1,nmax)
c	  endif
c	enddo
c	CLOSE (100)
112     continue

CC--- Write the results in formatted output named rpa.dat ----------------------
      nc=nv
      if(irpa.eq.2)nc=2*nv
      open(unit=11,status='unknown',file='rpa.dat',form='formatted')
      if((ngph.eq.1.and.ngpp.eq.1).or.freq(1).le.0.00001) then
      do i=1,nv1
       write(11,*) iorb,ispin,isflip,freq(i)
       write(11,725) 
     & bel_is(i),sq_BEL_em(i),bel_is(i)/s0_is,bel_iv(i)/s0_iv
       write(11,*) nv1
725    format(2x,4(e12.5,3x))
       do j=1,nv1
       ip1=ipp(j)
       iq1=iq(ip1)
       nn1=nn(ip1)
       ll1=ll(ip1)
       jj1=lj(ip1)
       ncode1 = (-1)**iq1*(nn1*10000+ll1*100+jj1)
       in1=inn(j)
       iq2=iq(in1)
       nn2=nn(in1)
       ll2=ll(in1)
       jj2=lj(in1)
       ncode2 = (-1)**iq2*(nn2*10000+ll2*100+jj2)
       write(11,726)j,ncode1,ncode2,a(j,i),b(j,i)
726    format(3(2x,i12),3x,2(e12.5,2x))
       enddo
      end do
      endif
      close(39)
      close(40)
      close(11)

      RETURN
      END

C**********************************************************************

      FUNCTION GAMMA(y)

C     if x>1 --> Stirling formula for ln(x!)
C     if x<1 --> see Hastings : Approx for Digital Computers 
C                (Princeton Univ. Press, 1955), or Arfkin p.465

      pi=2.d0*dsin(1.d0)
      gamma=1.0
      X=Y
      if (x.le.0) stop ' ERROR in GAMMA : X<0 !! '
      if (x.gt.1.0) then 
         gamma=sqrt(2.*pi)*x**(x+0.5)*exp(-x)*
     &         exp(1./(12.*x)-1./(360.*x**3)-1./(1260.*x**5))
      else
         gamma=1-0.577191652*x+0.988205891*x**2-0.897056937*x**3
     &         +0.918206857*x**4-0.756704078*x**5+0.482199394*x**6
     &         -0.193527818*x**7+0.035868343*x**8
      endif
      gamma=gamma/y

      return
      end

C**********************************************************************

      FUNCTION GAMMACOMP(x,y)

C     Absolute value of complex Euler Gamma Function |Gamma(x+iy)|
C     (exercise 10.1.19 pag.458 Arfkin)   (from E.Ormand)

      fact=1.
      I=-1
2     I=I+1
      a=float(i)
      OLD=FACT
      fact=fact/sqrt(1.+(y/(x+a))**2)
      DIFF=ABS(FACT-OLD)
      IF (DIFF.GT.1.E-6) goto 2
      Gammacomp=fact*gamma(1.+x)/x

      return
      end

C**********************************************************************

      SUBROUTINE subphme(i_t3,Ia,Ib,Ic,Id,Jt,nr,dr,res)

      implicit double precision (a-h,o-z)
c     real*4 sla,slb,slc,sld,sja,sjb,sjc,sjd,sjt,sj,c1,c2,c3,c4,
c    &       sl,slk1,slk2,sll,slka,slkb,slkc,slkd
      include 'param.qrpa'
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1calfv(7),cbetv(7) 
      COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common /bcons/ icons,iovert
      dimension hh(9),vv(7,3)

      res=0.d0               ! the result is prepared
      do i=1,9
         hh(i)=0.d0
      enddo
 
c---- Isospin part:
c---- the quantities iq* are tau_z(*)      
      iqa=1-2*iq(ia)
      iqb=1-2*iq(ib)
      iqc=1-2*iq(ic)
      iqd=1-2*iq(id)
      
      if((iqa.eq.iqd).and.(iqb.eq.iqc)) then
       elias=1.d0
       eliav=float(iqa*iqb)
      else
       elias=0.d0
       eliav=0.d0
      endif                                                                
      
      if(((iqa+iqd).eq.0).and.((iqb+iqc).eq.0)) then
       eliav=eliav+1.d0+float(iqa*iqc)
      endif

      ieliav0=(iqa-iqb)**2+(iqb-iqc)**2+(iqc-iqd)**2
      if(ieliav0.eq.0) then
       ieliav1=2*iqa
      else
       ieliav1=0
      endif

      do 1039 j=1,nnn
       f0(j,3)=f0(j,1)*elias+eliav*f0(j,2)+float(ieliav1)*corr_res(j)
       g0(j,3)=g0(j,1)*elias+eliav*g0(j,2)
 1039 continue
C ------ Matrix Element calculation --------------------------------------------
      call rdint1(i_t3,ia,ib,ic,id,nr,dr)        
      djp1=2.d0*jt+1.d0
      la=ll(ia)
      lb=ll(ib)
      lc=ll(ic)
      ld=ll(id)
      ja=lj(ia)
      jb=lj(ib)
      jc=lj(ic)
      jd=lj(id)
      sla=float(la)
      slb=float(lb)
      slc=float(lc)
      sld=float(ld)
      sja=float(ja)/2.d0
      sjb=float(jb)/2.d0
      sjc=float(jc)/2.d0
      sjd=float(jd)/2.d0
      sjt=float(jt)
      hh(8)=ri0*yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/djp1 
c   1 continue
      ladm=iabs(la-ld)
      ladp=la+ld
      lbcm=iabs(lb-lc)
      lbcp=lb+lc
      j1m=iabs(jt-1)
      j1p=jt+1
      l1=max0(ladm,lbcm,j1m)
      l2=min0(ladp,lbcp,j1p)
      x9=0.d0
      if(l1.gt.l2) go to 2
      do 3 l=l1,l2
      x9=x9+tjl(la,ja,ld,jd,jt,l)*tjl(lc,jc,lb,jb,jt,l)
3     continue
2     hh(9)=ris*x9/djp1                         
      l1=max0(ladm,lbcm)
      l2=min0(ladp,lbcp)
      if(l1.gt.l2) go to 32
      do i1=1,7
         do i2=1,3
            vv(i1,i2)=0.d0
         enddo
      enddo
      cofs=(ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0)
      cofs=dsqrt(cofs)*6.d0
      is=la+lc+(jb+jd)/2+1
      si=1.d0-2.d0*mod(is,2)
      call sixj(sja,sla,0.5d0,sld,sjd,sjt,c1)
      call sixj(sjc,slc,0.5d0,slb,sjb,sjt,c2)
      c1=c1*c2
      cof0=si*cofs*c1/6.d0
      ilpt=1
      ilgd=3
      do il=ilpt,ilgd              
         l=jt-2+il
         if(l.lt.0) go to 31
         if((l-l1)*(l-l2).gt.0) go to 31
         vv(1,il)=rx1*redy(la,l,ld)*redy(lc,l,lb)/(2.d0*l+1.d0) 
31    enddo
      do 8 iup=2,5
      go to (9,9,10,11,30), iup    
    9 ka=ia
      kb=ib
      kc=ic
      kd=id
      go to 13
   10 ka=id
      kb=ic
      kc=ib
      kd=ia
      go to 13
   11 ka=ia
      kb=ic
      kc=ib
      kd=id
      go to 13
   30 ka=ib
      kb=id
      kc=ia
      kd=ic
   13 lka=ll(ka)
      lkb=ll(kb)
      lkc=ll(kc)
      lkd=ll(kd)
      slka=float(lka)
      slkb=float(lkb)
      slkc=float(lkc)
      slkd=float(lkd)
      call rdint2(ka,kb,kc,kd,nr,dr)
      xlkab=(2.d0*lka+1.d0)*(2.d0*lkb+1.d0)
      xlkab=dsqrt(xlkab)

      do 14 il=ilpt,ilgd                            
      l=jt-2+il
      sl=float(l)
      if(l.lt.0) go to 14
      if((l-l1)*(l-l2).gt.0) go to 14
      sgn1=1.d0-2.d0*mod(l,2)
      if(iup.eq.4.or.iup.eq.5) sgn1=1.d0

      do 15 ik1=1,2
      k1=lka-3+2*ik1
      if(k1.lt.0) go to 15
      slk1=float(k1)
      call troisj(slk1,slka,1.d0,0.d0,0.d0,0.d0,c1)
      cofk1=(1.d0-2.d0*mod(k1,2))*dsqrt(2.d0*k1+1.d0)*c1

      do 16 ik2=1,2
      k2=lkb-3+2*ik2
      if(k2.lt.0) go to 16
      slk2=float(k2)
      call troisj(slk2,slkb,1.d0,0.d0,0.d0,0.d0,c2)
      cofk2=(1.d0-2.d0*mod(k2,2))*dsqrt(2.d0*k2+1.d0)*c2
      ll1m=iabs(l-1)
      ll2m=iabs(lkd-k1)
      ll3m=iabs(lkc-k2)
      ll1p=l+1
      ll2p=lkd+k1
      ll3p=lkc+k2
      llm=max0(ll1m,ll2m,ll3m)
      llp=min0(ll1p,ll2p,ll3p)
      if(llm.gt.llp) go to 16

      do 17 l1l=llm,llp
      sll=float(l1l)
      call sixj(sll,1.d0,sl,slka,slkd,slk1,c1)
      call sixj(sll,1.d0,sl,slkb,slkc,slk2,c2)
      c3=c1*c2
      vv(iup,il)=vv(iup,il)+sgn1*xlkab*cofk1*cofk2*c3*
     &           redy(lkd,l1l,k1)*redy(lkc,l1l,k2)*rx2(ik1,ik2)
   17 continue                           ! I rx2(*,*)
   16 continue
   15 continue
   14 continue
    8 continue

      do 18 iup=6,7
      ka=ia
      kb=ib
      kc=ic
      kd=id
      if(iup.eq.6) go to 19
      ka=ib
      kb=ia
      kc=id
      kd=ic
   19 lka=ll(ka)
      lkb=ll(kb)
      lkc=ll(kc)
      lkd=ll(kd)
      slka=float(lka)
      slkb=float(lkb)
      slkc=float(lkc)
      slkd=float(lkd)
      call rdint2(ka,kd,kb,kc,nr,dr)
      xlkad=(2.d0*lka+1.d0)*(2.d0*lkd+1.d0)
      xlkad=dsqrt(xlkad)

      do 20 il=ilpt,ilgd
      l=jt-2+il
      sl=float(l)
      if(l.lt.0) go to 20
      if((l-l1)*(l-l2).gt.0) go to 20
      cofl=redy(lkb,l,lkc)/(2.d0*l+1.d0)

      do 21 ik1=1,2
      k1=lka-3+2*ik1
      if(k1.lt.0) go to 21
      slk1=float(k1)
      call troisj(slk1,slka,1.d0,0.d0,0.d0,0.d0,c1)
      cofk1=dsqrt(2.d0*k1+1.d0)*c1

      do 22 ik2=1,2
      k2=lkd-3+2*ik2
      if(k2.lt.0) go to 22
      slk2=float(k2)
      call troisj(slk2,slkd,1.d0,0.d0,0.d0,0.d0,c2)
      cofk2=dsqrt(2.d0*k2+1.d0)*c2
      call sixj(slk2,slk1,sl,slka,slkd,1.d0,c3)
      cofk12=redy(k1,l,k2)*c3
      vv(iup,il)=vv(iup,il)+xlkad*cofl*cofk1*cofk2*cofk12*rx2  ! II rx2
     &(ik1,ik2)
   22 continue
   21 continue
   20 continue
   18 continue
      do 23 il=ilpt,ilgd
         l=jt-2+il
         if(l.lt.0) go to 23
         sl=float(l)
         dlp1=2.d0*l+1.d0
         call neufj(sja,sla,0.5d0,sjd,sld,0.5d0,sjt,sl,1.d0,c3)
         call neufj(sjc,slc,0.5d0,sjb,slb,0.5d0,sjt,sl,1.d0,c4)
         c4=c3*c4
            do 24 iup=1,7
               hh(iup)=hh(iup)+cofs*dlp1*c4*vv(iup,il)
   24       continue
   23 continue
      do iup=1,7
         cbet(iup)=elias*cbets(iup)+eliav*cbetv(iup)
         calf(iup)=elias*calfs(iup)+eliav*calfv(iup)  
         hh(iup)=cbet(iup)*hh(iup)+calf(iup)*cof0*vv(iup,2)
      enddo
   32 continue
C_TEST
c     hh(9)=0.d0
      do iup=1,9
         res=res+hh(iup)             ! Here the result !!
         if(ia.eq.69.and.id.eq.1.and.ic.eq.39.and.ib.eq.7)then
          write(2,*) 'iup,hh=',iup,hh(iup)
c         if(iup.eq.9)stop
         endif 
      enddo
C_TEST
c     res=hh(8)
c     res=hh(8)+hh(9)

c1000 continue
      Iexp=(-la-lb+lc+ld)/2          ! Phase Convention II is assumed !!
      res = Ipot(Iexp)*res

      return
      end
    
C**********************************************************************

      SUBROUTINE RDINT1(i_t3,IA,IB,IC,ID,NR,DR)

c     produce ri0,ris,rx1

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)

      ri0=0.D0
      ris=0.D0
      rx1=0.d0
      xl=ll(ia)*(ll(ia)+1)+ll(ib)*(ll(ib)+1)+ll(ic)*(ll(ic)+1)
      xl=xl+ll(id)*(ll(id)+1) 
      do 1 j=1,nr
         x1=j*dr
         x2=x1*x1
         c=wf(j,ia)*wf(j,ib)*wf(j,ic)*wf(j,id)*dr/x2
         ri0=ri0 + float(i_t3)*f0(j,3)*c + float(1-i_t3)*f1(3)*c           
         ris=ris + float(i_t3)*g0(j,3)*c + float(1-i_t3)*g1(3)*c
         c1=-xl*c/x2+2.d0*dr*((dwf(j,ia)*wf(j,ib)+wf(j,ia)*dwf(j,ib))
     &   *wf(j,ic)*wf(j,id)+(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id))*
     &   wf(j,ia)*wf(j,ib))/(x1*x2)
         c1=c1-2.d0*dr*(dwf(j,ia)*(dwf(j,ib)*wf(j,ic)*wf(j,id)+wf(j,
     &   ib)*(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id)))+wf(j,ia)*(dwf(
     &   j,ib)*(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id))+wf(j,ib)*dwf(
     &   j,ic)*dwf(j,id)))/x2
         rx1=rx1+c1      
    1 continue  
                                                                               
      return
      end

C**********************************************************************
 
      SUBROUTINE RDINT2(IA,IB,IC,ID,NR,DR)

c     produce rx2(*,*)

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      dimension ua(2),ub(2)

      common/brint/ri0,ris,rx1,rx2(2,2)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)
      common/bfg/f0(nnn,3),g0(nnn,3),corr_res(nnn),f1(3),g1(3)


      do 2 k1=1,2
         do 6 k2=1,2
            rx2(k1,k2)=0.d0
    6    continue
    2 continue
      do 1 j=1,nr
         x1=j*dr
         x2=x1*x1
         do 3 k=1,2
            sk=3.d0-2.d0*k
            ua(k)=dwf(j,ia)+sk*(ll(ia)+k-1)*wf(j,ia)/x1
            ub(k)=dwf(j,ib)+sk*(ll(ib)+k-1)*wf(j,ib)/x1
    3    continue
         do 4 k1=1,2
            do 5 k2=1,2
               c23=ua(k1)*ub(k2)*wf(j,ic)*wf(j,id)
               rx2(k1,k2)=rx2(k1,k2)+c23*dr/x2
    5       continue
    4    continue
    1 continue  
      return

      end

C**********************************************************************

      FUNCTION REDY(L1,L2,L3)


      implicit double precision (a-h,o-z)
c     real*4 x1,x2,x3,xm,c3j

      redy=0.d0
      l13m=iabs(l1-l3)
      l13p=l1+l3
      if((l2-l13m)*(l2-l13p).gt.0) go to 1
      l123=l13p+l2
      l123=mod(l123,2)
      if(l123.ne.0) go to 1
      qp=4.d0*3.14159265d0
      x4=(2.d0*l1+1.d0)*(2.d0*l2+1.d0)*(2.d0*l3+1.d0)
      x4=x4/qp
      sg=1.d0-2.d0*mod(l1,2)
      x1=float(l1)
      x2=float(l2)
      x3=float(l3)
      xm=0.
      call troisj(x1,x2,x3,xm,xm,xm,c3j)
      redy=sg*dsqrt(x4)*c3j
    1 return

      end

C**********************************************************************

      FUNCTION TJL(LA,JA,LB,JB,J,L) 

ccc   ja and jb are twice the true values;
ccc   other variables have their true values

      implicit double precision (a-h,o-z)

      tjl=0.D0
      z=0.D0  
      qp=4.D0*3.14159265D0  
      fj=float(j)   
      fl=float(l)   
      jd=2*j
      lp=la+lb  
      lm=iabs(la-lb)
      ll=(l-lp)*(l-lm)  
      if(ll.gt.0) go to 1   
      ll=la+lb+l-2*((la+lb+l)/2)
      if(ll.ne.0) go to 1   
      jp=(ja+jb)/2  
      jm=iabs(ja-jb)/2  
      jj=(j-jp)*(j-jm)  
      if(jj.gt.0) go to 1   
      ljp=l+j   
      ljm=iabs(l-j) 
      lj=(1-ljp)*(1-ljm)
      if(lj.gt.0) go to 1   
      if(l.ne.j) go to 2
      ig=(ja+jb)/2+j
      cc=(jd+1.D0)/(fj*(fj+1.D0))   
      z=-0.5D0*dsqrt(cc)*(jb+1.D0+(1.D0-2.D0*mod(ig,2))*(ja
     &+1.D0)) 
      go to 3   
    2 ig=lb+(jb+1)/2
      cc=(la-0.5D0*ja)*(ja+1.D0)+(lb-0.5D0*jb)*(jb+1.D0)
      sg=1.D0-2.D0*mod(ig,2)
      if(l.eq.(j+1)) z=sg*(cc+fl)/dsqrt(fl)  
      if(l.eq.(j-1)) z=sg*(cc-fl-1.d0)/dsqrt(fl+1.d0)
    3 c=(ja+1.D0)*(jb+1.d0)/((jd+1.d0)*qp)
      xja=ja/2.d0
      xjb=jb/2.d0
      xjd=j
      jpha=(ja+jb)/2+1
      jpha=1-2*mod(jpha,2)  
      tjl=(1.d0-2.d0*mod(la,2))*dsqrt(c)*z*jpha
     **cofcg(xja,xjb,xjd,0.5d0,-0.5d0,0.d0)
    1 return

      end   

C**********************************************************************

      FUNCTION RMAT_M1(LA,JA,LB,JB,IQ) 

ccc   ja and jb are twice the true values;
ccc   other variables have their true values

      implicit double precision (a-h,o-z)
      data glp,gln,gsp,gsn /1.d0,0.d0,5.586,-3.826/

      rmat_m1=0.d0
      if(la.ne.lb)go to 1
      xja=float(ja)/2.d0
      xjb=float(jb)/2.d0
      if(iq.eq.0)then
       gl=gln
       gs=gsn
      else 
       gl=glp
       gs=gsp
      end if
      if(ja.ne.jb)go to 2
      prefac=dsqrt((xja+1.d0)/xja)
      isign=(xja-float(la))/dabs(xja-float(la))
      schmidt=xja*(gl+isign*(gs-gl)/float(2*la+1))
      rmat_m1=prefac*schmidt
      go to 1
    2 continue
      rmat_m1=(gs-gl)*sqrt(float(la+1)/float(2*la+1))
    1 return
      end       

C**********************************************************************

      SUBROUTINE VPPT3(Ia,Ib,Ic,Id,iphase,Jt,nr,dr,res)
c     iphase = 1 qq; iphase = -1 qq'.

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1calfv(7),cbetv(7) 
      common/brint3/ris3
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      Common /Basis1/focc(nsp),delta(nsp),spe(nsp)

c     write(49,*) 'x3 = ',x3

C ------ Matrix Element calculation --------------------------------------------
      call rdint3(ia,ib,ic,id,nr,dr)        
c     write(49,*) 'rdint3 = ',ris3
      PI = 3.14159265358d0          ! Pi for Pitagora...
      djp1=2.d0*jt+1.d0
      la=ll(ia)
      lb=ll(ib)
      lc=ll(ic)
      ld=ll(id)
      ja=lj(ia)
      jb=lj(ib)
      jc=lj(ic)
      jd=lj(id)
      xja= float(ja)/2.d0
      xjb= float(jb)/2.d0
      xjc= float(jc)/2.d0
      xjd= float(jd)/2.d0
      xjt= float(jt)
     
      c1 = -(1.d0/12.d0)*t3*ris3/(djp1*4*pi)
c     write(49,*) c1
      c2= sqrt(float(ja+1))*sqrt(float(jb+1))
     &           *sqrt(float(jc+1))*sqrt(float(jd+1))
      c1=c1*c2
c     write(49,*) ja,jb,jc,jd,c2
c     write(49,*) c2

      if(iphase.eq.1)then 

       iexp_qq = (ja-jc)/2 + lb + ld + 1
       iphase_qq = ipot(iexp_qq)
       iphase_qq_2 = ipot(la+lb+jt)
       res = iphase_qq*cofcg(xja,xjb,xjt,0.5d0,-0.5d0,0.d0)*
     & cofcg(xjc,xjd,xjt,0.5d0,-0.5d0,0.d0)*(iphase_qq_2 + 1)*
     & (1-x3)
c      write(49,*) res

      else if(iphase.eq.-1)then
       
       iexp_qqp = (jb-jd)/2 + lb + ld + 1
       iphase_qqp = ipot(iexp_qqp)
       iphase_qqp_2 = ipot(la+lb+jt)
       iphase_j = ipot((ja+jb+jc+jd)/2)
       res = ((1+x3)*cofcg(xja,xjb,xjt,0.5d0,0.5d0,1.d0)*
     & cofcg(xjc,xjd,xjt,0.5d0,0.5d0,1.d0) + iphase_qqp*
     & cofcg(xja,xjb,xjt,0.5d0,-0.5d0,0.d0)*
     & cofcg(xjc,xjd,xjt,0.5d0,-0.5d0,0.d0)*(iphase_qqp_2 - x3))
     & *iphase_j

      endif

      res= res * c1
c     write(49,*) res

c     cc1= -(1.d0/24.d0)*t3*ris3*iphase/djp1
c     
c     res1= 1.d0/2.d0 * float(1+isgn) *
c    &   ((1.d0+x3)*tjl(la,ja,lb,jb,jt,jt) * tjl(lc,jc,ld,jd,jt,jt) 
c    &   - 3.d0*(1.d0-x3)*yl(la,ja,lb,jb,jt) * yl(lc,jc,ld,jd,jt))
c    &   + 1.d0/2.d0 * float(1-isgn) *
c    &   (1.d0+x3)*(tjl(la,ja,lb,jb,jt,jt-1) * tjl(lc,jc,ld,jd,jt,jt-1) 
c    &   + tjl(la,ja,lb,jb,jt,jt+1) * tjl(lc,jc,ld,jd,jt,jt+1))
c
c     res1= res1 * cc1
c
c     if (dabs(res1-res).gt.0.00001) then
c       STOP
c     endif
       
      Iiexp=(-la-lb+lc+ld)/2          
c     Phase Convention II is assumed !!
      res = Ipot(Iiexp)*res
c     res1 = Ipot(Iiexp)*res1

c     write(49,*) res

      return
      end


C**********************************************************************

      SUBROUTINE RDINT3(IA,IB,IC,ID,NR,DR)

c     produce ris3

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      common/brint3/ris3
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      COMMON /Sk/T0,T1,T2,T3,X0,X1,X2,X3,W,alfe
      common/betalf/calf(7),cbet(7),calfs(7),cbets(7),
     1calfv(7),cbetv(7) 
      COMMON/DENS/dt(nnn),dn(nnn),dp(nnn)

      ris3=0.D0
      do 1 j=1,nr
         x1=j*dr
         x2=x1*x1
         ris3=ris3 + wf(j,ia)*wf(j,ib)*wf(j,ic)*wf(j,id)
     &          *(dt(j))**(alfe) * dr/x2
    1 continue  
                                                                               
      return
      end

C**********************************************************************

      subroutine qwg_nonag(b,ns,ls,js,iq) 
      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      dimension u(nnn,nmx),up(nnn,nmx),d(5,5),vu(nnn),wu(nnn)   
      dimension a(nmx,nmx),v(nmx)
c     ,tt1(nmx),tt2(nmx),tt3(nmx)
      dimension e(nmx),vec(nmx),nord(nmx)
      common/bwg/wg(nnn,nmx),ewg(nmx),xor(nmx),dwg(nnn,nmx)   
      common/hf1/del
      common/hf/nmax,nocc,nunocc,norb
      common/bpot/xmn(nnn),xmp(nnn),vvn(nnn),vvp(nnn),yn(nnn),yp(nnn),
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)   
      data d(1,1),d(1,2),d(1,3),d(1,4),d(1,5)/-50.d0,96.d0,-72.d0,
     132.d0,-6.d0/
      data d(2,1),d(2,2),d(2,3),d(2,4),d(2,5)/-6.d0,-20.d0,36.d0,
     1-12.d0,2.d0/ 
      data d(3,1),d(3,2),d(3,3),d(3,4),d(3,5)/2.d0,-16.d0,0.d0,
     116.d0,-2.d0/   
      data d(4,1),d(4,2),d(4,3),d(4,4),d(4,5)/-2.d0,12.d0,-36.d0,
     120.d0,6.d0/  
      data d(5,1),d(5,2),d(5,3),d(5,4),d(5,5)/6.d0,-32.d0,72.d0,
     1-96.d0,50.d0/ 

      if(ns.gt.nmx)then
       write(2,*) ' '
       write(2,*) '>>> PARAMETER NMX MUST BE INCREASED '
       stop
      end if

      h = del
      cof = 1.d0 
      nri = nmax
      dr = h 

      eps=(1.0d-015)*ns*ns  
      nn2=nri-2 
      l=ls  
      b3=b*b*b*1.772454d0  
      db3=b3
      do 1 n1=1,ns   
      n=n1-1
      do 2 j=1,nri   
      r=h*j 
      x=r/b 
      dx=x 
      xsq=-2.d0*x*x   
      y=-x*x/2.d0 
      dy=y
      dxsq=xsq
      dbr=r
      bh=24.d0*h
      dbh=bh
      rho=1.d0
      vnl=1.d0
      if(n) 101,101,102 
  102 xa=1.d0 
      k=2*l+1   
      nss=n1
      ni=0  
      do 103 ii=1,n  
      k=k+2 
      ni=ni+1   
      nss=nss-1 
      sni=float(nss)/(ni*float(k))
      dni=sni
      xa=xa*dni*dxsq
      sk=float(k)/(2.d0*ni)
      dk=sk
      vnl=vnl+xa
      rho=dk*rho
  103 continue
  101 dpl=1.d0
      xsl=1.d0
      if(l) 104,104,105 
  105 k=1   
      do 106 ii=1,l  
      k=k+2 
      xk=k  
      dk=xk
      dpl=2.d0*dpl/dk
      xsl=xsl*dx
  106 continue 
  104 arg=rho*dpl/db3
      xa=2.d0*xsl*dbr*dsqrt(arg) 
      u(j,n1)=xa*vnl*dexp(dy) 
    2 continue  
      do 3 j=1,nri   
      k=3   
      if(j.lt.3) k=j
      if(j.gt.nn2) k=j-nri+5
      sum=0.d0
      do 4 ii=1,5
      jj=j+ii-k 
      dbd=d(k,ii)
      sum=sum+dbd*u(jj,n1)
    4 continue  
      up(j,n1)=sum/(dbh)
    3 continue  
    1 continue  
      xll=l*(l+1.d0)  
      xjj=0.5d0*(0.25d0*js*(js+2.d0)-xll-0.75d0)  
      go to (5,6),iq
    5 do 7 j=1,nri   
      x2=h*j*h*j
      wu(j)=xmn(j)  
      vu(j)=xmn(j)*xll/x2+vvn(j)+xjj*vsn(j)
    7 continue  
      go to 8   
    6 do 9 j=1,nri   
      x2=h*j*h*j
      wu(j)=xmp(j)  
      vu(j)=xmp(j)*xll/x2+vvp(j)+xjj*vsp(j)+vc(j)
    9 continue
    8 do 10 i1=1,ns  
      do 11 i2=i1,ns 
      aa=0.d0 
      do 12 j=1,nri  
      dwu=wu(j)
      dvu=vu(j)
      aa=aa+dwu*up(j,i1)*up(j,i2)+dvu*u(j,i1)*u(j,i2)
   12 continue  
      aa=aa*h
      a(i1,i2)=(aa) 
      a(i2,i1)=a(i1,i2) 
   11 continue  
   10 continue  
c 110 format(18(1x,e12.5))
      call diasym(a,nmx,ns,v,e,vec,nord)

      do 20 n=1,ns   
      xnorm=0.d0  
      sgn=0.d0
      dy=0.d0
      ewg(n)=v(n) 
      do 50 ik=1,ns  
      dy=dy+a(ik,n)**2
      as=a(ik,n)  
      sgn=sgn+as*u(4,ik)
   50 continue
      xor(n)=dy
      sign=1.d0   
      if(sgn.lt.0.d0) sign=-1.d0
      do 21 j=1,nri  
      xx=0.d0 
      xyx=0.d0
      do 22 ii=1,ns  
      xyx=xyx+a(ii,n)*up(j,ii)
      xx=xx+a(ii,n)*u(j,ii)
   22 continue  
      wg(j,n)=xx*sign   
      dwg(j,n)=sign*xyx
      xnorm=xnorm+(wg(j,n)**2)*h
   21 continue
      xnorm=dsqrt(xnorm) 
      do 23 j=1,nri  
      dwg(j,n)=dwg(j,n)/xnorm
      wg(j,n)=wg(j,n)/xnorm
   23 continue 
   20 continue  
      return
      end
 
      SUBROUTINE DIASYM(A,NP,N,D,E,VEC,NORD)
ccc   diagonalization of a real, symmetric matrix A
ccc   cfr. Numerical Recipes, for this routine and related ones
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(NP,NP),D(NP),E(NP),NORD(NP),VEC(NP)
      CALL TRED2(A,N,NP,D,E)
      CALL TQLI(D,E,N,NP,A)
      CALL REORDERING(D,A,N,NP,VEC,NORD)
      RETURN
      END

      SUBROUTINE REORDERING(E,A,N,NP,VEC,NORD)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION NORD(NP),E(NP),VEC(NP),A(NP,NP)
CCC   ORDERING OF STATES
      DO 224 I=1,N
      NORD(I)=I
  224 continue 
      DO 225 I=1,N
      K=I   
      X=E(I)  
      DO 226 J=1,N  
      VEC(J)=A(J,I)
  226 continue  
      DO 227 II=I,N
      IF(X-E(II)) 227,227,228 
  228 K=II  
      X=E(II) 
      DO 229 J=1,N  
      VEC(J)=A(J,II)
  229 continue 
  227 CONTINUE  
      KK=NORD(K)
      NORD(K)=NORD(I)   
      NORD(I)=KK
      E(K)=E(I) 
      DO 230 J=1,N  
      A(J,K)=A(J,I)
  230 continue   
      E(I)=X  
      DO 231 J=1,N  
      A(J,I)=VEC(J)
  231 continue  
  225 CONTINUE  
      RETURN
      END

      function yll(l1,l2,lambda)
      implicit double precision (a-h,o-z)
      yll=0.d0
      sl1=float(l1)
      sl2=float(l2)
      slambda=float(lambda)
      call troisj(sl1,slambda,sl2,0.d0,0.d0,0.d0,c)
      if(dabs(c).le.1e-6)go to 1
      yll=dsqrt(2.d0*float(l2)+1.d0)*ipot(l1)*c
      fact=dsqrt((2.d0*float(lambda)+1.d0)/16.d0/datan(1.d0))
      yll=yll*fact
    1 return
      end

C**********************************************************************

      SUBROUTINE subphme_Coulomb(Ia,Ib,Ic,Id,Jt,nr,dr,res)

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bii/ia1,ib1,ic1,id1
      common/bjj/jtc
      common/bnn/xxmin,xxmax
      data e2 /1.43986/
      external fintc,phi1,phi2

      jtc=jt
c     write(23,*) 'jtc=',jtc
      res = 0.d0
      if(iq(ia).ne.1.or.iq(ib).ne.1.or.iq(ic).ne.1.or.iq(id).ne.1)
     &go to 1
      la=ll(ia)
      lb=ll(ib)
      lc=ll(ic)
      ld=ll(id)
      ja=lj(ia)
      jb=lj(ib)
      jc=lj(ic)
      jd=lj(id)

      doubleint=0.d0

      DO 30 I1=1,NR
      R1=DR*FLOAT(I1)
      DO 40 I2=1,NR
      R2=DR*FLOAT(I2)
      rmax=max(r1,r2)
      rmin=min(r1,r2)
      PS4=wf(i1,ia)*wf(i2,ib)*wf(i2,ic)*wf(i1,id)
      doubleint=doubleint+PS4*rmin**jt/rmax**(jt+1)
  40  CONTINUE 
  30  CONTINUE 
      doubleint=doubleint*dr*dr

c     xxmin=dr
c     xxmax=dr*float(nr-1)
c     absac=1.0d-4 
c     ifail=0
c     ia1=ia
c     ib1=ib
c     ic1=ic
c     id1=id
c     call d01daf(xxmin,xxmax,phi1,phi2,fintc,absac,doubleint,
c    &            npout,ifail)      

      res=e2*16.d0*atan(1.d0)*doubleint*
     &yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/float(2*jt+1)**2
      Iexp=(-la-lb+lc+ld)/2    ! Phase Convention II is assumed !!
      res=Ipot(Iexp)*res
c     write(95,31) ia,ib,ic,id,yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)
c    &,doubleint,res
c  31 format(4(i3),3(2x,e12.5))
    1 return
      end 

C**********************************************************************

      function fintc(r1,r2)
      implicit real*8(a-h,o-z)
      include 'param.qrpa'
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bii/ia1,ib1,ic1,id1
      common/bjj/jtc
      external wf_int

c     write(23,*) 'In fintc: r1,r2 = ',r1,r2

      rmax=max(r1,r2)
      rmin=min(r1,r2)
      
      ps4=wf_int(r1,r2)
      produ=ps4*rmin**jtc/rmax**(jtc+1)
      fintc=produ
c     write(23,*) 'In fintc: ',ps4,r1,r2,jtc 
c     write(23,*) 'In fintc: ',produ
c     fintc=rmin**jtc/rmax**(jtc+1)*dexp(-r1**2/5.d0)*dexp(-r2**2/5.d0)
c     fintc=1.d0

      return
      end

C**********************************************************************

      function phi1(y)
      implicit real*8(a-h,o-z)
      common/bnn/xxmin,xxmax

      dummy=y
      phi1=xxmin

      return
      end

C**********************************************************************

      function phi2(y)
      implicit real*8(a-h,o-z)
      common/bnn/xxmin,xxmax

      dummy=y
      phi2=xxmax

      return
      end

C**********************************************************************
 
      function wf_int(x,y)
      implicit real*8 (a-h,o-z)
      include 'param.qrpa'
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common /hf1/del
      common/bii/ia1,ib1,ic1,id1
      
c     write(2,*) 'In wf_int: '
c     write(2,*) ia1,ib1,ic1,id1
c     write(2,*) x,y
      dx=del
      ind_x=(x-dx)/dx+1.0001
      resto_x=(dx*float(ind_x+1)-x)/dx
      ind_y=(y-dx)/dx+1.0001
      resto_y=(dx*float(ind_y+1)-y)/dx
      wf_int1=wf(ind_x,ia1)*resto_x+wf(ind_x+1,ia1)*(1.d0-resto_x)
      wf_int2=wf(ind_x,id1)*resto_x+wf(ind_x+1,id1)*(1.d0-resto_x)
      wf_int3=wf(ind_y,ic1)*resto_y+wf(ind_y+1,ic1)*(1.d0-resto_y)
      wf_int4=wf(ind_y,ib1)*resto_y+wf(ind_y+1,ib1)*(1.d0-resto_y)
c     if(ia1.eq.1) then
c      write(25,*) x,ind_x,resto_x
c      write(25,*) wf(ind_x,ia1),wf(ind_x+1,ia1),wf_int1 
c      write(24,*) x,wf_int1
c     end if
c     if(id1.eq.1) write(24,*) x,wf_int2
c     if(ic1.eq.1) write(24,*) y,wf_int3
c     if(ib1.eq.1) write(24,*) y,wf_int4
      wf_int=wf_int1*wf_int2*wf_int3*wf_int4
    1 return
      end

C**********************************************************************

      SUBROUTINE subphme_Coulomb_exchange(Ia,Ib,Ic,Id,Jt,nr,dr,res)

      implicit double precision (a-h,o-z)
      include 'param.qrpa'
      Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),
     &               iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      COMMON /DENS/ dt(nnn),dn(nnn),dp(nnn)
      data e2 /1.43986/

      res = 0.d0

      if(iq(ia).ne.1.or.iq(ib).ne.1.or.iq(ic).ne.1.or.iq(id).ne.1)
     &go to 1

      pi=4.d0*datan(1.d0)
      c=-(e2/3.d0)*(3.d0/pi)**(1.d0/3.d0)

      la=ll(ia)
      lb=ll(ib)
      lc=ll(ic)
      ld=ll(id)
      ja=lj(ia)
      jb=lj(ib)
      jc=lj(ic)
      jd=lj(id)
 
      djp1=2.d0*jt+1.d0

      ri=0.d0
      DO 40 I1=1,NR
      R2=DR*FLOAT(I1)
      r2=r2*r2
      PS4=wf(i1,ia)*wf(i1,ib)*wf(i1,ic)*wf(i1,id)
      ri=ri+dr*PS4/r2*dp(i1)**(-2.d0/3.d0)
  40  CONTINUE 
c     write(23,*) 'ri=',ri
      ri=ri*c      
c     write(23,*) 'c=',c

      res=ri*yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/djp1
c     write(23,*) yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt),djp1

      Iexp=(-la-lb+lc+ld)/2    ! Phase Convention II is assumed !!
      res=Ipot(Iexp)*res
c     write(95,31) ia,ib,ic,id,yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)
c    &,doubleint,res
c  31 format(4(i3),3(2x,e12.5))
c     write(23,*) Iexp,res
    1 return
      end 

C**************************************************************************

        SUBROUTINE SUBPHME_SO(J,ia,ib,ic,id,SOphme)

c       phme <(ad)J|v_so|(cb)J> is obtained from ppmes <(ab)J'|v_so|(dc)J'> 
c       via Pandya relation. Chex case:
c       <(p1,n1)J|v_so|(p2,n2)J>  <---  <(p1,n2)J'|v_so|(n1,p2)J'> 

        include 'param.qrpa'                                       ! nsp
        implicit real*8 (a-h,o-z)
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),           ! ll,lj,iq
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)   
        common /bres_so/ w2,w2p,c_so,ri_so 
          
cp        write(*,*) 'check entro in subphme_so'
cp        write(*,*) J,ia,ib,ic,id,VphJ

        ita=1-2*iq(ia)
        itb=1-2*iq(ib)
        itc=1-2*iq(ic)
        itd=1-2*iq(id)

        sJ=dfloat(J)

        la=ll(ia)
        lb=ll(ib)
        ja=lj(ia)                  !twice the true value
        jb=lj(ib)
 
        sla=dfloat(la)
        slb=dfloat(lb)
        sja=dfloat(ja)*0.5d0       !real value
        sjb=dfloat(jb)*0.5d0   

        lc=ll(ic)
        ld=ll(id)
        jc=lj(ic)
        jd=lj(id)     

        slc=dfloat(lc)
        sld=dfloat(ld)
        sjc=dfloat(jc)*0.5d0
        sjd=dfloat(jd)*0.5d0

        JTprime_min=max(dabs(sja-sjb),dabs(sjc-sjd))
        JTprime_max=min(sja+sjb,sjc+sjd)

cp        write(*,*)'dabs(sja-sjb),dabs(sjc-sjd),sja+sjb,sjc+sjd
cp     &  (J1)extr='
cp        write(*,*)dabs(sja-sjb),dabs(sjc-sjd),sja+sjb,sjc+sjd
cp        write(*,*)JTprime_min, JTprime_max
cp        write(*,*)'_____________________________________________'

        VphJ=0.0d0
        do J1 =JTprime_min,JTprime_max
           sJ1=dfloat(J1)

           VppJ1=0.0d0
           c6j=0.0d0
           call subppme_so(J1,ia,ib,id,ic,VppJ1)
           call sixj(sja,sjb,sJ1,sjc,sjd,sJ,c6j)
cp           write(*,*)'sja,sjb,sJ1,sjc,sjd,sJ,c6j'
cp           write(*,*)sja,sjb,sJ1,sjc,sjd,sJ,c6j
           VphJ =VphJ+(2*J1+1)*C6j*ipot((jc+jd)/2-J1)*VppJ1   ! Pandya

cp           write(*,*) J1,2*J1+1
cp           write(*,*) C6j,ipot((jc+jd)/2-J1)
cp           write(*,*) VppJ1
cp           write(*,*) 'parz', VphJ
cp           write(*,*) 


        enddo

C---Isospin Matrix Elements

        if((ita.eq.itd).and.(itb.eq.itc)) then   ! no chex  
         cis=1.d0                               
         civ=float(ita*itb)                             ! (pp)(nn)           
        else 
         cis=0.d0                              ! prep. chex
         civ=0.d0
        endif                                                                

        if(((ita+itd).eq.0).and.((itb+itc).eq.0)) then     ! chex   (pn)(pn)
         civ=civ+1.d0+float(ita*itc)                   !         del.(pn)(np)
        endif

C---
        SOphme=0.5d0*(2.d0*w2+w2p)*cis*VphJ+0.5d0*w2p*civ*VphJ

        iexp=(-la-lb+lc+ld)/2          ! Phase Convention II is assumed !!
        SOphme=Ipot(Iexp)*SOphme

        RETURN
        END


C*******************************************************************************************


        SUBROUTINE SUBPPME_SO(JT,Ia,Ib,Ic,Id,SOppme)


c       Calculate two body Spin-Orbit particle-particle matrix element
c       <(a1,b2)J|v_so|(c1,d2)J|>.
c       The interaction is in the form (cfr. D.Vautherin and D.M.Brink PRC5(1972),626)
c       v_so = W0*i/4(sigma1+sigma2)[(grad1'-grad2') x delta(1-2)(grad1-grad2)].
c       The calculation is performed separating spin dependent part and orbital dependent part
c       of the interaction.
c
c       (to obtain two body Spin-Orbit particle-hole matrix elements via 
c        inverse Pandya relation)

        include 'param.qrpa'                                            ! nsp
        implicit real*8 (a-h,o-z)
c       common/rad/nmax,dstep
        Common /Basis/ lev(nsp),nn(nsp),ll(nsp),lj(nsp),                ! ll,lj
     &                 iq(nsp),JJ(ncf),It(ncf),ipp(ncf),inn(ncf)
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del                                                 ! del
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
        common /bres_so/ w2,w2p,c_so,ri_so 


c       open(unit=120,status='unknown',name='SOppme.dat') 


        w0=w2
        w0p=w2p
        dstep=del

c        write(2,*)'W0 in the subroutine=',W0
c        write(*,*)'del=',del

cp        write(*,*)'----------------------------------------------------'
cp        write(*,*)'entering  in SO_PPMATRIXELEMENTS'
cp        write(*,*)JT,Ia,Ib,Ic,Id,soppme
     
        sJT=dfloat(Jt)
 
        la=ll(ia)
        lb=ll(ib)
        ja=lj(ia)                  !twice the true value
        jb=lj(ib)
 
        sla=dfloat(la)
        slb=dfloat(lb)
        sja=dfloat(ja)*0.5d0       !real value
        sjb=dfloat(jb)*0.5d0   

        lc=ll(ic)
        ld=ll(id)
        jc=lj(ic)
        jd=lj(id)     

        slc=dfloat(lc)
        sld=dfloat(ld)
        sjc=dfloat(jc)*0.5d0
        sjd=dfloat(jd)*0.5d0

cp      write(*,*)sla,slb,slc,sld
cp      write(*,*)sja,sjb,sjc,sjd
 
c-----calculating the only total Spin dependent terms----------------------

                                                   ! spin reduced matrix element 
        depS1S2 = -1.0d0*3.0d0*dsqrt(6.0d0)*2.0d0  ! is non vanishing only for S1=S2=1;
                                                   ! working with S=sigma/2.

c-----calculating the total orbital angular momentum dependent terms---(setted S1=S2=1)---
        

c       write(1002,*)'calculating depL1L2:'

        depL1L2 = 0.0d0                      
        L1min = max(iabs(la-lb),iabs(1-JT))
        L1max = min(la+lb,1+JT)

c       if(L1min.gt.L1max) goto 22       ! cmq autom.
c       if(L1min.gt.L1max) then
c    &  write(1002,*)'triangular conditions on L1 failing',L1min,L1max
c       goto 22  
c       endif      

        do L1 = L1min,L1max     ! final
           sL1 = dfloat(L1)
   
           L2min = max(iabs(lc-ld),iabs(1-JT),iabs(L1-1))
           L2max = min(lc+ld,1+JT,L1+1)

c          if(L2min.gt.L2max) goto 22  
c          if(L2min.gt.L2max) then
c    &     write(1002,*)'triangular conditions on L2 failing',L2min,L2max
c          goto 22
c          endif      
        
           do L2 = L2min,L2max    ! initial
              sL2 = dfloat(L2)

              C9j_ab = 0.0d0      !initializing
              C9j_cd = 0.0d0
              C6j = 0.0d0
              reduxL1L2 = 0.0d0

              call neufJ(sla,0.5d0,sja,slb,0.5d0,sjb,
     &                                 sL1,1.0d0,sJT,C9j_ab) 

              call neufJ(slc,0.5d0,sjc,sld,0.5d0,sjd,
     &                                 sL2,1.0d0,sJT,C9j_cd)

              call sixJ(sL2,1.0d0,sJT,1.0d0,sL1,1.0d0,C6j)

              call redmat_orb(ia,ib,ic,id,la,lb,lc,ld,
     &                                      L1,L2,reduxL1L2)


              depL1L2 = depL1L2 + dsqrt(2.0d0*L1+1.0d0)*
     &                   dsqrt(2.0d0*L2+1.0d0)*
     &                   ipot(L2)* C9j_ab* C9j_cd* C6j*    
     &                   reduxL1L2
              

           enddo
        enddo


c-----assembling the expression--------------------------------------------


        SOppme =  1.d0/(4.0d0)*ipot(JT)*    
     &            dsqrt(ja+1.0d0)*dsqrt(jb+1.0d0)* 
     &            dsqrt(jc+1.0d0)*dsqrt(jd+1.0d0)*
     &            depS1S2 * depL1L2
            
  
cp        write(*,*) 'soppme',JT,soppme

        RETURN
        END  
        

C--------------------------------------------------------------------------
C                              subSUBROUTINES
C--------------------------------------------------------------------------

C   Oss. non ho potuto isolare azione grad. Pensare per event altra struttura


        SUBROUTINE REDMAT_ORB(ia,ib,ic,id,la,lb,lc,ld,Lf,Li,reduxLfLi)

        implicit real*8 (a-h,o-z)


c       write(1002,*)'Calculating orbital reduced matrix element'       

        pi=4.0d0*datan(1.0d0)
  
        do i=-1,1,2
           do j=-1,1,2

           C_abcd = 0.0d0
           C_badc = 0.0d0
           S_ac = 0.0d0
           S_bd = 0.0d0

           call subCoupl(i,j,Lf,Li,la,lb,lc,ld,C_abcd)
           call subCoupl(i,j,Lf,Li,lb,la,ld,lc,C_badc)
           call subIntRad(i,j,ia,ib,ic,id,la,lc,S_ac)
           call subIntRad(i,j,ib,ia,id,ic,lb,ld,S_bd) 

c          write(1002,*)'C_abcd,S_ac', C_abcd,S_ac
c          write(1002,*)'C_badc,S_bd', C_badc,S_bd


           reduxLfLi = reduxLfLi + dsqrt(6.0d0)/(4.0d0*pi)*   
     &                 dsqrt(2.0d0*Lf+1.0d0)*dsqrt(2.0d0*Li+1.0d0)*
     &                 ipot(Lf+Li)*ipot(1+(i+j)/2)*
     &                 2.0d0*(C_abcd*S_ac + ipot(Lf+Li)* C_badc*S_bd)

      
           enddo
        enddo

         
        RETURN
        END


C..........................................................................


        SUBROUTINE subCOUPL(i,j,Lf,Li,I1,I2,I3,I4,C)
                        ! pensata per lavorare con le quantita' che servono, cioe' intere

        implicit real*8 (a-h,o-z)


c       write(1002,*)'entering subCOUPL'

        si=dfloat(i)
        sj=dfloat(j)

        sLf=dfloat(Lf)
        sLi=dfloat(Li)            !gia' calc prima. Magari inserire condiz. solo nelle 
                                  !rout 9C, 3J , 6J (se già non ci sono in quelle nuove) 
                                  !per lavorare con quantità reali, senza doverselo ricordare a
                                  !priori

        sI1=dfloat(I1)
        sI2=dfloat(I2)
        sI3=dfloat(I3)
        sI4=dfloat(I4)

        s1=dfloat(I1+i)
        s3=dfloat(I3+j)

        goto 11   

      
        kmin=iabs(I1-I3)
        kmax=I1+I3

        do k=kmin,kmax      !``s'' 
           sk=dfloat(k)
            
           lmin=max(iabs(k-1),iabs((I1+i)-(I3+j)),iabs(I2-I4))
           lmax=min(k+1,(I1+i)+(I3+j),I2+I4)
            
           do l=lmin,lmax   !``lambda''
              sl=dfloat(l)

              C9j_1 = 0.0d0
              C9j_2 = 0.0d0             
              C3j_1 = 0.0d0
              C3j_2 = 0.0d0

              call neufj(sl,sk,1.0d0,s1,sI1,1.0d0,s3,sI3,1.0d0,C9j_1)
              call neufj(sl,sk,1.0d0,sI2,sI1,sLf,sI4,sI3,sLi,C9j_2)
              call troisj(s1,s3,sl,0.0d0,0.0d0,0.0d0,C3j_1)
              call troisj(sI2,sI4,sl,0.0d0,0.0d0,0.0d0,C3j_2)

c       comodo definire una nuova funzione hat(j)=dsqrt(2.0d0*j+1.0d0) !!

              C  = C + ipot(l+I1+I2)*(2.0d0*l+1.0d0)*(2.0d0*k+1.0d0)*
     &             C9j_1* C9j_2* C3j_1* C3j_2* 
     &             dsqrt(2.0d0*(I1+i)+1.0d0)*dsqrt(2.0d0*I2+1.0d0)*
     &             dsqrt(2.0d0*(I3+j)+1.0d0)*dsqrt(2.0d0*I4+1.0d0)


            enddo
        enddo


11      lmin=max(iabs((I1+i)-(I3+j)),iabs(I2-I4))
        lmax=min((I1+i)+(I3+j),I2+I4)
                    
c       if(lmin.gt.lmax) goto 24

        do l=lmin,lmax   !``lambda''
           sl=dfloat(l)

           C3j_1 = 0.0d0
           C3j_2 = 0.0d0
           call troisj(s1,s3,sl,0.0d0,0.0d0,0.0d0,C3j_1)
           call troisj(sI2,sI4,sl,0.0d0,0.0d0,0.0d0,C3j_2)

           kmin=max(iabs(I1-I3),iabs(l-1))
           kmax=min(I1+I3,l+1)
                       
c          if(kmin.gt.kmax) goto 23
 
           do k=kmin,kmax   !``s''
              sk=dfloat(k)

              C9j_1 = 0.0d0
              C9j_2 = 0.0d0             
              call neufj(sl,sk,1.0d0,s1,sI1,1.0d0,s3,sI3,1.0d0,C9j_1)
              call neufj(sl,sk,1.0d0,sI2,sI1,sLf,sI4,sI3,sLi,C9j_2)
             
c       comodo definire una nuova funzione hat(j)=dsqrt(2.0d0*j+1.0d0) !!


              C  = C + ipot(l+I1+I2)*(2.0d0*l+1.0d0)*(2.0d0*k+1.0d0)*
     &             C9j_1* C9j_2* C3j_1* C3j_2* 
     &             dsqrt(2.0d0*(I1+i)+1.0d0)*dsqrt(2.0d0*I2+1.0d0)*
     &             dsqrt(2.0d0*(I3+j)+1.0d0)*dsqrt(2.0d0*I4+1.0d0)           

            
           enddo
        enddo


24      RETURN
        END


C..........................................................................


        SUBROUTINE subIntRad(i,j,I1,I2,I3,I4,l1,l3,Su)
        
        implicit real*8 (a-h,o-z)
        include 'param.qrpa'                 !containing nsp,ncf,nnn,(nmx?)
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
c        common/rad/nmax,dstep
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del   

c       write(1002,*)'entering subIntRad'


        dstep=del

        do ir=1,nmax

ct         wf(ir,I1)=(ir*dstep)**2                 ! testing radial part
ct         wf(ir,I2)=(ir*dstep)**2                 !
ct         wf(ir,I3)=(ir*dstep)**2                 !
ct         wf(ir,I4)=(ir*dstep)**2                 !

           R_I1 = 0.0d0
           R_I3 = 0.0d0
           radp=dstep*ir
           call subOpRad(ir,i,I1,l1,R_I1)
           call subOpRad(ir,j,I3,l3,R_I3)


           Su = Su + R_I1*R_I3*wf(ir,I2)*wf(ir,I4)*dstep
                 
        enddo

ct         Su_TEST = Su/(R_I1*R_I3)                ! testing radial part
ct         write(1002,*)'S_TEST',Su_TEST           !
     

        RETURN
        END


C..........................................................................


        SUBROUTINE subOpRad(ir,k,I,l,R)      !meglio definirla come funzione?

        implicit real*8 (a-h,o-z)
        include 'param.qrpa'                 !containing nsp,ncf,nnn,(nmx?)
        common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
c       common/rad/nmax,dstep  
        common /hf/nmax,nocc,nunocc,norb
        common /hf1/del
    
        dstep=del     

        radp=dstep*ir

ct      dwf(ir,I)=2*radp                        ! testing radial part
ct      wf(ir,I)=radp**2                        !


        R = (dwf(ir,I)/radp-1/(radp**2)*wf(ir,I)+
     &       ipot((k+1)/2)*(l+(1-k)/2)*wf(ir,I)/(radp**2))*
     &       dsqrt(l+(k+1.0d0)/2.0d0)

       
        RETURN
        END

