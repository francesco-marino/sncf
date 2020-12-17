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
      open(unit=39,status='unknown',name='reduced_EM_sp_old.dat')
      write(39,*) iorb,isflip,ispin
      open(unit=40,status='unknown',name='reduced_EM_phon.dat')
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
