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

      open(unit=390,status='unknown',name='reduced_EM_sp.dat')
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
                rdint=rdint+wf(ipoint,ix1)*wf(ipoint,iy1)
     &                *rdist**iorb*del                        !!! non ci va ispin??
                enddo

               cccc = rdint*yl(lx2,jx2,ly2,jy2,iorb)*Isgn     !!! non ci va ispin??
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

