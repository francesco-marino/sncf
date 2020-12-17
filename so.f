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

