
CCC   ****************************************************************
CCC   DISCRETE RPA IN H.O. BASIS WITH SKYRME INTERACTION             *
CCC   PLUS COUPLING WITH CONTINUUM STATES AND DOORWAY STATES         *
CCC   by Nguyen Van Giai and Gianluca Colo`                          *
CCC   **************************************************************** 

      subroutine gr0_j 
      implicit real*8 (a-h,o-z)
      parameter (nnn=180,nmx=12,nsp=300,ncf=680,ncf2=2*ncf,nft=3)
      parameter (npa=270,nho=38,ncht=7,nbins_j=223)
ccc   nnn=number of points of the radial grid for vectors read from HF code
ccc   nmx=greater or equal than nosc
ccc   nsp,ncf=greater or equal than the number of single particle states
ccc           and the number of particle-hole configurations, respectively
ccc   nft=maximum number of phonons included in doorway states
ccc   npa,nho=greater or equal to the numbers of particle and hole states
      character*5 skyrme
      character*80 file_hf,file_chanl,file_doo,file_sig1, 
     1file_sig2,file_td,file_indices,file_phonons,file_iso 
      character*80 file_phonons_td(nft) 
ccc   Checking standards: does this real*4 disturb ? CHECK ALL REAL*4 
      real*4 secv,tempo
      complex*16  slf,slfminus,slfpr1,slfpr2,slfminusa,slfminusb
      complex*16  wdself(ncf,ncf),wdcross(ncf,ncf)
      complex*16  wself,wcross,slfa,slfb,nvs1,nvs1minus
      common/bafa0/afa0
      common/biso/it1
      common/bisor/tt(ncf,ncf),tt1(ncf,ncf),tt2(ncf,ncf),
     1tt3(ncf,ncf),tt4(ncf,ncf)
      common/bptl/ipartl,nchanl(2),ichanl(ncht),nchanlt,
     1nos1,npo1,nossp,nossn,nholes
      common/bptlr/qtrnsf,spect(ncht),erin(ncht)
      common/blecp1/t0,t1,t2,t3,t13,alfe,x0,x1,x2,x3,anucl,znucl
      common/bwu1/nconf,ntot,nconf1,nconf2,nconf3,nconft
      common/bwu1r/sqbel(ncf),charge(ncf)
      common/bwu2/ih(ncf),ip(ncf),ecf(ncf),nig(ncf,2),nqqx(50,2),
     1nqqy(50,2),ngr(2),irpab,jtot,itot,ipi,ichex1,irad1
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      common/bpas/h,cof
      common/bnri/nri
      common/bnrir/dr
      common/pot/u(nnn),rmass,eta,xka,xmsm(nnn),sig,
     1sgne,dreg(nnn),dirg(nnn),c2d,s2d,sigma
      common/bdiam/aa(ncf,ncf),bb(ncf,ncf),evr(ncf),vecr(ncf2,ncf)
      common/bvari/iprint,nctr1,nctr2,nctr3,nctr4
      common/bvari1/icf1,icf2
      common/bspread/ionestate,inoint,ein1,ein2,ein1t,ein2t
      common/new/itd,ireal,idw
      common/bvin/aai(ncf,ncf),bbi(ncf,ncf)
      common/bvi2/aaiminus(ncf,ncf),bbiminus(ncf,ncf)
      common/bsum/sm1(2),s0(2),s1(2),s3(2)
      common/bsh/ishort2,idoo
      common/betalf/calf(7),cbet(7)
      common/bfile1/file_hf
      common/bfile2/file_doo
      common/bfile3/file_phonons_td
      common/bfile4/file_td
      common/bhosc/ahosc(2,0:12,1:28,nmx,nmx)
      common/bwg/wg(nnn,nmx),ewg(nmx),xor(nmx),dwg(nnn,nmx)   
      common/bfg/f0(nnn,2),g0(nnn,2) 
      common/bpot/xmn(nnn),xmp(nnn),vn(nnn),vp(nnn),yn(nnn),yp(nnn),
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)   
      common/blecp2/dt(nnn),taut(nnn),dt1(nnn),dt2(nnn),dp(nnn),dn(nnn) 
      common/bjost/nr
      common/bjostr/refp,aifp,refd,aifd
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/brintw/riw,riwminus,fw(nnn),fwminus(nnn)
      common/bvari12/nctrp1,nctrn1,nctr01,nctrp2,nctrn2,nctr02
      common/bvari34/nctrp3,nctrn3,nctr03,nctrp4,nctrn4,nctr04      
      common/bself/rslf12,aslf12,rslf34,aslf34
      common/bfoti/nfo(nft),lfo(nft),jfo(nft),isfo(nft),lifo(nft),
     1itrans(nft),ifon,ntotfo,iformat
      common/bfotr/omgfo(nft),alphasq(nft) 
      common/bdens/rintpp(npa,npa,nft,2),rinthh(nho,nho,nft,2),rintph
     1(npa,nho,nft,2)
      common/bove/ove(nsp,nsp)
      common/btt/tiso
      common/bisa_last/ein0 
      common/b_forceparam/skyrme
      common/bmax/nmaxt
      COMMON/CTR4/ECFMIN,ECFMAX,EPHOCUT,S0CUT,ESPMIN,ESPMAX
      COMMON/CTR2/IBASIS,ISPIN,IORB,IPAR,isflip,irad,ichex
      common/b_to_do/ichf,ihf,irpa
      common/hf/nmax,nocc,nunocc,norb
      common/b_occup/nos,npo,ioc01,ideg1,ioc02,ideg2,ioc03,ideg3,
     &ioc04,ideg4,ioc05,ideg5,ioc06,ideg6

      dimension nord(ncf2),vec(ncf2),tjl1(15),opr(nnn),rmat(ncf),   
     1bel(ncf2)
      dimension vint(ncf,ncf),expoint(nnn),drho(nnn),
     1rsqrho(nnn),rrho(nnn)
      dimension f1(2),g1(2),iikk(2),ngrg(2)
      dimension td(nnn)
      dimension imult_l(100),imult_j(100)
      dimension ss(nbins_j)
      dimension zxave(nbins_j,1,1),zyave(nbins_j,1,1)
      dimension enesup(nbins_j),eneinf(nbins_j)
      dimension anorma(nbins_j),ssstate(ncf)
      dimension focc(nsp)
      dimension ampli_l(0:1,0:10)
  100 format(16i3)  
  101 format(6e12.5)
 1010 format(6(1x,e12.5))
  103 format(2x,'Number of s.p. states i=',i3,' is larger than', 
     1' reserved dimension nsp=',i3)
  104 format(2x,'The norm is negative for state number',i3)
  105 format(2x,'Single particle states',//3x,'i',4x,'n',4x,'l',4x,'j', 
     14x,'lt',6x,'e',13x,'focc')
  106 format(1x,i3,3i5,'/2',i3,3x,e12.5,2x,e12.5)  
  107 format(//2x,'p-h configurations',//5x,'i',4x,'ih',4x,'ip',6x,   
     1'eph',8x,'b(el)') 
  108 format(2x,i4,2x,i4,2x,i4,3x,7e12.5)  
  109 format(/2x,'Moments   m-1=',e12.5,3x,'m0=',e12.5,3x,'m1=', 
     1e12.5,3x,'m3=',e12.5) 
  110 format(//2x,'state',5x,'energy',8x,'b(el)',6x,'fract(newsr)', 
     12x,'fract(ewsr)',2x,'charge',6x,'<T**2>') 
  111 format(3x,i3,3x,6e13.5,2x,6e13.5)   
  114 format(60('*'))
 1140 format(2x,'emin,emax= ',e12.5) 
  115 format(e12.5,2i3)
  116 format(2x,'cof,afa0,nri= ',f5.2,2x,f5.2,2x,i4)
  117 format(2x,'jtot,ltot,itot,ichex,irad,irpa= ',6i2)
  118 format(2x,'iread,ipartl= ',2i2)
  119 format(2x,'nos,npo,nossp,nossn,nosc,nholes= '
     1,5i3,i2) 
  120 format(2x,'itd,idw= ',3i2)
  122 format(2x,'iprint,icf1,icf2= ',i2,2i4)
  123 format(2x,'cutoff= ',e12.5)
  124 format(2x,'ifon,ntotfo,iformat= ',i2,i4,i2)
  125 format(2x,'omg,ga,ireal= ',f6.3,2x,f6.3,1x,i2)
  126 format(2x,'qtrnsf,nchanl(1),nchanl(2)= ',f7.4,2x,2i3)
  127 format(2x,'(ichanl(k3),k3=1,nchanlt)= ',10i3)
  128 format(2x,'ein0,nene,dene= ',f7.4,1x,i3,2x,f7.4)
  129 format(2x,'isw= ',i2)
  130 format(2x,'ionestate,inoint= ',2i3)
  131 format(2x,'ein1,ein2= ',e12.5,3x,e12.5)
 1310 format(2x,'ein1t,ein2t= ',e12.5,3x,e12.5)
  132 format(2x,'FOR THE NUCLEUS (A,Z) = ',f6.1,2x,f5.1)
  137 format(2x,'ishort1,ishort2= ',i2,i2/)
  247 format(1x,i4,1x,e12.5,1x,e12.5)
  248 format(10(1x,e12.5))
  717 format(2x,'ik= ',i3/)
  725 format(2e12.5,3i4)
 1070 format(5x,'Phonons coupled to single particle states'//,
     15x,'i',2x,'l',1x,'lt',2x,'e',11x,'perc',11x,'alphasq',
     24x,'itr'//)
 1071 format(2x,i4,i3,i3,i3,i3,1x,e12.5)
 1072 format(2x,i4,i3,i3,1x,e12.5,2x,e12.5)
 1200 format(2x,'nconf,nconft= ',2i4//)
 1201 format(2x,'nconf1,nconf2,nconft= ',3i4//)
 5562 format(2x,'afa0=',2e15.5,3x,'nosc=',2i4,3x,'irpa=',2i4)
 5563 format(2x,'ngrg(1),ngrg(2),nc =', 3i5)
 7792 format(5x,e16.8,5x,e16.8)
 7793 format(5x,e12.5,5x,e12.5,'   self-energy: ph vibrations')
 7794 format(5x,e12.5,5x,e12.5,'  induced int.: ph vibrations')
 7795 format(5x,e12.5,5x,e12.5,'   self-energy: pair vibrations')
 7796 format(5x,e12.5,5x,e12.5,'  induced int.: pair vibrations')

c     tempo=secnds(0.0)
ccc   Constants (units used are fm for length and MeV for energy) ****
      pi=4*datan(1.d0)
      sq4pi=dsqrt(4.d0*pi)
      hbc=197.3d0 
      dmshb=0.04823d0 
      bmas=0.5d0*dmshb*hbc*hbc

	write(2,*)
	write(2,*) 'Entering gr0_j'
	write(2,*)

c     open(unit=6,file='gr0.out',status='unknown')
      open(unit=7,file='pour_dw98.out',status='unknown')
ccc   Reading and expansion of s.p. states over h.o. basis ***********
      
c     read (5,*),emin,emax
      emin=ecfmin
	emax=ecfmax
c     read (5,*),cof,nri,afa0
      cof=1.d0
	nri=nmaxt-2
c     afa0 defined below
c     read (5,*),jtot,ltot,ispin,irad
      jtot=ispin
	ltot=iorb
	ispin=isflip
      ipi=1-2*mod(ltot,2) 
      if(jtot.ge.iabs(ltot-ispin).and.jtot.le.(ltot+ispin))go to 12
      write(2,*) '>>> J IS NOT CONSISTENT WITH L,S: '
      write(2,*) 'J, L, S = ',jtot,ltot,ispin 
      stop 
   12 irad1=irad
c     read (5,*),itot,ichex 
      itot=1
	ichex=1
	ichex1=ichex
c     read (5,*),irpa
      irpab=irpa
c     read (5,*),iread,ipartl
      iread=0
	ipartl=0
c     read (5,*),nos,npo,nossp,nossn,nosc,nholes 
      nosb=nos
	npob=npo
      nosc=12
	nholes=1
      if(nosc.gt.nmx)then
      write(2,*) ' PARAMETER NMX MUST BE INCREASED'
      stop
      end if
c     read (5,*),itd,idw
      itd=0
	idw=0
c     read (5,*),iprint,icf1,icf2
      iprint=0
	icf1=0
	icf2=0
c     read (5,*),ishort1,ishort2
      ishort1=1
	ishort2=1
      icomplete = 0 

      nr=nri/2
      if(cof.lt.0.001d0) call wsax
      if(cof.gt.0.001d0) call lecpot
C     hbom = 41 * (anucl**(-0.33333333d0)) 
      hbom=41.d0 * anucl**(-1.d0/3.d0)
C     afa0=HBOM*MC2/(HBC**2)
	afa0=hbom*939.d0/197.d0/197.d0
C_TEST
C     afa0=0.17
C_END_TEST
      dr=h
      if(cof.lt.0.001d0) rmass=(anucl+1.d0)/anucl   
      if(cof.gt.0.001d0) rmass=(anucl-1.d0)/anucl   
      tiso=(anucl-2.d0*znucl)/2.d0
      if(iread.eq.0) ipartl=0
      iprtl=ipartl+1
      if(dabs(afa0).gt.1e-8)bosc=1.d0/afa0  
      if(dabs(afa0).gt.1e-8)bosc=dsqrt(bosc)   
      lxmy=0
      write(7,*) '   '
      write(7,*) '*************** ---> INSERT IN ILECT = 1 ************'
      write(7,*) '*************** ---> AFTER THE LOGICAL FFFF *********'
      write(7,*) '   '
      icount_mul=0
      open(unit=99,status='old',file='lines')
      rewind(99)
      do 1 i=1,nos
      if(i.eq.1.and.iprint.lt.0)write(9,*) 'afa0=',afa0 
      if(dabs(afa0).gt.1e-8)then  
       read (99,*),nn(i),ll(i),lj(i),lt(i)
       if(iprint.lt.0)write(9,*) nn(i),ll(i),lj(i),lt(i)
       if(i.eq.1)then
        write(2,*) '  '
        write(2,*) ' Calculation with a HO basis: afa0 = ',afa0
	  write(2,*)
       end if
       if(i.eq.1)go to 1203
       do ii=1,icount_mul
        if(ll(i).eq.imult_l(ii).and.lj(i).eq.imult_j(ii))go to 1202
       end do
 1203  icount_mul=icount_mul+1
       imult_l(icount_mul)=ll(i)
       imult_j(icount_mul)=lj(i)
 1202  continue
       noc=nn(i)
       lp=ll(i)
       jp=lj(i)
       lt1=lt(i)+1  
c      write(2,*) 'bosc,nosc,lp,jp,lt1,noc = ',bosc,nosc,lp,jp,lt1,noc
C      if(i.eq.1) then
C       do jpr=1,nri
C       write(29,*) jpr,xmp(jpr)
C       end do
C      end if
       call qwg(bosc,nosc,lp,jp,lt1,noc) 
       ehf(i)=ewg(noc)
c      write(2,*) 'i,ehf(i) = ',i,ehf(i)
C_TEST
C       write(2,*) '  '
C       write(2,*) ' >>> TEST: USING GILLET-VINH MAU S.P. ENERGIES'
C       if(i.eq.1.or.i.eq.3)ehf(i)=-35.d0
C       if(i.eq.2.or.i.eq.4)ehf(i)=-18.72d0
C_END_TEST
       if(iprint.lt.0)write(9,*) noc,lp,jp,lt(i),ehf(i)
       do 300 j=1,nri
        dwf(j,i)=dwg(j,noc)
        wf(j,i)=wg(j,noc)
  300  continue
       if(iprint.lt.0)write(9,*) (wf(j,i),j=1,nri)
      else
       read (99,*),nb,lb,ljb,ltb
       if(i.eq.1)then
        write(2,*) '  '
        write(2,*) ' Calculation in a box: radius = ',float(nri)*dr
	  write(2,*)
       end if
       read(100,*) ehf(i),nn(i),ll(i),lj(i),lt(i)
       read(100,*)(wf(j,i),j=1,nri)
       read(100,*)(dwf(j,i),j=1,nri)
C_TEST
C       write(2,*) '  '
C       write(2,*) ' >>> TEST: USING GILLET-VINH MAU S.P. ENERGIES'
C       if(i.eq.1.or.i.eq.3)ehf(i)=-35.d0
C       if(i.eq.2.or.i.eq.4)ehf(i)=-18.72d0
C_END_TEST
       noc=nn(i)
       lp=ll(i)
       jp=lj(i)
       lt1=lt(i)+1
      end if
c     write(7,301) nn(i),ll(i),(lj(i)-2*ll(i)),lt(i)
c     write(7,302) anucl**(-1.d0/6.d0)
      xnor(i)=xor(noc) 
      lxmy=max0(lxmy,ll(i)) 
    1 continue  
      nos1=nos+1
      npo1=npo+1
      i=nos 
  200 read (99,*,end=201),lp,jp,ltp,n1,n2
c     if(lp.lt.0) go to 201
      ltp1=ltp+1
      lxmy=max0(lxmy,lp)
      ibb=(jp-2*lp+3)/2
      nqqy(lp+1,ibb)=n1
      nqqx(lp+1,ibb)=n2
      do 202 ii= n1,n2
      i=i+1
      if(i.le.nsp) go to 203
      write(2,103) i,nsp
  203 nn(i)=ii
      ll(i)=lp
      lj(i)=jp
      lt(i)=ltp
      if(dabs(afa0).gt.1e-8)then
       do ii3=1,icount_mul
        if(ll(i).eq.imult_l(ii3).and.lj(i).eq.imult_j(ii3))go to 2030
       end do
       icount_mul=icount_mul+1
       imult_l(icount_mul)=ll(i)
       imult_j(icount_mul)=lj(i)
 2030  continue
       call qwg(bosc,nosc,lp,jp,ltp1,ii)
       ehf(i)=ewg(ii)
       do 204 j=1,nri
       dwf(j,i)=dwg(j,ii)
       wf(j,i)=wg(j,ii)
  204  continue
C_TEST
C       if(i.eq.6.or.i.eq.10)ehf(i)=-4.95d0
C       if(i.eq.5.or.i.eq.9)ehf(i)=-1.86d0
C       if(i.eq.7.or.i.eq.11)ehf(i)=-1.1d0
C       if(i.eq.8.or.i.eq.12)ehf(i)=3.39d0
C_END_TEST
      else
       read(101,*) ehf(i)
       read(101,*)(wf(j,i),j=1,nri)
       read(101,*)(dwf(j,i),j=1,nri)
C_TEST
C       if(i.eq.6.or.i.eq.10)ehf(i)=-4.95d0
C       if(i.eq.5.or.i.eq.9)ehf(i)=-1.86d0
C       if(i.eq.7.or.i.eq.11)ehf(i)=-1.1d0
C       if(i.eq.8.or.i.eq.12)ehf(i)=3.39d0
C_END_TEST
      end if
      if(iprint.lt.0)write(9,*) (wf(j,i),j=1,nri+2)
c     write(7,301) nn(i),ll(i),(lj(i)-2*ll(i)),lt(i)
c     write(7,302) anucl**(-1.d0/6.d0)
  202 continue
      go to 200
  201 ntot=i
      do it=1,2
       do ii=1,icount_mul
        do jj=1,nosc
        icounts=jj+(ii-1)*nosc+(it-1)*nosc*icount_mul
        write(7,301) jj,imult_l(ii),(imult_j(ii)-2*imult_l(ii)), 
     &  icounts,it-1
        write(7,302) anucl**(-1.d0/6.d0)
        end do
       end do
      end do     
  301 format(3i5,7x,'/',1x,i4,i3)
  302 format(f10.5)
      lxmy=lxmy+1
      write(7,*) '   '
      if(iprint.lt.0)stop
      if(nos.gt.nho)then
      write(2,*) ' PARAMETER NHO MUST BE INCREASED'
      stop
      end if
      npaeff=ntot-nos
      if(npaeff.gt.npa)then
      write(2,*) ' PARAMETER NPA MUST BE INCREASED'
      stop
      end if
	close(99)

c     read (5,*),cutoff
      cutoff=0.00001
c     read (5,*),reno 
      reno=1.d0
c     read (5,*),ifon,ntotfo,iformat 
      ifon=0
	ntotfo=1
	iformat=0
      if(ntotfo.gt.nft)then
      write(2,*) ' PARAMETER NFT MUST BE INCREASED'
      stop
      end if

c     do ifile=1,ntotfo  
c     read (5,'(a)'),file_phonons_td(ifile) 
c     write(2,'(a)') file_phonons_td(ifile)
c     end do 

c     read (5,*),omg,ga,ireal 
      omg=1.
	ga=0.
	ireal=0
c     read (5,*),qtrnsf,nchanl(1),nchanl(2)
      qtrnsf=0.d0
	nchanl(1)=0
	nchanl(2)=0
      nchanlt=nchanl(1)+nchanl(2)
      if(nchanlt.gt.ncht)then
      write(2,*) ' PARAMETER NCHT MUST BE INCREASED'
      stop
      end if
c     read (5,*),ein0,nene,dene
c     read (5,*),isw
c     read (5,*),ionestate,inoint
c     read (5,*),ein1,ein2
c     read (5,*),ein1t,ein2t
c     read (5,*),i_out_dw98
c     read (5,*) ndelta
c     read (5,*),deltae,nbins,amplibin
c     read (5,*) ipartial
c     if(ipartial.eq.1)read (5,*) ioc,ideg

      do i=1,ntot
      if(i.le.nos) focc(i)=1.d0
      if(i.gt.nos) focc(i)=0.d0
      if(ipartial.eq.1.and.i.eq.ioc)focc(i)=float(ideg)/float(lj(i)+1)
      end do 

      write(2,114)
      if(irpa.eq.1)write(2,*) ' TDA CALCULATION WITH FORCE: ',
     1SKYRME
      if(irpa.eq.2)write(2,*) ' RPA CALCULATION WITH FORCE: ',
     1SKYRME
      if(irpa.eq.3)write(2,*) ' UNPERTURBED CALCULATION WITH FORCE: ',
     1SKYRME
      write(2,132) anucl,znucl
      write(2,114)
      write(2,*) ' >>>Make sure that NNN is the same as in the HF code'
      write(2,*) '          '
      write(2,1140) emin,emax 
      write(2,116) cof,afa0,nri
      write(2,117) jtot,ltot,itot,ichex,irad,irpa
      write(2,118) iread,ipartl
      write(2,119) nos,npo,nossp,nossn,nosc,nholes
      write(2,120) itd,idw
      write(2,122) iprint,icf1,icf2
      write(2,123) cutoff
      write(2,124) ifon,ntotfo,iformat
      write(2,125) omg,ga,ireal
      write(2,126) qtrnsf,nchanl(1),nchanl(2)
      write(2,128) ein0,nene,dene
      write(2,129) isw
      write(2,130) ionestate,inoint
      write(2,131) ein1,ein2
      write(2,1310) ein1t,ein2t
      write(2,137) ishort1,ishort2

      open(unit=13,status='old',file='chanl_imr_0.dat')
      do i=1,nchanlt
      read (13,*),ichanl(i),spect(i),erin(i)
      write(2,*)  'ichanl,spect,erin=',ichanl(i),spect(i),erin(i)
      end do
      close(13)

c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> Up to here ',secv,' sec. are needed.'
c     write(2,*) '          '

ccc   ---Explanations of read-ins: ------------------------------------
ccc   cof=0. for Woods-Saxon,1. for HF potential                      !
ccc   afa0=1./(b**2)=hbom*mc2/(hbc**2) and hbom=41*A**(-0.333333)     !
ccc   cutoff=cut-off on strength of real RPA basis vectors            !
ccc   qtrnsf=transferred momentum, in fm^-1                           !
ccc   nri=number of points used in radial integrals                   !
ccc   jtot,itot,ipi=spin,isospin,parity of RPA states                 !
ccc   ichex=1 selects charge-exchange states                          !
ccc   irad (if =1) gives radial dep. of operat.                       ! 
ccc   ispin (if =1) gives spin-flip                                   !
ccc   noss*=first hole considered (nossp for protons, nossn for n.)   !
ccc   n1,n2=number of nodes of first and last unocc. state(lp,jp,ltp) !
ccc   nos=number of occupied s.p. states                              !
ccc   irpa=1 if TDA,=2 if RPA,=3 if HF                                !
ccc   nosc=number of harmonic oscillator shells                       !
ccc   nn,ll,lj,lt,ehf=quantum numbers n,l,2*j,charge(0 or 1),         !
ccc                   s.p. energy                                     !
ccc   wf=corresponding w.f. normalized to \integral wf*wf*dr=1        !
ccc   iread=0 :calculates and punches results of discrete RPA         !
ccc   iread=1 :reads directly results of discrete RPA(main states)    !
ccc   ipartl=0 : strength function calculation                        !
ccc   ipartl=1 : the same + decay branching ratios calculation        !
ccc   nchanl(i)=number of decay channels                              !
ccc   isw=1: calculates the coupling of RPA states with p-h-phon. cnf.!
ccc   isw=2: calculates the coupling with the continuum               !
ccc   isw=3: calculates both couplings                                !
ccc   ireal=1: the calculation is complete                            !
ccc   ireal=2: the real part of W_down is skipped                     !
ccc   ifon=0: no spreading due to phonon couling is calculated        !
ccc   ifon=1: schematic transition densities are produced             !
ccc   ifon=2: t.d. stored in files are used                           !
ccc   iformat=1: t.d. are in files with r and td(r) on each line      !
ccc   iformat=2: t.d. are in files with (td(ir),ir=1,nri)             !
ccc   nholes=0: in selfph we use all holes                            !
ccc   nholes=1:              the holes used for ph configurations     !
ccc   iprint<0: calculation is incomplete                             !
ccc   iprint=0: details are not printed, calc. is complete            !
ccc   iprint=1: while calculating SIGMA, only one matrix element is   !
ccc             considered, and everything is printed out             !
ccc   iprint=2: same as before, but when building matrices in SPREAD  !
ccc   ishort1,ishort2: to speed up the calculation (iread=0 or 1)     !
ccc   ishort1=1: isospin calculation is neglected                     !
ccc   ishort2=1: check of bi-orthogonality and isospin are neglected  !
ccc   emin,emax: cut-offs in the energy on the p-h configurations     !
ccc   idw=1: data for the histogram of doorway st. are printed (un.14)!
ccc   -----------------------------------------------------------------

c     tempo=secnds(0.0)
ccc   Building parameters of residual interaction ********************
ccc   f0 is a matrix: f0(j,1) contains the isoscalar part and
ccc   f0(j,2) the isovector; the same for g0.
ccc   [ This was used to check with AUERBACH ]
ccx   t_au = -934.d0
ccx   x_au = 0.5d0 
ccx   f1(1) = 0.75d0*t_au
ccx   g1(1) = (0.5d0*x_au-0.25d0)*t_au 
ccx   do 60 j=1,nri
ccx   f0(j,1) = f1(1)
ccx   g0(j,1) = g1(1)
ccx60 continue
      f1(1)=0.75d0*t0   
      g1(1)=-0.25d0*t0*(1.d0-2.d0*x0)  
      do 6 j=1,nri   
      f0(j,1)=f1(1)+(t3/48.d0)*(dt(j)**alfe)*(3.d0*(alfe+1.d0)*(alfe+2. 
     1d0)+alfe*(1.d0-alfe)*(1.d0+2.d0*x3)*(((dn(j)-dp(j))/dt(j))**2))  
      g0(j,1)=g1(1)-(t3/24.d0)*(1.d0-2.d0*x3)*(dt(j)**alfe)
    6 continue
ccx   f1(2) = (-0.5d0*x_au-0.25d0)*t_au 
ccx   g1(2) = -0.25d0*t_au
ccx   do 61 j=1,nri
ccx   f0(j,2) = f1(2)
ccx   g0(j,2) = g1(2)
ccx61 continue
      f1(2)=-0.25d0*t0*(1.d0+2.d0*x0)  
      g1(2)=-0.25d0*t0
      do 8 j=1,nri   
      f0(j,2)=f1(2)-(t3/24.d0)*(1.d0+2.d0*x3)*(dt(j)**alfe)
      g0(j,2)=g1(2)-(t3/24.d0)*(dt(j)**alfe)
    8 continue  
      it1=itot+1
      go to (4,5), it1
    4 vv1=3.d0*t1/8
      vv2=t2*(5.d0+4.d0*x2)/4.d0
      go to 7
    5 vv1=-t1*(1.d0+2.d0*x1)/8.d0
      vv2=t2*(1.d0+2.d0*x2)/4.d0
    7 continue
   70 if (itot.eq.1) go to 72
      avv1=-3.d0*t1/32.d0
      avv2=-2.d0*avv1
      avv4=-0.25d0*t2*(x2+5.d0/4.d0)
      bvv1=-0.125d0*t1*(x1/2.d0-0.25d0)
      bvv2=-2.d0*bvv1
      bvv4=-0.25d0*t2*(x2/2.d0+0.25d0)
      go to 73
   72 avv1=(-t1/8.d0)*(-x1/2.d0-0.25d0)
      avv2=-2.d0*avv1
      avv4=-0.25d0*t2*(x2/2.d0+0.25d0)
      bvv1=-0.125d0*t1*(-0.25d0) 
      bvv2=-2.d0*bvv1
      bvv4=-0.25d0*t2*0.25d0
   73 calf(1)=avv1
      calf(2)=avv2
      calf(3)=avv2
      cbet(1)=bvv1
      cbet(2)=bvv2
      cbet(3)=bvv2
      do 500 iy=4,7
      iyy=1-2*(iy/6)
      calf(iy)=iyy*avv4
      cbet(iy)=iyy*bvv4
c     write(*,*) 'avv4,iyy,calf=',avv4,iyy,calf(iy)
  500 continue  
      l1g=iabs(jtot-1)  
      l2g=jtot+1
      l1=ltot    
      if(ltot.eq.0) l1=2
      funct=1.d0
      if(ltot.eq.0) funct=sq4pi
      do 234 j=1,nri 
      xj=dr*j
      go to (501,502),iprtl
  501 opr(j)=sq4pi*float(1-irad) + float(irad)*funct*xj**l1
      go to 234
  502 xj=xj*qtrnsf
      opr(j)=fa(jtot,xj)/xj   
  234 continue
ccc   Building configurations   **************************************   
      i=0   
      ik=0
  750 ik=ik+1
      do 205 ihh=1,ntot
      if (focc(ihh).le.0.01) go to 205
c     if (ihh.lt.nossp) go to 205
c     if (ihh.ge.npo1.and.ihh.lt.nossn) go to 205
      if (ichex.eq.1.and.ik.eq.1.and.lt(ihh).eq.1) go to 205
      if (ichex.eq.1.and.ik.eq.2.and.lt(ihh).eq.0) go to 205
      l1=ll(ihh)
      j1=lj(ihh)
      lt1=lt(ihh)   
      do 206 ipp=1,ntot
      if (focc(ipp).gt.0.99) go to 206   
      if (ipp.eq.ihh) go to 206
      etest=ehf(ipp)-ehf(ihh)  
      if(etest.le.emin.or.etest.ge.emax) go to 206 
      l2=ll(ipp)
      j2=lj(ipp)
      lt2=lt(ipp)   
      if(iabs(lt1-lt2).ne.ichex) go to 206  
      l12=l1+l2 
      l12=1-2*mod(l12,2)
      if(l12.ne.ipi) go to 206 
      j12m=iabs(j1-j2)/2
      j12p=(j1+j2)/2
      if((jtot-j12m)*(jtot-j12p)) 207,207,206
  207 i=i+1 
      ih(i)=ihh 
      ip(i)=ipp 
      ecf(i)=ehf(ipp)-ehf(ihh)  
  206 continue  
  205 continue 
      if(ichex.eq.0)go to 704 
      go to(701,702),ik
  701 nconf1=i
      nconf3=nconf1+1
      go to 703
  702 nconft=i
      nconf2=nconft-nconf1
  703 if(ik.le.1)go to 750
      go to 705
  704 nconf=i   
      nconft=2*nconf
  705 write(2,114)  
      write(2,105)  
      do 222 i=1,ntot
  222 write(2,106) i,nn(i),ll(i),lj(i),lt(i),ehf(i),focc(i) 
      istop=0
      if(ichex.eq.0.and.nconf.gt.ncf) istop=1
      nct = nconf1 + nconf2 
      if(ichex.eq.1.and.nct.gt.ncf) istop=1
      if(istop.eq.1)then 
      write(2,*) '   '
      write(2,*) ' PARAMETER NCF MUST BE INCREASED'
      stop
      end if 
      write(2,114)
      if(ichex.eq.0)write(2,1200) nconf,nconft
      if(ichex.eq.1)write(2,1201) nconf1,nconf2,nconft
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To build ph interaction and ph configurations '
c    1,secv,' sec. are needed.'
c     write(2,*) '          '
      if(iread.eq.1)go to 1000 
      if(ifon.ne.1.and.ishort1.eq.0)go to 9208
      if(ifon.ne.1.and.ishort1.eq.1)go to 9221
c     tempo=secnds(0.0)
ccc   Calculation of radial integrals for coupling with phonons   ****
ccc   we use different shapes for schematic transition densities:
ccc   d(rho)/dr (1), (r**2)*rho (2), r*rho (3).
ccc   cfr. Nucl.Phys. A327 (1979) 397; the alpha's are recalculated by 
ccc   program trans.
 9021 rho0=0.150225d0
      diffu=0.49d0
      rad=6.8d0   ! Ready for Pb208
      do 11 k=1,nri          
      expoint(k)=dexp((k*dr-rad)/diffu)
      drho(k)=-(rho0/diffu)*(expoint(k)/((1+expoint(k))**2))
      rr=k*dr
      rsqrho(k)=(rho0*((k*dr)**2))/(1+expoint(k))
      rrho(k)=(rho0*(k*dr))/(1+expoint(k))
   11 continue
ccc   the calculation starts here  
      do 9201 ip1=nos1,ntot
      do 9202 ip2=nos1,ntot
      ripp10=0.
      ripp11=0.
      ripp20=0.
      ripp31=0.
      do 9203 k=1,nri
      cpp10=wf(k,ip1)*wf(k,ip2)*f0(k,1)*drho(k)*dr
      cpp11=wf(k,ip1)*wf(k,ip2)*f0(k,2)*drho(k)*dr
      cpp20=wf(k,ip1)*wf(k,ip2)*f0(k,1)*rsqrho(k)*dr
      cpp31=wf(k,ip1)*wf(k,ip2)*f0(k,2)*rrho(k)*dr
      ripp10=ripp10+cpp10
      ripp11=ripp11+cpp11
      ripp20=ripp20+cpp20
      ripp31=ripp31+cpp31
 9203 continue
      k1=ip1-nos
      k2=ip2-nos
      rintpp(k1,k2,1,1)=ripp10
      rintpp(k1,k2,1,2)=ripp11
      rintpp(k1,k2,2,1)=ripp20
      rintpp(k1,k2,3,2)=ripp31
 9202 continue
 9201 continue
      do 9204 ih1=1,nos
      do 9205 ih2=1,nos
      rihh10=0.
      rihh11=0.
      rihh20=0.
      rihh31=0.
      do 9206 k=1,nri
      chh10=wf(k,ih1)*wf(k,ih2)*f0(k,1)*drho(k)*dr
      chh11=wf(k,ih1)*wf(k,ih2)*f0(k,2)*drho(k)*dr
      chh20=wf(k,ih1)*wf(k,ih2)*f0(k,1)*rsqrho(k)*dr
      chh31=wf(k,ih1)*wf(k,ih2)*f0(k,2)*rrho(k)*dr
      rihh10=rihh10+chh10
      rihh11=rihh11+chh11
      rihh20=rihh20+chh20
      rihh31=rihh31+chh31
 9206 continue
      k1=ih1
      k2=ih2
      rinthh(k1,k2,1,1)=rihh10
      rinthh(k1,k2,1,2)=rihh11
      rinthh(k1,k2,2,1)=rihh20
      rinthh(k1,k2,3,2)=rihh31
 9205 continue
 9204 continue
      do 9230 ip1=nos1,ntot
      do 9231 ih1=1,nos
      riph10=0.
      riph11=0.
      riph20=0.
      riph31=0.
      do 9232 k=1,nri
      cph10=wf(k,ip1)*wf(k,ih1)*f0(k,1)*drho(k)*dr
      cph11=wf(k,ip1)*wf(k,ih1)*f0(k,2)*drho(k)*dr
      cph20=wf(k,ip1)*wf(k,ih1)*f0(k,1)*rsqrho(k)*dr
      cph31=wf(k,ip1)*wf(k,ih1)*f0(k,2)*rrho(k)*dr
      riph10=riph10+cph10
      riph11=riph11+cph11
      riph20=riph20+cph20
      riph31=riph31+cph31
 9232 continue
      kp=ip1-nos
      kh=ih1
      rintph(kp,kh,1,1)=riph10
      rintph(kp,kh,1,2)=riph11
      rintph(kp,kh,2,1)=riph20
      rintph(kp,kh,3,2)=riph31
 9231 continue
 9230 continue
      if(ishort1.eq.1)go to 9221 
ccc   Overlap integrals (used for the exp. value of T^2) *************
 9208 do 9219 i=1,ntot
      do 9220 j=i,ntot
      ove(i,j)=0.d0
      ove(j,i)=0.d0
 9220 continue
 9219 continue
      do 9216 i=1,ntot
      do 9217 j=i,ntot
      do 9218 k=1,nri
      ove(i,j)=ove(i,j)+wf(k,i)*wf(k,j)*dr
 9218 continue
      ove(j,i)=ove(i,j)
 9217 continue
 9216 continue      
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To calculate overlap integrals ',secv
c    1,' sec. are needed.'
c     write(2,*) '          '
ccc   Building matrices aa and bb   **********************************
 9221 continue
c     tempo=secnds(0.0)
      do 135 i1=1,ncf
      do 135 i2=1,ncf
      aa(i1,i2)=0.d0
      bb(i1,i2)=0.d0
  135 continue
      ncond=1
      ik=1
      if(ichex.eq.0)go to 706
      ik=0
  707 ik=ik+1
      go to(708,709), ik
  708 ncond=1
      nconf=nconf1
      go to 706
  709 ncond=nconf3
      nconf=nconft
  706 sm1(ik)=0.d0
      s0(ik)=0.d0 
      s1(ik)=0.d0 
      s3(ik)=0.d0
      nctr1=0
      nctr2=0
      nctr3=0
      nctr4=0
      nctr5=0
      nctr6=0
      nctr7=0
      nctr8=0
      nctrp1=0
      nctrn1=0
      nctr01=0
      nctrp2=0
      nctrn2=0
      nctr02=0
      nctrp3=0
      nctrn3=0
      nctr03=0
      nctrp4=0
      nctrn4=0
      nctr04=0
      idoo=0
      if (ik.gt.1) go to 133
      if (ifon.eq.0) go to 121
      write(2,1070)
      if (ifon.eq.2) go to 92081
2080  open(unit=31,status='old',file='phonons.dat')
      do 10 k=1,ntotfo
      read(31,*) nfo(k),lfo(k),lifo(k),omgfo(k),
     1alphasq(k),itrans(k)
      write(2,*) nfo(k),lfo(k),lifo(k),omgfo(k),
     1alphasq(k),itrans(k)
   10 continue
      close(31)
      write(2,*) '>>> THE REST OF THE CODE IS NO MORE CONSISTENT '
      write(2,*) '>>> WITH THE OPTION IFON=1 '
      stop
      go to 121
92081 continue
      open(unit=31,status='old',file='phonons.dat')
      do 92082 k=1,ntotfo
      read(31,*) nfo(k),lfo(k),jfo(k),isfo(k),lifo(k),omgfo(k)
      write(2,1071) nfo(k),lfo(k),jfo(k),isfo(k),lifo(k),omgfo(k) 
      alphasq(k)=1.d0
92082 continue
      close(31)
      call transition
  121 write(2,114)  
  133 write(2,107)  
      do 208 i1=ncond,nconf
      if(iprint.eq.1.and.i1.ne.icf1) go to 208
      ip1=ip(i1)
      ih1=ih(i1)
      lp1=ll(ip1)   
      lh1=ll(ih1)   
      jp1=lj(ip1)   
      jh1=lj(ih1)
      xp1=lp1*(lp1+1.d0)
      xh1=lh1*(lh1+1.d0)
      zlp1=lp1
      zlh1=lh1
      zjp1=0.5d0*jp1
      zjh1=0.5d0*jh1
      yp1=1.d0
      yh1=1.d0
      if(lp1.eq.0) yp1=0.d0
      if(lh1.eq.0) yh1=0.d0   
      yph1=1.d0
      k=(jtot-iabs(lp1-lh1))*(jtot-lp1-lh1)
      if(k.gt.0) yph1=0.d0
      yl1=yl(lp1,jp1,lh1,jh1,jtot) 
      tjl2=tjl(lp1,jp1,lh1,jh1,jtot,ltot)
      x=0.d0  
      geom=yl1
      if(ispin.eq.1)geom=tjl2 
      do 235 j=1,nri 
  235 x=x+opr(j)*wf(j,ip1)*wf(j,ih1)
      rmat(i1)=dr*x*geom
c     write(*,*) x*dr,geom
      x=rmat(i1)**2 
      sm1(ik)=sm1(ik)+x/ecf(i1) 
      s0(ik)=s0(ik)+x   
      s1(ik)=s1(ik)+x*ecf(i1)   
      s3(ik)=s3(ik)+x*(ecf(i1)**3)  
      str=0.d0     
ccc   this is the average of <v(induced)**2>
      do 209 i2=i1,nconf
      if(iprint.eq.1.and.i2.ne.icf2) go to 209
      ip2=ip(i2)
      ih2=ih(i2)
      jp2=lj(ip2)
      jh2=lj(ih2)
      aaa=0.d0
      vvv=0.d0
      if(irpa.eq.3) go to 5559
ccc   induce must be defined
ccc   induce=0: Skyrme p-h interaction
ccc   induce=1: imaginary part of induced p-h interaction. In this case
ccc   the result of phme subroutine must be multiplied by i (or -i ?) 
      call phme(ip1,ih2,ip2,ih1,jtot,nri,dr,aaar,aaarm,0)
      aaar=aaar
     1*dsqrt((focc(ih1)-focc(ip1)))*dsqrt((focc(ih2)-focc(ip2)))
      aaar=aaar*reno
ccc   isospin factor
      if (ichex.eq.1) then
      aaar= 2.d0*aaar
      aaarm=2.d0*aaarm
      vvvr= 2.d0*vvvr
      vvvrm=2.d0*vvvrm
      end if
c_TEST
c     filling for 58Ni
c     if(ip1.eq.65)aaar=aaar/dsqrt(2.d0)
c     if(ip2.eq.65)aaar=aaar/dsqrt(2.d0)
c_end_TEST
      vint(i1,i2)=aaar
      vint(i2,i1)=aaar
      aaa=aaar
      vvv=vvvr
      vvvm=vvvrm
 5559 continue
      if (i1.eq.i2) aaapr=aaa
      if (i1.eq.i2) aaa=aaa+ecf(i1)  
      slf=(0.d0,0.d0)
      slfminus=(0.d0,0.d0)
      slf1=(0.d0,0.d0)
      slf2=(0.d0,0.d0)
      slf3=(0.d0,0.d0)
      slf4=(0.d0,0.d0)
      slfa=(0.d0,0.d0)
      slfb=(0.d0,0.d0)
      slfminusa=(0.d0,0.d0)
      slfminusb=(0.d0,0.d0)
      nvs1=(0.d0,0.d0)
      nvs1minus=(0.d0,0.d0)
      if (ifon.eq.0) go to 7777
      if (lt(ip1).ne.lt(ip2)) go to 7777
ccx   if (i1.eq.i2) idt=1   ! Probably not used any more... VERIFY ! 
ccx   if (i1.ne.i2) idt=0
 5570 continue
      call selfph(ik,ip1,ih1,ip2,ih2,omg,ga,slfa,slfminusa)
      nvs1=slfa+slfb
      nvs1minus=slfminusa+slfminusb
c     write(*,*) 'OUTSIDE SELFPH; THIS IS SLF: ',slf
c     write(*,*) 'OUTSIDE SELFPH; THIS IS SLFMINUS: ',slfminus
 5572 continue
 7777 aa(i1,i2)=aaa
      aa(i2,i1)=aaa
 7785 if (i1.eq.i2) then 
       slfpr1=nvs1
       slfpr2=nvs1minus
      end if 
      aai(i1,i2)=dreal(nvs1)
      bbi(i1,i2)=dimag(nvs1)
      aaiminus(i1,i2)=dreal(nvs1minus)
      bbiminus(i1,i2)=dimag(nvs1minus)
      if(icomplete.eq.0.or.i1.eq.i2)then
      aai(i2,i1)=dreal(nvs1)
      bbi(i2,i1)=dimag(nvs1)
      aaiminus(i2,i1)=dreal(nvs1minus)
      bbiminus(i2,i1)=dimag(nvs1minus)
      end if
      if(icomplete.eq.1.and.i1.ne.i2)then
      call selfph(ik,ip2,ih2,ip1,ih1,omg,ga,slfa,slfminusa)
      nvs1=slfa+slfb
      nvs1minus=slfminusa+slfminusb
      aai(i2,i1)=dreal(nvs1)
      bbi(i2,i1)=dimag(nvs1)
      aaiminus(i2,i1)=dreal(nvs1minus)
      bbiminus(i2,i1)=dimag(nvs1minus)
      end if
      if (lt(ip1).ne.lt(ip2)) then
      rslf12=0.d0
      aslf12=0.d0
      rslf34=0.d0
      aslf34=0.d0
      go to 7784
      end if
 7784 wdself(i1,i2)=dcmplx(rslf12,aslf12)
      wdself(i2,i1)=dcmplx(rslf12,aslf12)
      wdcross(i1,i2)=dcmplx(rslf34,aslf34)
      wdcross(i2,i1)=dcmplx(rslf34,aslf34)
      if(iprint.eq.1)write(2,7793) wdself(i1,i2)
      if(iprint.eq.1)write(2,7794) wdcross(i1,i2)
      str=str+(aslf34**2)
 7787 if(ichex.eq.1)go to 210
      bbb=0.d0
      if(irpa.eq.1) go to 210
      if(irpa.eq.3) go to 210   
      jbb=(jp2+jh2)/2   
      sgn=1.d0-2.d0*mod(jbb,2)
ccc   induce must be defined
ccc   same comments as above (q.v.)
      call phme(ip1,ip2,ih2,ih1,jtot,nri,dr,bbb,bbbm,0)
      bbb=bbb
     1*dsqrt((focc(ih1)-focc(ip1)))*dsqrt((focc(ih2)-focc(ip2)))
      bbb=bbb*reno
      bbb=sgn*bbb 
      www=sgn*www
      wwwm=sgn*wwwm
      bb(i1,i2)=bbb 
      bb(i2,i1)=bbb 
      if(i1.eq.i2) bbbpr=bbb
  210 if (iprint.eq.1) go to 211
  209 continue
      if(ichex.eq.0)go to 211
      if(irpa.ne.2)go to 211
      if(ik.eq.2)go to 713
      do 711 i2=nconf3,nconft
      ip2=ip(i2)
      ih2=ih(i2)
      jp2=lj(ip2)
      jh2=lj(ih2)
      bbb=0.d0
      jbb=(jp2+jh2)/2   
      sgn=1.d0-2.d0*mod(jbb,2)
ccc   induce must be defined
ccc   same comments as above (q.v.)
      call phme(ip1,ip2,ih2,ih1,jtot,nri,dr,bbb,bbbm,0)
      bbb=bbb
     1*dsqrt((focc(ih1)-focc(ip1)))*dsqrt((focc(ih2)-focc(ip2)))
      bbb=bbb*reno
      bbb=sgn*bbb 
      www=sgn*www
      wwwm=sgn*wwwm
      if (ichex.eq.1) then
      bbb=2.d0*bbb
      www=2.d0*www
      wwwm=2.d0*wwwm
      end if
      bb(i1,i2)=bbb 
      if(i1.eq.i2) bbbpr=bbb
  212 if (iprint.eq.1) go to 211
  711 continue
  713 continue
      do 714 i2=1,nconf1
      bb(i1,i2)=bb(i2,i1)
  714 continue
  211 write(2,108) i1,ih1,ip1,ecf(i1),x,aaapr,slfpr1,slfpr2 
      if (iprint.eq.1) stop
      str=str/(float(nconf-i1+1))
  208 continue 
      write(2,109) sm1(ik),s0(ik),s1(ik),s3(ik) 
      write(2,114)
      do i1=1,3
       do i2=1,3
       e=0.d0
       if(i1.eq.i2)e=ecf(i1)
       aaapr=aa(i1,i2)-e
       write(2,2081) i1,i2,e,aaapr,bb(i1,i2)
 2081  format(1x,2i3,1x,3(e12.5,1x))
       end do
      end do
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To build matrices ',secv,' sec. are needed.'
c     write(2,*) '          '
      if(ifon.eq.0.or.iprint.eq.1)go to 2084 
ccc   The effect of interference among diagrams **********************
ccc   HERE    changed from 2080 to 3080  
 2080 continue
      wself=(0.d0,0.d0)
      wcross=(0.d0,0.d0)
      vav=0.d0
      denom=0.d0
      do 2090 k1=ncond,nconf
      do 2091 k2=ncond,nconf
      wself=wself+rmat(k1)*wdself(k1,k2)*rmat(k2)
      wcross=wcross+rmat(k1)*wdcross(k1,k2)*rmat(k2)
      vav=vav+rmat(k1)*vint(k1,k2)*rmat(k2)
 2091 continue
      denom=denom+rmat(k1)*rmat(k1)
 2090 continue
      wself=wself/denom
      wcross=wcross/denom
      vav=vav/denom
      aself=2*dimag(wself)
      across=2*dimag(wcross)
      write(2,114)
      write(2,*) 'Mean values with matrix elements of the excitation '
      write(2,*) 'operator as weighting factors: '
      write(2,*) 'Mean value of A(i,j) =   ',vav
      write(2,*) 'Mean value of Im Sigma(i,j) [self-energies] =  ',aself
      write(2,*) 'Mean value of Im Sigma(i,j) [cross diagrams] =   ',across
      write(2,114)
      write(2,*) 'nctr1,nctr2,nctr3,nctr4=',nctr1,nctr2,nctr3,nctr4
      write(2,*) 'nctr5,nctr6,nctr7,nctr8=',nctr5,nctr6,nctr7,nctr8
      write(2,*) 'nctrp1,nctrp2,nctrp3,nctrp4=',nctrp1,nctrp2,nctrp3,nctrp4
      write(2,*) 'nctrn1,nctrn2,nctrn3,nctrn4=',nctrn1,nctrn2,nctrn3,nctrn4
      write(2,*) 'nctr01,nctr02,nctr03,nctr04=',nctr01,nctr02,nctr03,nctr04
ccc   Diagonalization and normalization    ***************************
 2084 continue
      if(ichex.eq.1.and.ik.eq.1)go to 707
c     tempo=secnds(0.0)
      if(ichex.eq.1)nconf=nconft
      nc=nconf
      irpa1=1
      if(irpa.eq.2) irpa1=2
      if(irpa.eq.2) nc=2*nconf
      nc1=nc+1  
      write(2,*) 'Entering diagma'
      call diagma(nconf,irpa1)  
      write(2,*) 'Exiting diagma'
      do 218 i=1,nconf   
      if(evr(i).le.0.d0) go to 218
      xno=0.d0
      do 219 j=1,nc  
      sg=1.d0 
      if(j.gt.nconf) sg=-1.d0 
      xno=xno+sg*(vecr(j,i)**2) 
  219 continue
      if(xno.gt.0.d0) go to 220   
      write(2,104) i
      stop  
  220 xno=dsqrt(xno) 
      do 221 j=1,nc  
  221 vecr(j,i)=vecr(j,i)/xno   
  218 continue  
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To diagonalize and normalize ',secv
c    1,' sec. are needed.'
c     write(2,*) '          '
ccc   Ordering of states   *******************************************
c     tempo=secnds(0.0)
      do 224 i=1,nconf   
  224 nord(i)=i 
      do 225 i=1,nconf   
      k=i   
      x=evr(i)  
      do 226 j=1,nc  
  226 vec(j)=vecr(j,i)  
      do 227 ii=i,nconf  
      if(x-evr(ii)) 227,227,228 
  228 k=ii  
      x=evr(ii) 
      do 229 j=1,nc  
  229 vec(j)=vecr(j,ii) 
  227 continue  
      kk=nord(k)
      nord(k)=nord(i)   
      nord(i)=kk
      evr(k)=evr(i) 
      do 230 j=1,nc  
  230 vecr(j,k)=vecr(j,i)   
      evr(i)=x  
      do 231 j=1,nc  
  231 vecr(j,i)=vec(j)  
  225 continue  
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To order the states ',secv,' sec. are needed.'
c     write(2,*) '          '
      if(itot.eq.0.or.ishort1.eq.1) go to 2270
ccc   Matrix elements of T**2   **************************************
c     tempo=secnds(0.0)
      open(unit=50,status='new',file='isospin.dat') 
      tiso1=tiso*(tiso-1.d0)
      do 2252 nevr1=1,nconf
      do 2261 nevr2=nevr1,nconf
      tsum1=0.d0
      do 2250 i=1,npo
      do 2251 j=npo1,nos
      if (ll(i).ne.ll(j)) go to 2251
      if (lj(i).ne.lj(j)) go to 2251
      tsum1=tsum1 + (ove(i,j)**2)*(float(lj(i))+1.d0)
 2251 continue
 2250 continue
      xx=0.d0      
      do 2262 k=1,nconf
      xx=xx+vecr(k,nevr1)*vecr(k,nevr2)
 2262 continue
      tt1(nevr1,nevr2)=(znucl-tsum1)*xx         
ccc   term a1 of my notes (first copybook on IAR, p.9)
      tsum2=0.d0
      do 2253 n1=1,nconf
      ipp1=ip(n1)
      do 2254 n2=1,nconf
      ipp2=ip(n2)
      if (ih(n1).ne.ih(n2)) go to 2254
      if (ll(ipp1).ne.ll(ipp2)) go to 2254
      if (lj(ipp1).ne.lj(ipp2)) go to 2254
      do 2255 j=npo1,nos
      if (ll(ipp1).ne.ll(j)) go to 2255
      if (lj(ipp1).ne.lj(j)) go to 2255
      tc2=vecr(n1,nevr1)*vecr(n2,nevr2)*ove(ipp1,j)*ove(ipp2,j)
      tsum2=tsum2+tc2
 2255 continue
 2254 continue
 2253 continue
      tt2(nevr1,nevr2)=1.d0-tsum2          
ccc   term a2 of my notes (p.10)
      tsum3=0.d0
      do 2256 n1=1,nconf
      ihh1=ih(n1)
      do 2257 n2=1,nconf
      ihh2=ih(n2)
      if (ip(n1).ne.ip(n2)) go to 2257
      if (ll(ihh1).ne.ll(ihh2)) go to 2257
      if (lj(ihh1).ne.lj(ihh2)) go to 2257
      do 2258 i=1,npo
      if (ll(ihh1).ne.ll(i)) go to 2258
      if (lj(ihh1).ne.lj(i)) go to 2258
      tc3=vecr(n1,nevr1)*vecr(n2,nevr2)*ove(i,ihh1)*ove(i,ihh2)
      tsum3=tsum3+tc3
 2258 continue
 2257 continue
 2256 continue
      tt3(nevr1,nevr2)=tsum3          
ccc   term a3 of my notes (p.10)
      tsum4=0.d0
      do 2259 n1=1,nconf
      ipp1=ip(n1)
      ihh1=ih(n1)
      if (ll(ipp1).ne.ll(ihh1)) go to 2259
      if (lj(ipp1).ne.lj(ihh1)) go to 2259
      do 2260 n2=1,nconf
      ipp2=ip(n2)
      ihh2=ih(n2)
      if (ll(ipp2).ne.ll(ihh2)) go to 2260
      if (lj(ipp2).ne.lj(ihh2)) go to 2260
      tc4=vecr(n1,nevr1)*vecr(n2,nevr2)*ove(ipp1,ihh1)*ove(ipp2,
     1   ihh2)*dsqrt((float(lj(ipp1))+1.d0)*(float(lj(ipp2))+1.d0))
      tsum4=tsum4+tc4
 2260 continue
 2259 continue
      tt4(nevr1,nevr2)=tsum4          
ccc   term a4 of my notes (p.11)
      tt(nevr1,nevr2)=tt1(nevr1,nevr2)+tt2(nevr1,nevr2)+
     1                tt3(nevr1,nevr2)+tt4(nevr1,nevr2)
      if (nevr1.eq.nevr2) tt(nevr1,nevr2)=tt(nevr1,nevr2)+tiso1
      write(50,*) nevr1,nevr2,tt(nevr1,nevr2)
ccx   write(50,*) tt1(nevr1,nevr2),tt2(nevr1,nevr2),tt3(nevr1,nevr2),
ccx  1            tt4(nevr1,nevr2)
 2261 continue
 2252 continue
      close(50)
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To calculate expect. values of T**2 '
c    1,secv,' sec. are needed.'
c     write(2,*) '          '
ccc   Transition probabilities and charges   *************************
 2270 continue
c     tempo=secnds(0.0)
      open(unit=21,status='unknown',file='z+1.dat')
      open(unit=22,status='unknown',file='z-1.dat')
      do 715 ik=1,2
      sm1(ik)=0.d0
      s0(ik)=0.d0 
      s1(ik)=0.d0 
      s3(ik)=0.d0 
  715 continue
      if(jtot.gt.8)then
       stop '>>> INCREASE DIMENSION OF AMPLI_L'
      end if
      do 236 i=1,nconf   
      if(evr(i).le.0.d0) go to 236
      x=0.d0  
      if(ispin.eq.1.and.jtot.ne.0)then
       do il=jtot-1,jtot+1,1
       do is=0,1
       ampli_l(is,il)=0.d0
       end do
       end do
      end if
      charge(i)=0.d0
      do 237 ii=1,nconf  
      y=vecr(ii,i)*dsqrt(focc(ih(ii))-focc(ip(ii)))
      if(irpa.eq.2) then
       isignrel = 1
       if(ispin.eq.1.and.jtot.eq.ltot)isignrel = -1 
       y=y-isignrel*vecr(ii+nconf,i)*dsqrt(focc(ih(ii))-focc(ip(ii)))  
      end if 
      x=x+y*rmat(ii)
      sgn=1.d0
      if(ii.gt.nconf1)sgn=-1.d0
      charge(i)=charge(i)+sgn*(vecr(ii,i)**2+vecr(ii+nconf,i)**2)
      if(ispin.eq.0.or.jtot.eq.0) go to 237 
       do il=jtot-1,jtot+1,1
       do is=0,1
       xl=dfloat(il)
       xj=dfloat(jtot)
       xs=dfloat(is)
       xl1=dfloat(ll(ip(ii)))
       xl2=dfloat(ll(ih(ii)))
       xj1=dfloat(lj(ip(ii)))/2.d0
       xj2=dfloat(lj(ih(ii)))/2.d0
       call NEUFJ(xl,xs,xj,xl1,0.5d0,xj1,xl2,0.5d0,xj2,c9j) 
       ctr=dsqrt(dfloat(lj(ip(ii)))+1.d0)*
     1     dsqrt(dfloat(lj(ih(ii)))+1.d0)*
     1     dsqrt(dfloat(2*il+1))*dsqrt(dfloat(2*is+1))*
     1     c9j
       ctr=ctr**2
c      y=vecr(ii,i)*dsqrt(focc(ih(ii))-focc(ip(ii)))
       y=vecr(ii,i)
       if(irpa.eq.2) then
        isignrel = 1
c       if(ispin.eq.1.and.jtot.eq.ltot)isignrel = -1 
        y=y-isignrel*vecr(ii+nconf,i)  
       end if 
       y=y**2
       ampli_l(is,il)=ampli_l(is,il)+ctr*y
       end do
       end do
  237 continue
      sqbel(i)=x
      x=x*x 
      bel(i)=x  
      ik=1
      if(ichex.eq.1.and.charge(i).gt.0.d0)ik=1
      if(ichex.eq.1.and.charge(i).le.0.d0)ik=2
      sm1(ik)=sm1(ik)+x/evr(i)  
      s0(ik)=s0(ik)+x   
      s1(ik)=s1(ik)+x*evr(i)
      s3(ik)=s3(ik)+x*(evr(i)**3)   
      if(ispin.eq.0.or.jtot.eq.0) go to 236
      anorma2=0.d0
      do il=jtot-1,jtot+1,1
      do is=0,1
      anorma2=anorma2+ampli_l(is,il)
      end do
      end do
      anorma2=dsqrt(anorma2)
      write(70,2370) i,evr(i),ampli_l(0,jtot-1),ampli_l(0,jtot),
     1ampli_l(0,jtot+1)
      write(70,2370) i,evr(i),ampli_l(1,jtot-1),ampli_l(1,jtot),
     1ampli_l(1,jtot+1),anorma2
      ratio=ampli_l(1,jtot+1)/ampli_l(1,jtot-1)
      write(71,*) i,evr(i),charge(i),ratio
 2370 format(1x,i4,5(2x,e12.5))
  236 continue  
      ikfin=1
      if(ichex.eq.1)ikfin=2
      do 716 ik=1,ikfin
      write(2,114)  
      if(ichex.eq.1)write(2,717) ik
      write(2,109) sm1(ik),s0(ik),s1(ik),s3(ik) 
  716 continue
      if(ichex.eq.1)then
       write(2,114)
       write(2,*)
       write(2,*) 'M0(-1)-M0(+1)=',s0(1)-s0(2)
       write(2,*) 'M1(-1)+M1(+1)=',s1(1)+s1(2)
       write(2,*)
      end if
      write(2,114) 
      write(2,*) 'Analytic moments: '
      if(ichex.eq.0) then 
       write(2,*) '>>> This version does not print the m'
       go to 7160
      end if
      write(2,*) '...'  
 7160 write(2,114) 
      write(2,110)  
      do j=1,2
      iikk(j)=0  
      end do
      do 238 i=1,nconf   
      if(evr(i).lt.0.d0) go to 238
      ik=1
      if(ichex.eq.1.and.charge(i).le.0.d0)ik=2
      x=bel(i)/s0(ik)   
      iikk(ik)=iikk(ik)+1   
      indx=iikk(ik)
      nig(indx,ik)=i 
  239 continue  
      y=evr(i)*bel(i)/s1(ik)
      write(2,111) i,evr(i),bel(i),x,y,charge(i),tt(i,i)
      if(charge(i).gt.0.d0)write(21,2380) i,evr(i),bel(i)
      if(charge(i).lt.0.d0)write(22,2381) i,evr(i),bel(i),evr(i)+deltae
  238 continue 
 2380 format(1x,i4,1x,e12.5,1x,e12.5)
 2381 format(1x,i4,1x,e12.5,1x,e12.5,1x,e12.5)
      do j=1,2 
      ngr(j)=iikk(j)
      end do
  240 continue
      write(2,*) 'End of section in which results are written'
      if(itd.eq.0)go to 252
ccc   Transition densities   *****************************************
      open(unit=19,status='new',file='TD.dat')
      do 241 i=1,nconf
      if(evr(i).le.0.d0)go to 241
      do 242 j=1,nri
      td(j)=0.d0
  242 continue
  243 do 246 ii=1,nconf 
      ampli=vecr(ii,i)
      if(irpa.eq.2)ampli=ampli-vecr(ii+nconf,i)
      ipp=ip(ii)
      ihh=ih(ii)
      lp=ll(ipp)
      lh=ll(ihh)      
      jp=lj(ipp)
      jh=lj(ihh)      
      tdmat=yl(lp,jp,lh,jh,jtot)
      if(ispin.eq.1)tdmat=tjl(lp,jp,lh,jh,jtot,ltot)
      do 245 ir=1,nri
      xj=float(ir)*dr
      x2=xj**2
      td(ir)=td(ir)+(ampli*tdmat*wf(ir,ipp)*wf(ir,ihh))/x2 
  245 continue
  246 continue
      aint=0.d0
      do 250 ir=1,nri
      xj=float(ir)*dr
      x2=xj**2
      aint=aint+dr*x2*opr(ir)*td(ir) 
  250 continue
      areno=aint/sqbel(i)
      do 249 ir=1,nri
      td(ir)=td(ir)/areno
  249 continue
      aint=0.d0
      do 251 ir=1,nri
      xj=ir*dr
      x2=xj**2
      aint=aint+dr*x2*opr(ir)*td(ir) 
  251 continue
      write(19,247) i,charge(i),evr(i)
      write(19,*) areno,sqbel(i)
      do j=1,nri
      write(19,*) j*dr,td(j)
      end do
  241 continue
      close(19) 
ccc   Outputs   ******************************************************
ccc   Also the output for DW98 is included; see my notes.
ccc   IMPORTANT: in the calculation, 1 refer to the hole and 2 to the
ccc              particle. Since DW98 wants the indices in opposite
ccc              order, they are print out reversed.
  252 open(unit=11,status='unknown',file='rpa.dat')
      write(2,*) 'File rpa opened'
      if(ifon.ne.0)then 
       open(unit=15,status='new',file='sigma1.dat')
       open(unit=16,status='new',file='sigma2.dat')
      end if 
c     write(7,*) 'Writing the wf of the state n.',i_out_dw98
c     write(7,*) 'E = ',evr(i_out_dw98)
      write(7,*) ''
      write(7,*) '************* ---> INSERT IN ILECT = 6 **************'
      write(7,*) '************* ---> AFTER THE LINE AFTER FFFFF *******'
      dummy=0.d0

c_20_01_04: 1 bin=1 state !
c           ***************
      ng3=ngr(1)+ngr(2)
      nbins=0
c     write(2,*) 'New part',i_out_dw98
      do 2522 i=1,ng3
c     write(2,*) i,evr(i),charge(i)
      if(i_out_dw98.ne.1)go to 2522
      if(charge(i).gt.0.d0)go to 2522
c     (t,3He)=(n,p) is considered now, so we go to z-1
      if(evr(i).gt.40.d0)go to 2522
      nbins=nbins+1
c     write(2,*) nbins
      enesup(nbins)=evr(i)+deltae+0.001d0
      eneinf(nbins)=evr(i)+deltae-0.001d0
 2522 continue
      write(2,*) 'nbins= ',nbins
c_20_01_04_END

      if(nbins.gt.nbins_j)then
       write(2,*)
       write(2,*) '>>> Too many bins !'
       write(2,*) 
       stop
      end if
      write(70,*) nbins
      do istre=1,nbins
c_20_01_04
c      enesup(istre)=dfloat(istre)*amplibin
c      eneinf(istre)=enesup(istre)-amplibin
c_20_01_04_END
       write(70,*) istre,eneinf(istre),enesup(istre)
       ss(istre)=0.d0
       anorma(istre)=0.d0
       do jstre1=1,1
        do jstre2=1,1
         zxave(istre,jstre1,jstre2)=0.d0
         zyave(istre,jstre1,jstre2)=0.d0
        end do
       end do
      end do

      write(11,725) s0(1),s0(2),ngr(1),ngr(2),nc
      write(11,115) afa0,nosc,irpa
      ng3=ngr(1)+ngr(2)
      do 1001 i=1,ng3
      write(11,101) evr(i),sqbel(i),bel(i),charge(i)
      write(11,1010) (vecr(j,i),j=1,nc)
      if(i_out_dw98.ne.1)go to 1001
      if(charge(i).gt.0.d0)go to 1001
c     (t,3He)=(n,p) is considered now, so we go to z-1

      istre_true=100000
      do 10001 istre=1,nbins
c_TEST_48K
c     we put a state in the bin x,x+1 MeV if its energy in the daughter
c     lies in that interval. The energy in the daughter is evr(i)+deltae. 
c     deltae=-11.31 MeV for 48Ca-48K, deltae=+0.40 MeV for 58Ni-58Co.
       ene_test=evr(i)+deltae
       if(ene_test.gt.eneinf(istre).and.
     & ene_test.lt.enesup(istre))then
        istre_true=istre
        write(70,*) i,evr(i),evr(i)+deltae,istre_true
       end if
10001 end do
      if(istre_true.eq.10000)go to 1001

      anorma(istre_true)=anorma(istre_true)+bel(i)
      ssstate(i)=0.d0  

      do 2524 it1=1,2
      do 2525 ii1=1,icount_mul
      l1=imult_l(ii1)
      j1=imult_j(ii1)
      do 2526 jj1=1,nosc
      icounts1=jj1+(ii1-1)*nosc+(it1-1)*nosc*icount_mul

      do 2527 it2=1,2
      if(iabs(it1-it2).ne.ichex) go to 2527
      do 2528 ii2=1,icount_mul
      l2=imult_l(ii2)
      j2=imult_j(ii2)
      l12=l1+l2 
      l12=1-2*mod(l12,2)
      if(l12.ne.ipi) go to 2528
      j12m=iabs(j1-j2)/2
      j12p=(j1+j2)/2
      if((jtot-j12m)*(jtot-j12p).gt.0)go to 2528 
      do 2529 jj2=1,nosc
      icounts2=jj2+(ii2-1)*nosc+(it2-1)*nosc*icount_mul

      zx=0.d0
      zy=0.d0
      do 2523 j=1,nconf
      if((it1-1).ne.lt(ih(j)))go to 2523 
      if((it2-1).ne.lt(ip(j)))go to 2523 
      if(l1.ne.ll(ih(j)))go to 2523
      if(l2.ne.ll(ip(j)))go to 2523
      if(j1.ne.lj(ih(j)))go to 2523
      if(j2.ne.lj(ip(j)))go to 2523
      zx=zx+vecr(j,i)
     1*ahosc(it1,l1,j1,jj1,nn(ih(j)))*ahosc(it2,l2,j2,jj2,nn(ip(j)))
      zy=zy+vecr(j+nconf,i)*(-1)**((lj(ih(j))+lj(ip(j)))/2)
     1*ahosc(it1,l1,j1,jj1,nn(ih(j)))*ahosc(it2,l2,j2,jj2,nn(ip(j)))
 2523 continue

      if(dabs(zx).gt.1e-6)then
      zxave(istre_true,icounts2,icounts1)=
     1zxave(istre_true,icounts2,icounts1)+
     2zx*sqbel(i)
      if(istre_true.eq.16)then
c     write(70,2520) i,icounts2,icounts1
c     write(70,2521) zx,sqbel(i),zxave(istre_true,icounts2,icounts1)
      end if 
      ss(i)=ss(i)
     1+radho(jj1,l1,jj2,l2,ltot)*zx*yl(l2,j2,l1,j1,ltot)
c     write(70,777) icounts1,icounts2,yl(l2,j2,l1,j1,ltot),
c    1radho(jj1,l1,jj2,l2,ltot),
c    1yl(l2,j2,l1,j1,ltot)*radho(jj1,l1,jj2,l2,ltot),zx,ss
  777 FORMAT(1X,I2,2X,I2,7(2X,F10.6))
      end if
      if(irpa.eq.1) go to 2530
      if(dabs(zy).gt.1e-6)then
      zyave(istre_true,icounts1,icounts2)=
     1zyave(istre_true,icounts1,icounts2)+
     2zy*sqbel(i)
c     write(7,2520) icounts1,icounts2
c     write(7,2521) zy,dummy
      ss(i)=ss(i)
     1+radho(jj1,l1,jj2,l2,ltot)*zy*yl(l1,j1,l2,j2,ltot)
c     write(70,777) icounts2,icounts1,yl(l1,j1,l2,j2,ltot),
c    1radho(jj1,l1,jj2,l2,ltot),
c    1yl(l1,j1,l2,j2,ltot)*radho(jj1,l1,jj2,l2,ltot),zy,ss
      end if
 2530 if(ltot.ne.1)go to 2529

 2520 format(2i5)
 2521 format(2f10.5)

 2529 continue
 2528 continue
 2527 continue

 2526 continue
 2525 continue
 2524 continue

c     write(7,*) '   '
c     if(ltot.ne.1)write(7,*) '>>> Strength not tested'
c     if(ltot.eq.1)write(7,*) 'State n. ',i,'Strength = ',ss(i)**2

 1001 continue

      dummy1=1.d0

      do 1007 istre=1,nbins
       if(anorma(istre).lt.1e-5)go to 1007 
       anorma(istre)=dsqrt(anorma(istre))
       write(7,*) ''
       write(7,*) 'Energy bin: ',eneinf(istre),'-',enesup(istre)
       write(7,*) ''
       write(7,*) anorma(istre),anorma(istre)**2
       write(7,*) ''
       write(7,2531) dummy1,dummy1
 2531  format(2f10.5)
       anorma_new=0.d0
       do jstre1=1,1
        do jstre2=1,1
         zxave(istre,jstre2,jstre1)=
     1   zxave(istre,jstre2,jstre1)/anorma(istre)
         if(dabs(zxave(istre,jstre2,jstre1)).gt.1e-3)then
         write(7,2520) jstre2,jstre1 
         write(7,2521) zxave(istre,jstre2,jstre1),dummy 
         end if
         anorma_new=anorma_new+zxave(istre,jstre2,jstre1)**2
         zyave(istre,jstre1,jstre2)=
     1   zyave(istre,jstre1,jstre2)/anorma(istre)
         if(irpa.eq.2.and.dabs(zyave(istre,jstre1,jstre2)).gt.1e-3)then
         write(7,2520) jstre1,jstre2
         write(7,2521) zyave(istre,jstre1,jstre2),dummy
         end if
         anorma_new=anorma_new-zyave(istre,jstre1,jstre2)**2
        end do
       end do
       write(7,*) 'Norm: ',anorma_new
1007   continue

1008  continue

ccc   general comments for J.G.
      write(7,*) '   '
      write(7,*) '*************** ---> GENERAL COMMENTS ***************'
      write(7,*) '   '
      write(7,*) 'The operator used has L,S,J given by: '
      write(7,*) ltot,ispin,jtot
      write(7,*) 'Delta hbar-omega is: '
      write(7,*) ndelta

      if(ifon.eq.0)go to 136
      do 2400 i1=1,nconf
      do 2401 i2=i1,nconf
      write(15,7792) aai(i1,i2),bbi(i1,i2)
 2401 continue
 2400 continue
      do 24000 i1=1,nconf
      do 24001 i2=i1,nconf
      write(16,7792) aaiminus(i1,i2),bbiminus(i1,i2)
24001 continue
24000 continue
  136 continue
      close(11)
      if(ifon.ne.0)then 
       close(15)
       close(16)
      end if                                                              
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> To come to the end ',secv,' sec. are needed.'
c     write(2,*) '          '
      stop
 1000 continue
c     tempo=secnds(0.0)
      open(unit=11,status='old',file='rpa.dat')
      if(ifon.ne.0)then 
       open(unit=15,status='old',file='sigma1.dat')
       open(unit=16,status='old',file='sigma2.dat')
      end if 
      if(itot.ne.0.and.iprint.eq.0.and.ishort2.ne.1)open(unit
     1=50,status='old',file='isospin.dat')
      if(ichex.eq.1)nconf=nconft
      do 2402 i1=1,nconf
      do 2403 i2=i1,nconf
      if(ifon.eq.0)go to 2406
 2404 read(15,7792) aai(i1,i2),bbi(i1,i2)
      if (ireal.eq.2) aai(i1,i2)=0.d0
      aai(i2,i1)=aai(i1,i2)
      bbi(i2,i1)=bbi(i1,i2)
      read(16,7792) aaiminus(i1,i2),bbiminus(i1,i2)
      if (ireal.eq.2) aaiminus(i1,i2)=0.d0
      aaiminus(i2,i1)=aaiminus(i1,i2)
      bbiminus(i2,i1)=bbiminus(i1,i2)
 2406 if(itot.eq.0.or.iprint.ne.0.or.ishort2.eq.1)go to 2407
      read(50,*) idum1,idum2,tt(i1,i2)
      tt(i2,i1)=tt(i1,i2)
ccx   read(50,*) tt1(i1,i2),tt2(i1,i2),tt3(i1,i2),tt4(i1,i2)
 2407 continue
 2403 continue
 2402 continue
      if(itot.ne.0.and.iprint.eq.0.and.ishort2.ne.1)close(50)
      read(11,725) s0(1),s0(2),ngrg(1),ngrg(2),nc
      read(11,115) afa0p,noscp,irpap
      ngrf=nc
      if(irpa.eq.2) ngrf=nc/2
      if(dabs(afa0p-afa0).lt.1.d-06) go to 5560
 5561 write(2,5562) afa0,afa0p,nosc,noscp,irpa,irpap
      write(2,5563) ngrg(1),ngrg(2),nc
      stop
 5560 continue
      ngrg3=ngrg(1)+ngrg(2)
      if((nosc.ne.noscp).or.(irpa.ne.irpap).or.(ngrg3.ne.ngrf)) 
     1go to 5561
      do j=1,2
      iikk(j)=0  
      end do
      do 1003 i=1,ngrg3
      read(11,101) evr(i),sqbel(i),bel(i),charge(i)
      read(11,*) (vecr(j,i),j=1,nc)
 1003 continue
      close(11)
      if(ifon.ne.0)then 
       close(15)
       close(16)
      end if  
      open(unit=20,status='new',file='indices.dat')
      ik=0
 1006 ik=ik+1
      do 1005 i=1,ngrg3
      if(ichex.eq.1.and.ik.eq.1.and.charge(i).le.0.d0)go to 1005
      if(ichex.eq.1.and.ik.eq.2.and.charge(i).gt.0.d0)go to 1005
      x=bel(i)/s0(ik)
      if(x.lt.cutoff)go to 1005
      iikk(ik)=iikk(ik)+1
      indx=iikk(ik)
      nig(indx,ik)=i
      write(20,*) ik,indx,nig(indx,ik)
 1005 continue
      if(ichex.eq.1.and.ik.eq.1)go to 1006
      do j=1,2
      ngr(j)=iikk(j)
      end do
 1002 continue
      close(20)  
c     secv=secnds(tempo)
c     write(2,*) '          '
c     write(2,*) '>>> Up to calling of SPREAD ',secv,' sec. are needed.'
c     write(2,*) '          '
      do 1004 jene=1,nene
      ein=ein0+float(jene-1)*dene
      write(2,*) 'ichex = ',ichex
      call spread(ein,isw)
 1004 continue
      stop
      end

      subroutine coul(rho,eta,minl,maxl,fc,fcp,gc,gcp,accur,step)               
      implicit double precision (a-h,o-z)                                       
      double precision k,k1,k2,k3,k4,m1,m2,m3,m4                                
      dimension fc(1),fcp(1),gc(1),gcp(1)                                       
c*****Coulomb wavefunctions calculated at r=rho by the                          
c*****continued-fraction method off steed minl,maxl are actual l-values         
c*****see Barnett Feng Steed and Goldfarb Computer Physics Commun 1974          
      pace = step                                                               
      acc  = accur                                                              
      if(pace.lt.100.d+00) pace = 100.d+00                                      
      if(acc.lt.1.0d-15.or.acc.gt.1.0d-6) acc=1.0d-6                           
      r=rho                                                                   
      ktr=1                                                                     
      lmax=maxl                                                                 
      lmin1=minl+1                                                              
      xll1=minl*lmin1                                                           
      eta2=eta*eta                                                              
      turn = eta + dsqrt(eta2 + xll1)                                           
      if(r.lt.turn.and.dabs(eta).ge.1.d-6)ktr = -1                              
      ktrp = ktr                                                                
      go to 2                                                                   
    1 r=turn                                                                  
      tf=f                                                                    
      tfp=fp                                                                   
      lmax=minl                                                                 
      ktrp=1                                                                    
    2 etar=eta*r                                                                
      rho2=r*r                                                                
      pl=lmax+1                                                                 
      pmx=pl+0.5d+00                                                            
c*****continued fraction for fp(maxl)/f(maxl) xl is f xlprime is fp           
      fp = eta/pl +pl/r                                                        
      dk=etar*2.d+00                                                            
      del=0.0d+00                                                              
      d=0.0d+00                                                              
      f=1.0d+00                                                              
      k=(pl*pl - pl + etar)*(2.0d+00*pl - 1.0d+00)                           
      if(pl*pl+pl+etar.ne.0.d+00) go to 3                                       
      r=r + 1.0d-6                                                           
      go to 2                                                                   
    3 h = (pl*pl + eta2)*(1.0d+00 - pl*pl)*rho2                               
      k = k + dk + pl*pl*6.0d+00                                              
      d = 1.0d+00/(d*h + k)                                                   
      del=del*(d*k - 1.0d+00)                                                  
      if(pl.lt.pmx) del = -r*(pl*pl + eta2)*(pl + 1.0d+00)*d/pl                 
      pl = pl + 1.0d+00                                                         
      fp = fp + del                                                             
      if(d.lt.0d+00) f=-f                                                       
      if(pl.gt.2.d+04) go to 11                                                 
      if(dabs(del/fp).ge.acc) go to 3                                           
      fp = f*fp                                                                
      if(lmax.eq.minl) go to 5                                                 
      fc (lmax+1)=f                                                             
      fcp(lmax+1)=fp                                                            
c*****downward recursion to minl for f and fp, arrays gc,gcp are storage        
      l=lmax                                                                  
      do 4 lp =lmin1,lmax                                                       
      pl=l                                                                      
      gc(l+1) = eta/pl + pl/r                                                  
      gcp(l+1) = dsqrt(eta2 + pl*pl)/pl                                         
      fc (l) = (gc(l+1)*fc(l+1) + fcp(l+1))/gcp(l+1)                          
      fcp(l) = gc(l+1)*fc(l)  - gcp(l+1)*fc(l+1)                              
    4 l=l - 1                                                                  
      f= fc (lmin1)                                                            
      fp= fcp(lmin1)                                                            
    5 if(ktrp.eq.-1) go to 1                                                    
c*****repeat for r = turn if rho lt turn                                        
c*****now obtain p + i.q for minl from continued fraction (32)                  
c*****real arithmetic to facilitate conversion to ibm using real*8              
      p  =0.0d+00                                                               
      q  =r-eta                                                                 
      pl =0.0d+00                                                               
      ar =-(eta2 + xll1)                                                        
      ai = eta                                                                 
      br =2.0d+00*q                                                             
      bi =2.0d+00                                                               
      wi =2.0d+00*eta                                                           
      dr =  br/(br*br + bi*bi)                                                  
      di = -bi/(br*br + bi*bi)                                                  
      dp = -(ar*di + ai*dr)                                                     
      dq =  (ar*dr - ai*di)                                                     
    6 p  =  p + dp                                                              
      q  =  q + dq                                                              
      pl = pl + 2.0d+00                                                         
      ar = ar + pl                                                              
      ai = ai + wi                                                              
      bi = bi + 2.0d+00                                                         
      d  = ar*dr - ai*di + br                                                   
      di = ai*dr + ar*di +bi                                                    
      t = 1.0d+00/(d*d+di*di)                                                   
      dr = t*d                                                                  
      di = -t*di                                                                
      h = br*dr - bi*di - 1.0d+00                                               
      k = bi*dr + br*di                                                         
      t = dp*h - dq*k                                                           
      dq = dp*k + dq*h                                                          
      dp = t                                                                    
      if(pl.gt.4.6d+04) go to 11                                                
      if(dabs(dp)+dabs(dq).ge.(dabs(p)+dabs(q))*acc) go to 6                    
      p = p/r                                                                   
      q = q/r                                                                   
c**** solve for fp,g,gp and normalise f at l=minl                               
      g = (fp - p*f)/q                                                          
      gp = p*g - q*f                                                            
      w = 1.0d+00/dsqrt(fp*g - f*gp)                                            
      g = w*g                                                                   
      gp = w*gp                                                                 
      if (ktr.eq.1) go to 8                                                     
      f = tf                                                                    
      fp = tfp                                                                  
      lmax=maxl                                                                
c**** Runge-Kutta integration of g(minl) and gp(minl) inwards from turn         
c****             see Fox and Mayers 1968 pg 202                                
      if(rho.lt.0.2d+00*turn) pace = 999.d+00                                   
      r3=1.d+00/3.d+00                                                          
      h = (rho - turn)/(pace + 1.0d+00)                                         
      h2 = 0.5d+00*h                                                            
      i2 = ifix(real(pace + 0.001))                                            
      etah = eta*h                                                              
      h2ll = h2*xll1                                                            
      s = (etah + h2ll/r  )/r  - h2                                             
    7 rh2=r + h2                                                               
      t = (etah + h2ll/rh2)/rh2 - h2                                            
      k1 = h2*gp                                                                
      m1 =  s*g                                                                 
      k2 = h2*(gp + m1)                                                         
      m2 =  t*(g  + k1)                                                         
      k3 =  h*(gp + m2)                                                         
      m3 =  t*(g  + k2)                                                         
      m3 =     m3 + m3                                                          
      k4 = h2*(gp + m3)                                                         
      rh = r + h                                                                
      s  = (etah + h2ll/rh )/rh  - h2                                           
      m4 =  s*(g + k3)                                                          
      g  = g  + (k1 + k2 + k2 + k3 +k4)*r3                                      
      gp = gp + (m1 + m2 + m2 + m3 + m4)*r3                                     
      r = rh                                                                    
      i2 = i2 - 1                                                               
      if(i2.ge.0) go to 7                                                       
      w = 1.0d+00/(fp*g - f*gp)                                                 
c***  upward recursion from gc(minl) and gcp(minl),stored values are r,s        
c***  renormalise fc,fcp for each l-value                                       
    8 gc (lmin1) = g                                                            
      gcp(lmin1) = gp                                                           
      if(lmax.eq.minl) go to 10                                                 
      do 9 l = lmin1,lmax                                                       
      t      = gc(l+1)                                                          
      gc (l+1) = (gc(l)*gc (l+1) - gcp(l))/gcp(l+1)                             
      gcp(l+1) = gc(l)*gcp(l+1) - gc(l+1)*t                                     
      fc (l+1) = w*fc (l+1)                                                     
    9 fcp(l+1) = w*fcp(l+1)                                                     
      fc (lmin1) = fc (lmin1)*w                                                 
      fcp(lmin1) = fcp(lmin1)*w                                                 
      return                                                                    
   10 fc (lmin1) = w*f                                                          
      fcp(lmin1) = w*fp                                                         
      return                                                                    
   11 w = 0.0d+00                                                               
      g = 0.0d+00                                                               
      gp = 0.0d+00                                                              
      go to 8                                                                   
      end                                                                       

      function fa(l,x)
      implicit real*8 (a-h,o-z)
      dimension fef(0:20)
      fef(0)=dsin(x)
      fef(1)=dsin(x)/(x)-dcos(x)
      if(l.ne.0) go to 1
      fa=fef(0)
      go to 3
    1 if(l.ne.1) go to 2
      fa=fef(1)
      go to 3
    2 continue
      do 4 le=2,l
      le1=le-1
      le2=le-2
      fef(le)=fef(le1)*(2.d0*le-1.d0)/(x)-fef(le2)
    4 continue
      fa=fef(l)
    3 return
      end

      function ga(l,x)
      implicit real*8 (a-h,o-z)
      dimension gg(0:20)
      gg(0)=dcos(x)
      gg(1)=dcos(x)/(x)+dsin(x)
      if(l.ne.0) go to 1
      ga=gg(0)
      go to 3
    1 if(l.ne.1) go to 2
      ga=gg(1)
      go to 3
    2 continue
      do 4 le=2,l
      le1=le-1
      le2=le-2
      gg(le)=gg(le1)*(2.d0*le-1.d0)/(x)-gg(le2)
    4 continue
      ga=gg(l)
    3 return
      end

      subroutine lecpot                                                         
      implicit real*8 (a-h,o-z)
      parameter (nnn=180)
      character*80 file_hf 
      common/bfile1/file_hf                                                         
      common/bpot/xmn(nnn),xmp(nnn),vn(nnn),vp(nnn),yn(nnn),yp(nnn),            
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)                               
      common/bpas/h,cof                                                         
      common/bnri/nri
      common/blecp1/t0,t1,t2,t3,t13,alfe,x0,x1,x2,x3,anucl,znucl
      common/blecp2/dt(nnn),taut(nnn),dt1(nnn),dt2(nnn),dp(nnn),dn(nnn)         
      dimension dal(3),mash(3)                                                  
c 208 format(6e12.5)                                                            
c inconsistent with HF ?
  208 format(6e12.5)                                                            
  209 format(15i3)
      open(unit=10,status='old',file='fort.10')   
      read(10,*) ndummy
      read(10,208) t0,t1,t2,t3,t13,x0         
      read(10,208) alfe,x3,x1,x2,w2                  
      read(10,208) anucl,znucl,h,(dal(i),i=1,3)   
      read(10,209)(mash(i),i=1,3)                 
      read(10,209) i1,i2,i3,i4,i5                 
      read(10,208) (xmn(j),j=1,nri),(xmp(j),j=1,nri) 
      read(10,208) ( vn(j),j=1,nri),( vp(j),j=1,nri) 
      do 400 j=1,nri    
      vn(j)=cof*vn(j)                                                           
  400 vp(j)=cof*vp(j)                                                           
      read(10,208) (yn(j),j=1,nri),(yp(j),j=1,nri)   
      read(10,208) (vsn(j),j=1,nri),(vsp(j),j=1,nri) 
      read(10,208)(vc(j),j=1,nri)                    
      read(10,208)(dt(j),j=1,nri),(taut(j),j=1,nri)  
      read(10,208) (dt1(j),j=1,nri),(dt2(j),j=1,nri) 
      read(10,208)(dp(j),j=1,nri),(dn(j),j=1,nri)   
      n2=nri-2                                                                  
      wnd(1)=(-11.d0*yn(1)+18.d0*yn(2)-9.d0*yn(3)+2.d0*yn(4))/
     1(6.d0*h)                    
      wpd(1)=(-11.d0*yp(1)+18.d0*yp(2)-9.d0*yp(3)+2.d0*yp(4))/
     1(6.d0*h)                    
      do 1 i=2,n2    
      wnd(i)=(-2.d0*yn(i-1)-3.d0*yn(i)+6.d0*yn(i+1)-yn(i+2))/
     1(6.d0*h)                   
      wpd(i)=(-2.d0*yp(i-1)-3.d0*yp(i)+6.d0*yp(i+1)-yp(i+2))/
     1(6.d0*h)                   
    1 continue
      wnd(nri-1)=wnd(n2)                                                        
      wnd(nri)=wnd(n2)                                                          
      wpd(nri-1)=wpd(n2)                                                        
      wpd(nri)=wpd(n2)                        
      close(10)                                   
      return
      end                                                                    

      subroutine numnl(l,jlj,sreg,sirg,defa) 
ccc   c2d=cos(2*delta),s2d=sin(2*delta),sigma=elastic cross-section
ccc   defa=(u*dw-du*w)*(hb/2m*) if negative energy                              
      implicit real*8 (a-h,o-z)
      parameter (nnn=180)
      common/pot/u(nnn),rmass,eta,xka,xmsm(nnn),sig, 
     1sgne,dreg(nnn),dirg(nnn),c2d,s2d,sigma  
      common/bpas/h,cof                                                         
      common/temp/i2
      common/tempr/xjp,omega                                                  
      common/bjost/nr
      common/bjostr/refp,aifp,refd,aifd
      common/bnri/nri
      common/bnrir/dr
      dimension fcd(100),fcpd(100),gcd(100),gcpd(100),f(200),g(200)
      dimension sreg(nnn),sirg(nnn)             
  100 format(/2x,'potentiel v-k**2'/20('*'))                                    
  101 format(10e12.5)                                                           
  102 format(/2x,'solution ksi(r)'/20('*'))                                     
  103 format(/2x,'solution u(r) non normee'/26('*'))                            
  104 format(2x,'facteur dsqrt(2m/(pi*hb2*k))=',e12.5/30('*'))                  
  105 format(2x,'wronskien=',e12.5/22('*'))                                     
  106 format(2x,'wronskien calcule'/20('*'))                                    
      pi=4*datan(1.d0)
      xka2=xka*xka*sgne                                                         
      spd=100.d0                                                                
      epsd=1.d-10                                                               
      n1=nri-1                                                                  
      n2=nri-2                                                                  
      do 1 j=1,nri                                                             
    1 u(j)=u(j)-xka2                                                            
      h2=h*h                                                                    
      s3=dsqrt(3.d0) 
      s12=dsqrt(12.d0) 
      sig=0.d0                                                                  
      sigl1=0.d0                                                                
      lmax=40                                                                   
      et2=eta*eta                                                               
      etad=eta                                                            
      l1=l+1                                                                    
      al=float(lmax-1)                                                          
      ll=lmax                                                                   
      c1=al+1.d0                                                                
      if(eta-1.d-36) 8,8,7                                                      
    7 phi2=c1*c1+et2                                                            
      phi=dsqrt(phi2)                                                           
      chi=atan(eta/c1)                                                          
      sigl=chi*(al+0.5d0)+eta*(dlog(phi)-1.d0)-(sin(chi)/12.d0-                 
     1(sin(3.d0*chi)/360.d0-(sin(5.d0*chi)/1260.d0-
     2(sin(7.d0*chi)/1680.d0-(sin(9.d0*chi)/1188.d0)/phi2)/
     3phi2)/phi2)/phi2)/phi                           
    9 sigl1=sigl-atan(eta/al)                                                   
      sigl=sigl1                                                                
      al=al-1.d0                                                                
      ll=ll-1                                                                   
      if(ll.eq.l1) sig=sigl                                                     
      if(ll.ge.2) go to 9                                                       
    8 continue                                                                  
      xl=float(l)                                                               
      a=1.d-07                                                                  
      c=u(2)-xl*(xl+1.d0)/(4.d0*h2)                                            
      sreg(1)=a*(h**l1)*(1.d0+c*h2/(4.d0*xl+6.d0))                              
      sreg(1)=sreg(1)*(1.d0-xl*(xl+1.d0)/12.d0-c*h2/12.d0)                      
      sreg(2)=a*((2.d0*h)**l1)*(1.d0+c*4.d0*h2/(4.d0*xl+6.d0))                  
      sreg(2)=sreg(2)*(1.d0-h2*u(2)/12.d0)                                      
    3 do 4 j=2,n1                                                            
    4 sreg(j+1)=(2.d0+12.d0*h2*u(j)/(12.d0-h2*u(j)))*sreg(j)
     1-sreg(j-1)                
      do 5 j=1,nri                                                           
    5 sreg(j)=sreg(j)/(1.d0-h2*u(j)/12.d0)                                      
      sreg(1)=sreg(1)/(1.d0-xl*(xl+1.d0)/12.d0-c*h2/12.d0)                      
      rm=h*(nri-3.d0)                                                           
      rkr=rm*xka                                                                
      ib=1                                                                      
      if((eta.gt.25.d0).and.(rkr.lt.10.d0)) ib=2     
      xim=rm/h+0.0001d0                                                        
      im=int(xim)                                                               
      xnorm=3.14159265d0*rmass*xka/0.04819d0                                  
      xnorm=1.d0/dsqrt(xnorm)                          
      im1=im-1                                                                  
      rm1=h*float(im1)                                                          
      rkr1=rm1*xka                                                              
      if(sgne) 30,30,31                                                         
   31 continue                                                                  
      do 99 kkk=1,16                                                          
      f(kkk)=0.d0                                                               
   99 g(kkk)=0.d0                                                               
      l2=l+2                                                                    
      if(l) 200,201,200                                                         
  201 l3=l                                                                      
      l4=l1                                                                     
      go to 202                                                                 
  200 l3=l-1                                                                    
      l4=l                                                                      
  202 continue                                                                  
      if(eta-1.d-06) 97,97,96                                                   
   97 lll=l3-1                                                                  
      do 95 kkk=1,2                                                           
      lll=lll+1                                                                 
      ll1=lll+1                                                                 
      f(ll1)=fa(lll,rkr)                                                 
   95 g(ll1)=ga(lll,rkr)                                                
      go to 90                                                                  
   96 go to (1001,1003),ib                                                      
 1003 eta=0.d0                                                                  
      et2=0.d0                                                                  
      go to 97                                                                  
 1001 rod=rkr                                                             
      call coul(rod,etad,l3,l4,fcd,fcpd,gcd,gcpd,epsd,spd)                      
      l31=l3+1                                                                  
      l32=l3+2                                                                  
      do 203 j=l31,l32                                                       
      f(j)=fcd(j)                                                         
  203 g(j)=gcd(j)                                                         
   90 continue                                                                  
      xf=f(l1)                                                                  
      xg=g(l1)                                                                  
      xl2=float(l*l)                                                            
      if(l) 10,11,10                                                            
   10 xfp=xka*(dsqrt(xl2+et2)*f(l)-(eta+xl2/rkr)*xf)/float(l)                  
      xgp=xka*(dsqrt(xl2+et2)*g(l)-(eta+xl2/rkr)*xg)/float(l)                   
      go to 12                                                                  
   11 xfp=xka*((eta+1.d0/rkr)*xf-dsqrt(1.d0+et2)*f(l2))                         
      xgp=xka*((eta+1.d0/rkr)*xg-dsqrt(1.d0+et2)*g(l2))                         
   12 wrsk=(xfp*xg-xf*xgp)/xka                                                  
      wrs=abs(wrsk-1.d0)                                                        
      if(wrs-1.d-04) 14,14,13                                                   
   13 continue                                                                  
      write(2,1000) i2,xjp,l,eta,rkr,omega                                        
 1000 format(2x,'i2=',i3,3x,'xjp=',f5.2,3x,'l=',i3,3x,'eta=',e12.5,             
     13x,'rkr=',e12.5,3x,'omega=',e12.5)                                        
      write(2,300) (f(kkk),kkk=1,6)                                               
      write(2,300) (g(kkk),kkk=1,6)                                               
      write(2,300) xka                                                            
  300 format(6e12.5)                                                           
      write(2,105) wrsk                                                           
   14 continue                                                                  
      dre=(-2.d0*sreg(im-1)-3.d0*sreg(im)+6.d0*sreg(im+1)-sreg(im+2))
     1/(6.d0*h)          
      xa=sreg(im)*xfp-dre*xf                                                    
      xb=dre*xg-sreg(im)*xgp                                                    
      xab=xa/xb                                                                 
      tg=2.d0*xa/(xb-xa*xab)                                                    
      defa=0.5d0*atan(tg)                                                       
      sg=1.d0                                                                   
      c2d=0.d0                                                                  
      s2d=0.d0                                                                  
      sigma=0.d0                                                                
      dxa=xa                                                              
      dxb=xb                                                              
      if(sgne.le.0.d0) go to 205                                                
      dxba=dxb**2+dxa**2                                                        
      dc2d=(dxb**2-dxa**2)/dxba                                                 
      ds2d=2.d0*dxa*dxb/dxba                                                    
      dqs=dsqrt(dxba)                                                           
      dcod=dxb/dqs                                                              
      dsid=dxa/dqs                                                              
      c2d=dc2d                                                            
      s2d=ds2d                                                            
      cod=dcod                                                            
      sid=dsid                                                            
      sigma=pi*(jlj+1.d0)*((1.d0-c2d)**2+s2d**2)/(4.d0*xka2)                    
  205 continue                                                                  
      if(xb.lt.0.d0) sg=-1.d0
      xlam=3.14159265d0*rmass*xka*(1.d0+xab**2)/0.04819d0                       
      xlam=sg*(xf+xab*xg)/(sreg(im)*dsqrt(xlam))                                
      go to 32                                                                  
   30 xlam=1.d0                                                                 
   32 continue                                                                  
      do 6 j=1,nri                                                           
    6 sreg(j)=sreg(j)*xlam*dsqrt(xmsm(j))                                       
ccc   Integration of irregular solution                                         
      if(sgne) 34,34,33                                                         
   33 continue                                                                  
      do 98 kkk=1,16                                                         
      f(kkk)=0.d0                                                               
   98 g(kkk)=0.d0                                                               
      if(eta-1.d-06) 94,94,93                                                   
   94 f(l1)=fa(l,rkr1)                                                  
      g(l1)=ga(l,rkr1)                                               
      go to 89                                                                  
   93 go to (1002,94),ib                                                        
 1002 rod=rkr1                                                            
      call coul(rod,etad,l,l,fcd,fcpd,gcd,gcpd,epsd,spd)                        
      f(l1)=fcd(l1)                                                       
      g(l1)=gcd(l1)                                                       
   89 continue                                                                  
      xno=(xa/xb)**2                                                            
      xno=1.d0/dsqrt(1.d0+xno)   
      sg=1.d0                  
      if(xb.lt.0.d0) sg=-1.d0
      sirg(im)=sg*xnorm*xno*(xg-xa*xf/xb)*(1.d0-h2*u(im)/12.d0)                 
      sirg(im1)=sg*xnorm*xno*(g(l1)-xa*f(l1)/xb)*(1.d0-h2*u(im1)/12.d0)         
      go to 35                                                                  
   34 sirg(im)=xnorm*(1.d0-h2*u(im)/12.d0)/exp(rkr)                             
      sirg(im1)=xnorm*(1.d0-h2*u(im1)/12.d0)/exp(rkr1)                          
   35 continue                                                                  
      do 15 j=im1,3,-1                                                       
   15 sirg(j-1)=(2.d0+12.d0*h2*u(j)/(12.d0-h2*u(j)))*sirg(j)-sirg(j+1)          
      do 16 j=im,n1                                                           
   16 sirg(j+1)=(2.d0+12.d0*h2*u(j)/(12.d0-h2*u(j)))*sirg(j)-sirg(j-1)          
      do 17 j=1,nri                                                           
   17 sirg(j)=dsqrt(xmsm(j))*sirg(j)*12.d0/(12.d0-h2*u(j))                      
      if(l) 20,19,20                                                            
   19 den=1.d0                                                                  
      den1=1.d0                                                                 
      go to 21                                                                  
   20 den=(2.d0*h)**l                                                           
      den1=h**l                                                                 
   21 afa=(1.d0+4.d0*c*h2/(2.d0-4.d0*xl))/den                                   
      afa=sirg(2)/afa                                                           
      sirg(1)=afa*(1.d0+c*h2/(2.d0-4.d0*xl))/den1                               
ccc   Test du wronskien                                                         
      if(l.le.100) go to 40                                                     
      u(1)=0.d0                                                                 
      do 18 j=2,n2                                                           
      f(j)=(-2.d0*sreg(j-1)-3.d0*sreg(j)+6.d0*sreg(j+1)-sreg(j+2))
     1/(6.d0*h)             
      g(j)=(-2.d0*sirg(j-1)-3.d0*sirg(j)+6.d0*sirg(j+1)-sirg(j+2))
     1/(6.d0*h)             
   18 u(j)=(sreg(j)*g(j)-sirg(j)*f(j))/xmsm(j)                                  
      write(2,106)                                                                
      write(2,101) (u(j),j=1,n2)                                                  
   40 continue                                                                  
ccc   Calcul des derivees.                                                      
      d=1.d0/(24.d0*h)                                                          
      do 41 i=1,nri                                                            
      if(i-1) 42,42,43                                                          
   42 dreg(i)=(-50.d0*sreg(1)+96.d0*sreg(2)-72.d0*sreg(3)+32.d0                 
     1*sreg(4)-6.d0*sreg(5))*d                                                 
      dirg(1)=(-50.d0*sirg(1)+96.d0*sirg(2)-72.d0*sirg(3)+32.d0                 
     1*sirg(4)-6.d0*sirg(5))*d                                                  
      go to 41                                                                  
   43 if(i-2) 44,44,45                                                          
   44 dreg(2)=(-6.d0*sreg(1)-20.d0*sreg(2)+36.d0*sreg(3)-12.d0                  
     1*sreg(4)+2.d0*sreg(5))*d                                                  
      dirg(2)=(-6.d0*sirg(1)-20.d0*sirg(2)+36.d0*sirg(3)-12.d0                  
     1*sirg(4)+2.d0*sirg(5))*d                                                  
      go to 41                                                                  
   45 if(i-nri) 47,46,46                                                        
   46 dreg(nri)=(50.d0*sreg(i)-96.d0*sreg(i-1)+72.d0*sreg(i-2)                 
     1-32.d0*sreg(i-3)+6.d0*sreg(i-4))*d                                        
      dirg(nri)=(50.d0*sirg(i)-96.d0*sirg(i-1)+72.d0*sirg(i-2)                  
     1-32.d0*sirg(i-3)+6.d0*sirg(i-4))*d                                        
      go to 41                                                                  
   47 if(i-nri+1) 49,48,48                                                      
   48 dreg(i)=(6.d0*sreg(i+1)+20.d0*sreg(i)-36.d0*sreg(i-1)                     
     1+12.d0*sreg(i-2)-2.d0*sreg(i-3))*d                                        
      dirg(i)=(6.d0*sirg(i+1)+20.d0*sirg(i)-36.d0*sirg(i-1)                     
     1+12.d0*sirg(i-2)-2.d0*sirg(i-3))*d                                        
      go to 41                                                                  
   49 dreg(i)=(2.d0*(sreg(i-2)-sreg(i+2))-16.d0*(sreg(i-1)-
     1sreg(i+1)))*d            
      dirg(i)=(2.d0*(sirg(i-2)-sirg(i+2))-16.d0*(sirg(i-1)-
     1sirg(i+1)))*d            
   41 continue                                                                  
      if(sgne) 111,111,110                                                      
  111 rwk=4.d0                                                                  
      iwk=int(rwk/h)                                                            
      defa=sreg(iwk)*dirg(iwk)-dreg(iwk)*sirg(iwk)                              
      defa1=rmass/(0.04819d0*xmsm(iwk))                                         
      defa=defa*defa1                                                           
  110 continue                                                                  
      coq=0.04819d0/(3.14159265d0*rmass*xka)                                    
      refp=(sirg(nr)*cod+sreg(nr)*sid)/dsqrt(coq*xmsm(nr))                      
      aifp=(sreg(nr)*cod-sirg(nr)*sid)/dsqrt(coq*xmsm(nr))                      
      sp1=dsqrt(xmsm(nr+1))                                                     
      sp2=dsqrt(xmsm(nr+2))                                                     
      sm1=dsqrt(xmsm(nr-1))                                                     
      sm2=dsqrt(xmsm(nr-2))                                                     
      dru=(2.d0*(sreg(nr-2)/sm2-sreg(nr+2)/sp2)-16.d0*(sreg(nr-1)
     1/sm1-sreg(nr+1)/sp1))*d                                                   
      diu=(2.d0*(sirg(nr-2)/sm2-sirg(nr+2)/sp2)-16.d0*(sirg(nr-1)         
     1/sm1-sirg(nr+1)/sp1))*d                                                   
      refd=(diu   *cod+dru     *sid)/dsqrt(coq)                                
      aifd=(dru   *cod-diu     *sid)/dsqrt(coq)                                
      return                                                                    
      end                                                                       

      subroutine phme(ia,ib,ic,id,jt,nr,dr,res,resm,induce)
      implicit real*8 (a-h,o-z)
c     real*4 sla,slb,slc,sld,sja,sjb,sjc,sjd,sjt,c1,c2,c3,c4,
c    1sl,slk1,slk2,sll,slka,slkb,slkc,slkd
      parameter (nsp=300,nnn=180)
      common/betalf/calf(7),cbet(7)
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/brintw/riw,riwminus,fw(nnn),fwminus(nnn)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      dimension hh(9),vv(7,3)
   52 format(5x,'res=',e12.5)
      res=0.d0
      djp1=2.d0*jt+1.d0
      la=ll(ia)
      lb=ll(ib)
      lc=ll(ic)
      ld=ll(id)
      ja=lj(ia)
      jb=lj(ib)
      jc=lj(ic)
      jd=lj(id)
      indu=induce+1
      go to (1000,2000), indu
 1000 continue
      do 1 i=1,9
    1 hh(i)=0.d0
      call rdint1_g(ia,ib,ic,id,nr,dr)

      sla=dfloat(la)
      slb=dfloat(lb)
      slc=dfloat(lc)
      sld=dfloat(ld)
      sja=dfloat(ja)/2.d0
      sjb=dfloat(jb)/2.d0
      sjc=dfloat(jc)/2.d0
      sjd=dfloat(jd)/2.d0
      sjt=dfloat(jt)
      hh(8)=ri0*yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/djp1
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
    3 continue
    2 hh(9)=ris*x9/djp1
      l1=max0(ladm,lbcm)
      l2=min0(ladp,lbcp)
      if(l1.gt.l2) go to 4
      do 5 i1=1,7
      do 6 i2=1,3
      vv(i1,i2)=0.d0
    6 continue
    5 continue
      cofs=(ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0)
      cofs=dsqrt(cofs)*6.d0
      is=la+lc+(jb+jd)/2+1
      si=1.d0-2.d0*mod(is,2)
      call sixj(sja,sla,0.5d0,sld,sjd,sjt,c1)
      call sixj(sjc,slc,0.5d0,slb,sjb,sjt,c2)
      c1=c1*c2
c     cof0=si*cofs*dble(c1)/6.d0
      cof0=si*cofs*c1/6.d0
ccc   calculation of vl(i) starts here.
      ilpt=1
      ilgd=3
      do 7 il=ilpt,ilgd
      l=jt-2+il
      if(l.lt.0) go to 7
      if((l-l1)*(l-l2).gt.0) go to 7
      vv(1,il)=rx1*redy(la,l,ld)*redy(lc,l,lb)/(2.d0*l+1.d0)
    7 continue
      do 8 iup=2,5
      go to (9,9,10,11,12), iup
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
   12 ka=ib
      kb=id
      kc=ia
      kd=ic
   13 lka=ll(ka)
      lkb=ll(kb)
      lkc=ll(kc)
      lkd=ll(kd)
      slka=dfloat(lka)
      slkb=dfloat(lkb)
      slkc=dfloat(lkc)
      slkd=dfloat(lkd)
      call rdint2_g(ka,kb,kc,kd,nr,dr)
      xlkab=(2.d0*lka+1.d0)*(2.d0*lkb+1.d0)
      xlkab=dsqrt(xlkab)
      do 14 il=ilpt,ilgd
      l=jt-2+il
      sl=dfloat(l)
      if(l.lt.0) go to 14
      if((l-l1)*(l-l2).gt.0) go to 14
      sgn1=1.d0-2.d0*mod(l,2)
      if(iup.eq.4.or.iup.eq.5) sgn1=1.d0
      do 15 ik1=1,2
      k1=lka-3+2*ik1
      if(k1.lt.0) go to 15
      slk1=dfloat(k1)
c     write(21,*) 'Outside TROISJ'                                       
c     write(21,*) 'we have call troisj(slk1,slka,1.,0.,0.,0.,c1)'
c     write(21,*) slk1,slka
      call troisj(slk1,slka,1.d0,0.d0,0.d0,0.d0,c1)
c     cofk1=(1.d0-2.d0*mod(k1,2))*dsqrt(2.d0*k1+1.d0)*dble(c1)
      cofk1=(1.d0-2.d0*mod(k1,2))*dsqrt(2.d0*k1+1.d0)*c1 
      do 16 ik2=1,2
      k2=lkb-3+2*ik2
      if(k2.lt.0) go to 16
      slk2=dfloat(k2)
      call troisj(slk2,slkb,1.d0,0.d0,0.d0,0.d0,c2)
c     cofk2=(1.d0-2.d0*mod(k2,2))*dsqrt(2.d0*k2+1.d0)*dble(c2)
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
      sll=dfloat(l1l)
      call sixj(sll,1.d0,sl,slka,slkd,slk1,c1)
      call sixj(sll,1.d0,sl,slkb,slkc,slk2,c2)
      c3=c1*c2
c     vv(iup,il)=vv(iup,il)+sgn1*xlkab*cofk1*cofk2*dble(c3)*
c   1redy(lkd,l1l,k1)*redy(lkc,l1l,k2)*rx2(ik1,ik2)
       vv(iup,il)=vv(iup,il)+sgn1*xlkab*cofk1*cofk2*c3*
     1redy(lkd,l1l,k1)*redy(lkc,l1l,k2)*rx2(ik1,ik2)
   17 continue
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
      slka=dfloat(lka)
      slkb=dfloat(lkb)
      slkc=dfloat(lkc)
      slkd=dfloat(lkd)
      call rdint2_g(ka,kd,kb,kc,nr,dr)
      xlkad=(2.d0*lka+1.d0)*(2.d0*lkd+1.d0)
      xlkad=dsqrt(xlkad)
      do 20 il=ilpt,ilgd
      l=jt-2+il
      sl=dfloat(l)
      if(l.lt.0) go to 20
      if((l-l1)*(l-l2).gt.0) go to 20
      cofl=redy(lkb,l,lkc)/(2.d0*l+1.d0)
      do 21 ik1=1,2
      k1=lka-3+2*ik1
      if(k1.lt.0) go to 21
      slk1=dfloat(k1)
      call troisj(slk1,slka,1.d0,0.d0,0.d0,0.d0,c1)
c     cofk1=dsqrt(2.d0*k1+1.d0)*dble(c1)
      cofk1=dsqrt(2.d0*k1+1.d0)*c1 
      do 22 ik2=1,2
      k2=lkd-3+2*ik2
      if(k2.lt.0) go to 22
      slk2=dfloat(k2)
      call troisj(slk2,slkd,1.d0,0.d0,0.d0,0.d0,c2)
c     cofk2=dsqrt(2.d0*k2+1.d0)*dble(c2)
      cofk2=dsqrt(2.d0*k2+1.d0)*c2
      call sixj(slk2,slk1,sl,slka,slkd,1.d0,c3)
c     cofk12=redy(k1,l,k2)*dble(c3)
      cofk12=redy(k1,l,k2)*c3 
      vv(iup,il)=vv(iup,il)+xlkad*cofl*cofk1*cofk2*cofk12*rx2
     1(ik1,ik2)
   22 continue
   21 continue
   20 continue
   18 continue
      do 23 il=ilpt,ilgd
      l=jt-2+il
      if(l.lt.0) go to 23
      sl=dfloat(l)
      dlp1=2.d0*l+1.d0
      call neufj(sja,sla,0.5d0,sjd,sld,0.5d0,sjt,sl,1.d0,c3)
      call neufj(sjc,slc,0.5d0,sjb,slb,0.5d0,sjt,sl,1.d0,c4)
      c4=c3*c4
      do 24 iup=1,7
c  24 hh(iup)=hh(iup)+cofs*dlp1*dble(c4)*vv(iup,il)
      hh(iup)=hh(iup)+cofs*dlp1*c4*vv(iup,il)
   24 continue   
   23 continue
      do 25 iup=1,7
      hh(iup)=cbet(iup)*hh(iup)+calf(iup)*cof0*vv(iup,2)
c     if(ia.eq.60.and.ic.eq.60.and.id.eq.1.and.ib.eq.1)then
c     if(iup.eq.4)then 
c     write(*,*) cbet(iup),calf(iup),cof0,vv(iup,2)
c     if(iup.eq.4)write(*,*) hh(iup) 
c     end if
c     end if
   25 continue
    4 continue
  102 continue
C-AUERBACH
c     do 26 iup=8,9
c  26 res=res+hh(iup)
      do 26 iup=1,9
      res=res+hh(iup)
   26 continue
      if(ia.eq.60.and.ic.eq.60.and.id.eq.1.and.ib.eq.1)then
      write(*,*) res
c     stop
      end if
C-AUERBACH-END
      resm=0.d0
      go to 3000
 2000 call rdintw(ia,ib,ic,id,nr,dr)
      res=riw      *yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/djp1
      resm=riwminus*yl(la,ja,ld,jd,jt)*yl(lc,jc,lb,jb,jt)/djp1
 3000 continue
      return
      end

      subroutine qwg(b,ns,ls,js,iq,noc) 
      implicit real*8 (a-h,o-z)
      parameter (nnn=180,nmx=12)
      dimension u(nnn,nmx),up(nnn,nmx),d(5,5),vu(nnn),wu(nnn)   
      dimension a(nmx,nmx),v(nmx),t1(nmx),t2(nmx),t3(nmx)
      dimension e(nmx),vec(nmx),nord(nmx)
      common/bhosc/ahosc(2,0:12,1:28,nmx,nmx)
      common/bwg/wg(nnn,nmx),ewg(nmx),xor(nmx),dwg(nnn,nmx)   
      common/bpas/h,cof
      common/bnri/nri
      common/bnrir/dr
      common/hsqr1/t1
      common/hsqr2/t2
      common/hsqr3/t3 
      common/bpot/xmn(nnn),xmp(nnn),vn(nnn),vp(nnn),yn(nnn),yp(nnn),
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
  103 rho=dk*rho
  101 dpl=1.d0
      xsl=1.d0
      if(l) 104,104,105 
  105 k=1   
      do 106 ii=1,l  
      k=k+2 
      xk=k  
      dk=xk
      dpl=2.d0*dpl/dk
  106 xsl=xsl*dx 
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
    4 sum=sum+dbd*u(jj,n1)  
    3 up(j,n1)=sum/(dbh)  
    1 continue  
      xll=l*(l+1.d0)  
      xjj=0.5d0*(0.25d0*js*(js+2.d0)-xll-0.75d0)  
      go to (5,6),iq
    5 do 7 j=1,nri   
      x2=h*j*h*j
      wu(j)=xmn(j)  
    7 vu(j)=xmn(j)*xll/x2+vn(j)+xjj*vsn(j)  
      go to 8   
    6 do 9 j=1,nri   
      x2=h*j*h*j
      wu(j)=xmp(j)  
    9 vu(j)=xmp(j)*xll/x2+vp(j)+xjj*vsp(j)+vc(j)
    8 do 10 i1=1,ns  
      do 11 i2=i1,ns 
      aa=0.d0 
c     if(i1.eq.i2)write(89,*) ls,js,i1
      do 12 j=1,nri  
      rrr=float(j)*h
      dwu=wu(j)
      dvu=vu(j)
      aa=aa+dwu*up(j,i1)*up(j,i2)+dvu*u(j,i1)*u(j,i2)
      if(i1.eq.i2)then
c     write(89,*) float(j)*dr,u(j,i1),hwf(i1-1,ls,b,rrr)
      end if
c     if(ls.eq.0.and.js.eq.1.and.iq.eq.2.and.i1.eq.1.and.i2.eq.6)then
c      write(99,111) dwu,up(j,i1),up(j,i2),u(j,i1),u(j,i2)
c      write(99,*) aa
c 111  format(5(1x,e11.5))
c     end if
   12 continue  
      aa=aa*h
      a(i1,i2)=(aa) 
      a(i2,i1)=a(i1,i2) 
   11 continue  
   10 continue  
c     write(99,*) 'Entering qwg with ls,js,iq = '
c     write(99,*) ls,js,iq
c     if(ls.eq.0.and.js.eq.1.and.iq.eq.2)then
c      do i1=1,ns
c      write(99,110) (a(i1,i2),i2=1,ns)
c      end do
c     end if
c 110 format(18(1x,e12.5))
      call diasym(a,nmx,ns,v,e,vec,nord)
c     if(ls.eq.0.and.js.eq.1.and.iq.eq.2)then
c     write(99,*) 'Diagonalization on the HO basis'
c     write(99,*) 'l,j = ',ls,js
      if(l.gt.12.or.js.gt.28)then
       write(*,*) '>>> INCREASE DIMENSIONS OF AHOSC'
       stop
      end if
      do n=1,ns
       if(n.eq.noc)then
c       write(99,*) 'Eigenvalue: ',v(n)
c       write(99,*) (a(ik,n),ik=1,ns)
        do j=1,ns
        ahosc(iq,ls,js,j,noc)=a(j,n)
        end do
c       write(98,*) 'Eigenvalue: ',v(n)
c       write(98,*) iq,ls,js,' last index ',n
c       write(98,*) (ahosc(iq,ls,js,jpr,n),jpr=1,ns)
       end if
      end do
      do 20 n=1,ns   
      xnorm=0.d0  
      sgn=0.d0
      dy=0.d0
      ewg(n)=v(n) 
      do 50 ik=1,ns  
      dy=dy+a(ik,n)**2
      as=a(ik,n)  
   50 sgn=sgn+as*u(4,ik)
      xor(n)=dy
      sign=1.d0   
      if(sgn.lt.0.d0) sign=-1.d0
      do 21 j=1,nri  
      xx=0.d0 
      xyx=0.d0
      do 22 ii=1,ns  
      xyx=xyx+a(ii,n)*up(j,ii)
   22 xx=xx+a(ii,n)*u(j,ii)  
      wg(j,n)=xx*sign   
      dwg(j,n)=sign*xyx
   21 xnorm=xnorm+(wg(j,n)**2)*h
      xnorm=dsqrt(xnorm) 
      do 23 j=1,nri  
      dwg(j,n)=dwg(j,n)/xnorm
   23 wg(j,n)=wg(j,n)/xnorm
      if(n.eq.noc)then
c      write(78,*) 'Eigenvalue: ',ewg(n)
       do 240 ik=1,ns
       ahosc2=0.d0
       do 230 j=1,nri 
       ahosc2=ahosc2+h*wg(j,n)*u(j,ik)
  230  continue
c      write(78,*) ik,ahosc2
       ahosc(iq,ls,js,ik,noc)=ahosc2
  240  continue
      end if
   20 continue  
      return
      end

      subroutine rdintw(ia,ib,ic,id,nr,dr)
      implicit real*8 (a-h,o-z)
      parameter (nnn=180,nsp=300) 
      common/brintw/riw,riwminus,fw(nnn),fwminus(nnn)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      riw=0.d0
      riwminus=0.d0
      do 1 j=1,nr
      x1=j*dr
      x2=x1*x1
      c=wf(j,ia)*wf(j,ib)*wf(j,ic)*wf(j,id)*dr/x2
      riw=riw+fw(j)*c
      riwminus=riwminus+fwminus(j)*c
    1 continue  
      return
      end

      subroutine rdint1_g(ia,ib,ic,id,nr,dr)
      implicit real*8 (a-h,o-z)
      parameter (nnn=180,nsp=300,ncf=680)
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      common/bfg/f0(nnn,2),g0(nnn,2)
      common/biso/it1
      common/bisor/tt(ncf,ncf),tt1(ncf,ncf),tt2(ncf,ncf),
     1tt3(ncf,ncf),tt4(ncf,ncf)
      common/bnri/nri
      dimension f0rdint(nnn),g0rdint(nnn)
      go to (4,5),it1
    4 do 2 j=1,nri
      f0rdint(j)=f0(j,1)
      g0rdint(j)=g0(j,1)
    2 continue
      go to 6
    5 do 3 j=1,nri
      f0rdint(j)=f0(j,2)
      g0rdint(j)=g0(j,2)
    3 continue
    6 continue      
      ri0=0.d0
      ris=0.d0
      rx1=0.d0
      xl=ll(ia)*(ll(ia)+1)+ll(ib)*(ll(ib)+1)+ll(ic)*(ll(ic)+1)
      xl=xl+ll(id)*(ll(id)+1)
      do 1 j=1,nr
      x1=j*dr
      x2=x1*x1
      c=wf(j,ia)*wf(j,ib)*wf(j,ic)*wf(j,id)*dr/x2
      ri0=ri0+f0rdint(j)*c
      ris=ris+g0rdint(j)*c
      c1=-xl*c/x2+2.d0*dr*((dwf(j,ia)*wf(j,ib)+wf(j,ia)*dwf(j,ib))
     1*wf(j,ic)*wf(j,id)+(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id))*
     2wf(j,ia)*wf(j,ib))/(x1*x2)
      c1=c1-2.d0*dr*(dwf(j,ia)*(dwf(j,ib)*wf(j,ic)*wf(j,id)+wf(j,
     1ib)*(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id)))+wf(j,ia)*(dwf(
     2j,ib)*(dwf(j,ic)*wf(j,id)+wf(j,ic)*dwf(j,id))+wf(j,ib)*dwf(
     3j,ic)*dwf(j,id)))/x2
      rx1=rx1+c1
    1 continue  
      return
      end

      subroutine rdint2_g(ia,ib,ic,id,nr,dr)
      implicit real*8 (a-h,o-z)
      parameter (nnn=180,nsp=300) 
      dimension ua(2),ub(2)
      common/brint/ri0,ris,rx1,rx2(2,2)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      common/bfg/f0(nnn,2),g0(nnn,2)
      do 2 k1=1,2
      do 6 k2=1,2
      rx2(k1,k2)=0.d0
    6 continue
    2 continue
      do 1 j=1,nr
      x1=j*dr
      x2=x1*x1
      do 3 k=1,2
      sk=3.d0-2.d0*k
      ua(k)=dwf(j,ia)+sk*(ll(ia)+k-1)*wf(j,ia)/x1
      ub(k)=dwf(j,ib)+sk*(ll(ib)+k-1)*wf(j,ib)/x1
    3 continue
      do 4 k1=1,2
      do 5 k2=1,2
      c23=ua(k1)*ub(k2)*wf(j,ic)*wf(j,id)
      rx2(k1,k2)=rx2(k1,k2)+c23*dr/x2
    5 continue
    4 continue
    1 continue  
      return
      end

      subroutine selfph (it,ip1,ih1,ip2,ih2,omg,ga,nvs,nvsminus)
      implicit real*8(a-h,o-z)
      character*80 file_doo
c     real*4 sp1,sh1,sjtot,sh2,sp2,slfo,seij 
      parameter (nnn=180,nmx=12,nsp=300,ncf=680,ncf2=2*ncf,nft=3) 
      parameter (npa=270,nho=38,ncht=7)
      complex*16 den1,den1minus,ctr1,slf1,ctr1minus,slf1minus
      complex*16 den2,den2minus,ctr2,slf2,ctr2minus,slf2minus
      complex*16 den3,den3minus,ctr3,slf3,ctr3minus,slf3minus
      complex*16 den4,den4minus,ctr4,slf4,ctr4minus,slf4minus
      complex*16 slfs,slfminuss
      complex*16 nvs,nvsminus
      common/bfile2/file_doo
      common/blecp1/t0,t1,t2,t3,t13,alfe,x0,x1,x2,x3,anucl,znucl
      common/new/itd,ireal,idw
      common/bptl/ipartl,nchanl(2),ichanl(ncht),nchanlt,
     1nos1,npo1,nossp,nossn,nholes
      common/bptlr/qtrnsf,spect(ncht),erin(ncht)
      common/blecp2/dt(nnn),taut(nnn),dt1(nnn),dt2(nnn),dp(nnn),dn(nnn) 
      common/bfoti/nfo(nft),lfo(nft),jfo(nft),isfo(nft),lifo(nft),
     1itrans(nft),ifon,ntotfo,iformat
      common/bfotr/omgfo(nft),alphasq(nft) 
      common/bdens/rintpp(npa,npa,nft,2),rinthh(nho,nho,nft,2),rintph
     1(npa,nho,nft,2)
      common/bself/rslf12,aslf12,rslf34,aslf34
      common/bvari/iprint,nctr1,nctr2,nctr3,nctr4
      common/bvari12/nctrp1,nctrn1,nctr01,nctrp2,nctrn2,nctr02
      common/bvari34/nctrp3,nctrn3,nctr03,nctrp4,nctrn4,nctr04      
      common/btt/tiso
      common/bsh/ishort2,idoo
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bqwf/nn(nsp),ll(nsp),lj(nsp),lt(nsp),ehf(nsp),xnor(nsp)
      common/bwu1/nconf,ntot,nconf1,nconf2,nconf3,nconft
      common/bwu1r/sqbel(ncf),charge(ncf)
      common/bwu2/ih(ncf),ip(ncf),ecf(ncf),nig(ncf,2),nqqx(50,2),
     1nqqy(50,2),ngr(2),irpab,jtot,itot,ipi,ichex1,irad1
      if(idw.eq.1)open(unit=14,status='new',file='doorway.dat') 
      csq=1.d0
      if(ichex.eq.0)go to 7
      if(it.eq.2)go to 7 
ccc   Isospin correction for charge-exchange modes 
      tiso1=tiso+1.d0
      if(irad.eq.0.and.jtot.eq.0.and.ipi.eq.1)then
      csq=1.d0/tiso1   
ccc   (IAR)
      else
      csq=(2.d0*tiso - 1.d0)/(2.d0*tiso + 1.d0)                 
ccc   (GTR,IMR)
      end if 
    7 n0=1
      ljp1=lj(ip1)
      ljp2=lj(ip2)
      ljh1=lj(ih1)
      ljh2=lj(ih2)
      llp1=ll(ip1)
      llp2=ll(ip2)
      llh1=ll(ih1)
      llh2=ll(ih2)
      sp1=float(ljp1)/2.0
      sp2=float(ljp2)/2.0
      sh1=float(ljh1)/2.0
      sh2=float(ljh2)/2.0
      sqjp1=float(ljp1+1)
      sqjh1=float(ljh1+1)
      sqjp2=float(ljp2+1)
      sqjh2=float(ljh2+1)
      sjtot=float(jtot)
ccc   ----- W1, corresponding to the effective mass of the particle --
      slf1=(0.d0,0.d0)
      slf1minus=(0.d0,0.d0)
      if (ih1.ne.ih2) go to 2
      if (ljp1.ne.ljp2) go to 2
      k1=ip1-nos
      k2=ip2-nos
      do 12 ia=nos1,ntot
      if (lt(ia).ne.lt(ip1)) go to 12
      if (lt(ia).ne.lt(ip2)) go to 12
      ljia=lj(ia)
      llia=ll(ia)
      ka=ia-nos
      do 14 n=1,ntotfo
      llfo=lfo(n)
      ljfo=jfo(n)
      isf=isfo(n)
      ltfo=lifo(n)+1
      itr=itrans(n) 
      if (ljfo.lt.iabs((ljp1-ljia)/2)) go to 11
      if (ljfo.gt.(ljp1+ljia)/2) go to 11
      if (ljfo.lt.iabs((ljp2-ljia)/2)) go to 11
      if (ljfo.gt.(ljp2+ljia)/2) go to 11
      if (mod(llfo+llp1+llia,2).eq.1) go to 11
      if (mod(llfo+llp2+llia,2).eq.1) go to 11
      if (ifon.eq.1) then
      ri1=rintpp(k2,ka,itr,ltfo)
      ri3=rintpp(k1,ka,itr,ltfo)
      end if
      if (ifon.eq.2) then
      ri1=rintpp(k2,ka,n,1)
      ri3=rintpp(k1,ka,n,1)
      end if
      den1r=omg-(omgfo(n)+ehf(ia)-ehf(ih1))
      den1=cmplx(den1r,0.5d0*ga)
      den1rminus=den1r-2*omg
      den1minus=cmplx(den1rminus,0.5d0*ga)
      if(isf.eq.0)sph=yl(llp1,ljp1,llia,ljia,llfo)**2
      if(isf.eq.1)sph=tjl(llp1,ljp1,llia,ljia,ljfo,llfo)**2
      ctr1=alphasq(n)*ri1*ri3*(sph/sqjp1)/den1
      actr1=dimag(ctr1)
      if (actr1.gt.1.0d-30) nctrp1=nctrp1+1
      if (actr1.lt.-1.0d-30) nctrn1=nctrn1+1
      if (iprint.ne.1) go to 16
      write(2,*) 'n,ip1,ip2,ih1,ih2,ia=',n,ip1,ip2,ih1,ih2,ia
      write(2,*) 'ri1=',ri1
      write(2,*) 'ri3=',ri3
      write(2,*) 'den1r=',den1r
      write(2,*) 'den1=',den1
      write(2,*) 'sph=',sph
      write(2,*) 'sqjp1=',sqjp1
      write(2,*) 'ctr1=',ctr1
   16 ctr1minus=alphasq(n)*ri1*ri3*(sph/sqjp1)/den1minus
      if (iprint.ne.1) go to 17
      write(2,*) 'den1rminus=',den1rminus
      write(2,*) 'den1minus=',den1minus
      write(2,*) 'ctr1minus=',ctr1minus
   17 slf1=slf1+ctr1
      slf1minus=slf1minus+ctr1minus
      edoo=omgfo(n)+ehf(ia)-ehf(ih1)
      if(idw.eq.1.and.it.eq.2)then
      idoo=idoo+1
      write(14,*) ia,n,ih1,edoo
      end if 
      go to 15
   11 nctr01=nctr01+1
   15 nctr1=nctr1+1
   14 continue
   12 continue
      slf1=slf1*csq
      slf1minus=slf1minus*csq
      if (iprint.ne.1) go to 2
      write(2,*) 'nctrp1,nctrn1,nctr01,nctr1=',nctrp1,nctrn1,nctr01,nctr1
      write(2,*) 'slf1=',slf1
      write(2,*) 'slf1minus=',slf1minus
ccc   ----- W2, corresponding to the effective mass of the hole ------
    2 slf2=(0.d0,0.d0) 
      slf2minus=(0.d0,0.d0)
      if (ip1.ne.ip2) go to 3
      if (ljh1.ne.ljh2) go to 3
      kh1=ih1
      kh2=ih2
      do 22 ib=n0,nos
c     if (nholes.eq.1.and.ib.lt.nossp) go to 22
c     if (nholes.eq.1.and.ib.ge.npo1.and.ib.lt.nossn) go to 22
      if (lt(ib).ne.lt(ih1)) go to 22
      if (lt(ib).ne.lt(ih2)) go to 22
      ljib=lj(ib)
      llib=ll(ib)
      kb=ib
      do 24 n=1,ntotfo
      llfo=lfo(n)
      ljfo=jfo(n)                                  
      isf=isfo(n)  
      ltfo=lifo(n)+1
      itr=itrans(n)
      if (ljfo.lt.iabs((ljh1-ljib)/2)) go to 21
      if (ljfo.gt.(ljh1+ljib)/2) go to 21
      if (ljfo.lt.iabs((ljh2-ljib)/2)) go to 21
      if (ljfo.gt.(ljh2+ljib)/2) go to 21
      if (mod(llfo+llh1+llib,2).eq.1) go to 21
      if (mod(llfo+llh2+llib,2).eq.1) go to 21
      if (ifon.eq.1) then
      ri2=rinthh(kh2,kb,itr,ltfo)
      ri4=rinthh(kh1,kb,itr,ltfo)
      end if
      if (ifon.eq.2) then
      ri2=rinthh(kh2,kb,n,1)
      ri4=rinthh(kh1,kb,n,1)
      end if
      den2r=omg-(omgfo(n)+ehf(ip1)-ehf(ib))
      den2=cmplx(den2r,0.5d0*ga)
      den2rminus=den2r-2*omg
      den2minus=cmplx(den2rminus,0.5d0*ga)
      if(isf.eq.0)sph=yl(llh1,ljh1,llib,ljib,llfo)**2
      if(isf.eq.1)sph=tjl(llh1,ljh1,llib,ljib,ljfo,llfo)**2
      ctr2=alphasq(n)*ri2*ri4*(sph/sqjh1)/den2
      actr2=dimag(ctr2)
      if (actr2.gt.1.0d-30) nctrp2=nctrp2+1
      if (actr2.lt.-1.0d-30) nctrn2=nctrn2+1
      if (iprint.ne.1) go to 26
      write(2,*) 'n,ip1,ip2,ih1,ih2,ib=',n,ip1,ip2,ih1,ih2,ib
      write(2,*) 'ri2=',ri2
      write(2,*) 'ri4=',ri4
      write(2,*) 'den2r=',den2r
      write(2,*) 'den2=',den2
      write(2,*) 'sph=',sph
      write(2,*) 'sqjh1=',sqjh1
      write(2,*) 'ctr2=',ctr2
   26 ctr2minus=alphasq(n)*ri2*ri4*(sph/sqjh1)/den2minus
      if (iprint.ne.1) go to 27
      write(2,*) 'den2rminus=',den2rminus
      write(2,*) 'den2minus=',den2minus
      write(2,*) 'ctr2minus=',ctr2minus
   27 slf2=slf2+ctr2
      slf2minus=slf2minus+ctr2minus
      go to 25
   21 nctr02=nctr02+1
   25 nctr2=nctr2+1
   24 continue
   22 continue
      slf2=slf2*csq
      slf2minus=slf2minus*csq
      if (iprint.ne.1) go to 3
      write(2,*) 'nctrp2,nctrn2,nctr02,nctr2=',nctrp2,nctrn2,nctr02,nctr2
      write(2,*) 'slf2=',slf2
      write(2,*) 'slf2minus=',slf2minus
    3 rslf12=dreal(slf1+slf2)
      aslf12=dimag(slf1+slf2)
ccc   ----- W3,W4 corresponding to vertex corrections ----------------
      slf3=(0.d0,0.d0)
      slf4=(0.d0,0.d0)
      slf3minus=(0.d0,0.d0)
      slf4minus=(0.d0,0.d0)
      kp1=ip1-nos
      kp2=ip2-nos
      kh1=ih1
      kh2=ih2     
      do 34 n=1,ntotfo
      llfo=lfo(n)
      slfo=float(llfo)
      ljfo=jfo(n)
      sjfo=float(ljfo)                                  
      isf=isfo(n)  
      ltfo=lifo(n)+1
      itr=itrans(n)
      if (ljfo.lt.iabs((ljp1-ljp2)/2)) go to 31
      if (ljfo.gt.(ljp1+ljp2)/2) go to 31
      if (ljfo.lt.iabs((ljh1-ljh2)/2)) go to 31
      if (ljfo.gt.(ljh1+ljh2)/2) go to 31
      if (mod(llfo+llp1+llp2,2).eq.1) go to 31
      if (mod(llfo+llh1+llh2,2).eq.1) go to 31
      if (ifon.eq.1) then
      ri1=rintpp(kp1,kp2,itr,ltfo)
      ri4=rinthh(kh1,kh2,itr,ltfo)
      end if
      if (ifon.eq.2) then
      ri1=rintpp(kp1,kp2,n,1)
      ri4=rinthh(kh1,kh2,n,1)
      end if
      den3r=omg-(omgfo(n)+ehf(ip1)-ehf(ih2))
      den4r=omg-(omgfo(n)+ehf(ip2)-ehf(ih1))
      den3=cmplx(den3r,0.5d0*ga)
      den4=cmplx(den4r,0.5d0*ga)
      den3rminus=den3r-2*omg
      den4rminus=den4r-2*omg
      den3minus=cmplx(den3rminus,0.5d0*ga)
      den4minus=cmplx(den4rminus,0.5d0*ga)
      iphase=1+((ljp2+ljh2)/2)+jtot+ljfo
      sgn=(-1)**(mod(iphase,2))
      call sixj(sp1,sh1,sjtot,sh2,sp2,sjfo,seij)
      if(isf.eq.0)then 
       sph1=yl(llp1,ljp1,llp2,ljp2,llfo)
       sph2=yl(llh2,ljh2,llh1,ljh1,llfo)
      else 
       sph1=tjl(llp1,ljp1,llp2,ljp2,ljfo,llfo)
       sph2=tjl(llh2,ljh2,llh1,ljh1,ljfo,llfo)
      end if 
      ctr3=(sgn*alphasq(n)*ri1*ri4*sph1*sph2*dble(seij))/den3
      ctr4=(sgn*alphasq(n)*ri1*ri4*sph1*sph2*dble(seij))/den4
      actr3=dimag(ctr3)
      actr4=dimag(ctr4)
      if (actr3.gt.1.0d-30) nctrp3=nctrp3+1
      if (actr3.lt.-1.0d-30) nctrn3=nctrn3+1
      if (actr4.gt.1.0d-30) nctrp4=nctrp4+1
      if (actr4.lt.-1.0d-30) nctrn4=nctrn4+1
      if (iprint.ne.1) go to 36
      write(2,*) 'n,ip1,ip2,ih1,ih2=',n,ip1,ip2,ih1,ih2
      write(2,*) 'ri1=',ri1
      write(2,*) 'ri4=',ri4
      write(2,*) 'den3r=',den3r
      write(2,*) 'den3=',den3
      write(2,*) 'den4r=',den4r
      write(2,*) 'den4=',den4
      write(2,*) 'sph1=',sph1
      write(2,*) 'sph2=',sph2
      write(2,*) 'seij=',seij
      write(2,*) 'ctr3=',ctr3
      write(2,*) 'ctr4=',ctr4
   36 ctr3minus=(sgn*alphasq(n)*ri1*ri4*sph1*sph2*dble(seij))/den3minus
      ctr4minus=(sgn*alphasq(n)*ri1*ri4*sph1*sph2*dble(seij))/den4minus
      if (iprint.ne.1) go to 37
      write(2,*) 'den3rminus=',den3rminus
      write(2,*) 'den3minus=',den3minus
      write(2,*) 'ctr3minus=',ctr3minus
      write(2,*) 'den4rminus=',den4rminus
      write(2,*) 'den4minus=',den4minus
      write(2,*) 'ctr4minus=',ctr4minus
   37 slf3=slf3+ctr3
      slf4=slf4+ctr4
      slf3minus=slf3minus+ctr3minus
      slf4minus=slf4minus+ctr4minus
      go to 35
   31 nctr03=nctr03+1
      nctr04=nctr04+1
   35 nctr3=nctr3+1
      nctr4=nctr4+1
   34 continue
      slf3=slf3*csq
      slf3minus=slf3minus*csq
      slf4=slf4*csq
      slf4minus=slf4minus*csq
      if (iprint.ne.1) go to 4
      write(2,*) 'nctrp3,nctrn3,nctr03,nctr3=',nctrp3,nctrn3,nctr03,nctr3
      write(2,*) 'slf3=',slf3
      write(2,*) 'slf3minus=',slf3minus
      write(2,*) 'nctrp4,nctrn4,nctr04,nctr4=',nctrp4,nctrn4,nctr04,nctr4
      write(2,*) 'slf4=',slf4
      write(2,*) 'slf4minus=',slf4minus
    4 continue
      rslf34=dreal(slf3+slf4)
      aslf34=dimag(slf3+slf4)
ccc   ----- Sum of the four terms ------------------------------------
    6 continue
      nvs=slf1+slf2+slf3+slf4 
      nvsminus=slf1minus+slf2minus+slf3minus+slf4minus 
      if(idw.eq.1)close(14) 
      return
      end

      subroutine spot(l,jd,lc,elab,zz)                                          
      implicit real*8 (a-h,o-z)
      parameter (nnn=180)  
      common/bpas/h,cof                                                         
      common/bpot/xmn(nnn),xmp(nnn),vn(nnn),vp(nnn),yn(nnn),yp(nnn),            
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)                               
      common/pot/u(nnn),rmass,eta,xka,xmsm(nnn),sig, 
     1sgne,dreg(nnn),dirg(nnn),c2d,s2d,sigma  
      common/bnri/nri
      common/bnrir/dr
      hb=rmass/0.04823d0                                                        
      sq=dsqrt(0.04823d0) 
      if(elab) 4,4,3                                                            
    3 continue                                                                  
      sgne=1.d0                                                                 
      eta=0.15805086d0*zz/dsqrt(elab)  
      if(lc.eq.0) eta=0.d0                                                      
      xka=sq*dsqrt(elab)             
      go to 5                                                                   
    4 eta=0.d0                                                                  
      sgne=-1.d0                                                                
      eminus=-elab                                                              
      xka=sq*dsqrt(eminus)           
    5 sgn2=sgne*xka*xka                                                         
      x1=float(l*(l+1))                                                         
      x2=float(jd*(jd+2)-4*l*(l+1)-3)/8.d0                                      
      ap=float(lc)                                                              
      an=1.d0-ap                                                                
      do 2 i=1,nri                                                              
      xx=h*h*float(i*i)                                                         
      y1=1.d0/(an*xmn(i)+ap*xmp(i))                                             
      y2=an*vn(i)+ap*vp(i)                                                      
      y3=an*vsn(i)+ap*vsp(i)                                                    
      y4=an*(yn(i)**2)+ap*(yp(i)**2)                                            
      y5=0.5d0*(an*wnd(i)+ap*wpd(i))                                            
      xmsm(i)=hb*(an/xmn(i)+ap/xmp(i))                                          
      y6=1.d0-xmsm(i)                                                           
      u(i)=x1/xx+y6*sgn2                                                        
    2 u(i)=u(i)+y1*(y2+x2*y3+ap*vc(i)-0.25d0*y1*y4+y5)                          
      return                                                                    
      end                                                                       

      subroutine transition
      implicit real*8(a-h,o-z)
      parameter (nnn=180,nmx=12,nsp=300,ncf=680,ncf2=2*ncf,nft=3) 
      parameter (npa=270,nho=38,ncht=7)
      character*80 file_phonons_td(nft) 
      common/bdens/rintpp(npa,npa,nft,2),rinthh(nho,nho,nft,2),rintph
     1(npa,nho,nft,2)
      common/bfg/f0(nnn,2),g0(nnn,2)
      common/bfile3/file_phonons_td 
      common/bfoti/nfo(nft),lfo(nft),jfo(nft),isfo(nft),lifo(nft),
     1itrans(nft),ifon,ntotfo,iformat
      common/bfotr/omgfo(nft),alphasq(nft) 
      common/bnri/nri
      common/bnrir/dr
      common/bptl/ipartl,nchanl(2),ichanl(ncht),nchanlt,
     1nos1,npo1,nossp,nossn,nholes
      common/bptlr/qtrnsf,spect(ncht),erin(ncht)
      common/bwf/wf(nnn,nsp),dwf(nnn,nsp)
      common/bwu1/nconf,ntot,nconf1,nconf2,nconf3,nconft
      common/bwu1r/sqbel(ncf),charge(ncf)
      dimension rho(nnn)
  248 format(10(1x,e12.5))

      do 10 i=1,ntotfo
      itau=lifo(i)+1
      isf=isfo(i) 
      j=50+i
      open(unit=j,status='old',file=file_phonons_td(i))
      if(iformat.eq.1)then
      do k=1,nri
      read(50+i,*) rrr,rho(k)
      end do
      else
      read(50+i,248) (rho(k),k=1,nri)
      end if
      close(j)
      do 1 ip1=nos1,ntot
      do 2 ip2=nos1,ntot
      ripp=0.
      do 3 k=1,nri
      vph=f0(k,itau)
      if(isf.eq.1)vph=g0(k,itau)
      cpp=wf(k,ip1)*wf(k,ip2)*vph*rho(k)*dr
      ripp=ripp+cpp
    3 continue
      k1=ip1-nos
      k2=ip2-nos
      rintpp(k1,k2,i,1)=ripp
   30 format(5x,i3,i3,i3,i3,i3,3x,e12.5)
    2 continue
    1 continue
      do 4 ih1=1,nos
      do 5 ih2=1,nos
      rihh=0.
      do 6 k=1,nri
      vph=f0(k,itau) 
      if(isf.eq.1)vph=g0(k,itau)
      chh=wf(k,ih1)*wf(k,ih2)*vph*rho(k)*dr
      rihh=rihh+chh
    6 continue
      k1=ih1
      k2=ih2
      rinthh(k1,k2,i,1)=rihh
    5 continue
    4 continue
      do 7 ip1=nos1,ntot
      do 8 ih1=1,nos
      riph=0.
      do 9 k=1,nri
      vph=f0(k,itau) 
      if(isf.eq.1)vph=g0(k,itau)
      cph=wf(k,ip1)*wf(k,ih1)*vph*rho(k)*dr
      riph=riph+cph
    9 continue
      k1=ip1-nos
      k2=ih1
      rintph(k1,k2,i,1)=riph
    8 continue
    7 continue
   10 continue
      return
      end

      subroutine wsax                                                           
      implicit real*8 (a-h,o-z)
      parameter (nnn=180)  
      common/bnri/nri
      common/bpot/xmn(nnn),xmp(nnn),vn(nnn),vp(nnn),yn(nnn),yp(nnn),            
     1vsn(nnn),vsp(nnn),vc(nnn),wnd(nnn),wpd(nnn)                               
      common/bpas/h,cof                                                         
      common/blecp1/t0,t1,t2,t3,t13,alfe,x0,x1,x2,x3,anucl,znucl
      common/blecp2/dt(nnn),taut(nnn),dt1(nnn),dt2(nnn),dp(nnn),dn(nnn)         
      dimension dal(3)                                                  
  100 format(2x,'Woods-Saxon parameters v0=',e12.5,5x,'r0=',                    
     1e12.5,5x,'a0=',e12.5)                                                     
  102 format(2x,'v0n,r0n,a0n,vson,rson,ason',3x,6e14.5)       
  103 format(2x,'v0p,r0p,a0p,vsop,rsop,asop',3x,6e14.5)
  104 format(2x,'coulomb radius',3x,e14.5)
  208 format(6e12.5)                                                            
  209 format(15i3)
      read *,t0,t1,t2,t3,t13,x0        
      read *,alfe,x3,x1,x2
      read *,anucl,znucl,h,(dal(i),i=1,3)
      read *,v0n,r0n,a0n,vson,rson,ason
      read *,v0p,r0p,a0p,vsop,rsop,asop
      read *,rcb
      write(2,102) v0n,r0n,a0n,vson,rson,ason
      write(2,103) v0p,r0p,a0p,vsop,rsop,asop
      write(2,104) rcb
      rcb3=rcb**3
      rcb2=rcb**2
      e2=(znucl-1.d0)*197.3d0/137.d0
      do 101 j=1,nri                                                       
      x=h*j                                                                     
      xmn(j)=(anucl+1.)/(anucl*0.04819)                                         
      xmp(j)=xmn(j)                                                             
      yn(j)=0.                                                                  
      yp(j)=0.
      x1=h*(j-1) 
      x2=h*(j+1)
      wnd(j)=0.d0
      wpd(j)=0.d0
      z=(x-r0n)/a0n      
      vn(j)=v0n/(1.d0+dexp(z))        
      z=(x-r0p)/a0p
      vp(j)=v0p/(1.d0+dexp(z))
      z1=(x1-rson)/ason
      z2=(x2-rson)/ason
      z=(1.d0/(1.d0+dexp(z2))-1.d0/(1.d0+dexp(z1)))/(2.d0*h)    
      vsn(j)=-vson*z/x
      z1=(x1-rsop)/asop
      z2=(x2-rsop)/asop
      z=(1.d0/(1.d0+dexp(z2))-1.d0/(1.d0+dexp(z1)))/(2.d0*h)
      vsp(j)=-vsop*z/x
      if(x.gt.rcb) go to 105
      vc(j)=(3.d0*rcb2-x*x)*e2/(2.d0*rcb3)
      go to 101
  105 vc(j)=e2/x               
  101 continue                                              
      return                                                                    
      end                                                                       

      function ws(energy,it)
      implicit real*8(a-h,o-z)
      it1=it+1
      go to (1,2), it1
    1 betas=6.8d0
      ep=-2.39d0
      rs2=23.5225d0
      mv4=1570410.d0
      ef=-5.65d0
      go to 3
    2 betas=9.35d0
      ep=-1.8d0
      rs2=40.0689d0
      mv4=13046616.d0
      ef=-5.9d0
    3 if(energy.ge.ef) ene=energy
      if(energy.lt.ef) then
      ene1=ef-energy
      ene =ef+ene1
      end if
      if(ene.lt.ef)then
      write(2,*) 'error in the simmetry of w'
      stop
      end if
      e2=(ene-ep)**2
      e4=e2**2
      ws=-betas*((e2/(e2+rs2))-(e4/(e4+mv4)))
      if(ene.lt.ep) ws=0.d0      
      return
      end      

      function wv(energy,it)
      implicit real*8(a-h,o-z)
      it1=it+1
      go to (1,2), it1
    1 gammav=8.18d0
      ep=-2.39d0
      mv4=1570410.d0
      ef=-5.65d0
      go to 3
    2 gammav=13.26d0
      ep=-1.8d0
      mv4=13046616.d0
      ef=-5.9d0
    3 if(energy.ge.ef) ene=energy
      if(energy.lt.ef) then
      ene1=ef-energy
      ene =ef+ene1
      end if
      if(ene.lt.ef)then
      write(2,*) 'error in the simmetry of w'
      stop
      end if
      e2=(ene-ep)**2
      e4=e2**2
      wv=-gammav*(e4/(e4+mv4))
      if(ene.lt.ep) wv=0.d0
      return
      end      

C-------------------------------------------------------------------
C ROUTINE TO SOLVE COMPLEX SYSTEM OF LINEAR EQUATIONS BY TRANFORMING
C TO A REAL SYSTEM SOLVED WITH THE ROUTINES OF NUMERICAL RECIPES.
C-------------------------------------------------------------------
C234567
      SUBROUTINE COMLIN(CA,N,NP,CB,CX)
      PARAMETER (NMAX=150,NP2=2*NMAX)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      DIMENSION CA(NP,NP),CB(NP),CX(NP)
      DIMENSION A(NP2,NP2),B(NP2),X(NP2),ALUD(NP2,NP2)
      DIMENSION INDX(NP2)

      IF (NP.GT.NMAX) stop  'Increase NMAX in COMPACK.'
      N2=2*N

      DO 10 I=1,N
      DO 20 J=1,N
      AR=DREAL(CA(I,J))
      AI=DIMAG(CA(I,J))
      ALUD(I,J)=AR
      ALUD(I,J+N)=-AI
      ALUD(I+N,J)=AI
      ALUD(I+N,J+N)=AR
 20   CONTINUE
      BR=DREAL(CB(I))
      BI=DIMAG(CB(I))
      X(I)=BR
      X(I+N)=BI
 10   CONTINUE
C MAKE A COPY OF A AND B FOR ITERATIVE IMPROVEMENT
      DO 11 I=1,N2
      DO 12 J=1,N2
      A(I,J)=ALUD(I,J)  
 12   CONTINUE
      B(I)=X(I)
 11   CONTINUE   
C CALL TO NUMERICAL RECIPES' ROUTINES
      CALL LUDCMP(ALUD,N2,NP2,INDX,D)
      CALL LUBKSB(ALUD,N2,NP2,INDX,X)
      CALL MPROVE(A,ALUD,N2,NP2,INDX,B,X)

      DO 30 I=1,N
      CX(I)=DCMPLX(X(I),X(I+N))
 30   CONTINUE

      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=150,TINY=1.D-30)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
 11     CONTINUE
        IF (AAMAX.EQ.0.D0) stop  'Singular matrix.'
        VV(I)=1.D0/AAMAX
 12   CONTINUE
      DO 19 J=1,N
        DO 14 I=1,J-1
          SUM=A(I,J)
          DO 13 K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
 13       CONTINUE
          A(I,J)=SUM
 14     CONTINUE
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
           DO 15 K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
 15        CONTINUE
           A(I,J)=SUM
           DUM=VV(I)*ABS(SUM)
           IF (DUM.GE.AAMAX) THEN
             IMAX=I
             AAMAX=DUM
           ENDIF
 16        CONTINUE
           IF (J.NE.IMAX) THEN
           DO 17 K=1,N
             DUM=A(IMAX,K)
             A(IMAX,K)=A(J,K)
             A(J,K)=DUM
 17        CONTINUE
           D=-D
           VV(IMAX)=VV(J)
          ENDIF
      INDX(J)=IMAX
      IF (A(J,J).EQ.0.D0) A(J,J)=TINY
      IF (J.NE.N) THEN
         DUM=1.D0/A(J,J)
         DO 18 I=J+1,N
           A(I,J)=A(I,J)*DUM
 18      CONTINUE
      ENDIF
 19   CONTINUE
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0) THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
 11         CONTINUE
         ELSE IF (SUM.NE.0) THEN
            II=I
         ENDIF
         B(I)=SUM
 12   CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
 13      CONTINUE
         B(I)=SUM/A(I,I)
 14   CONTINUE
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=150)
      DIMENSION A(NP,NP),ALUD(NP,NP),INDX(N),B(N),X(N),R(NMAX)
      DO 12 I=1,N
        SDP=-B(I)
        DO 11 J=1,N
          SDP=SDP+DBLE(A(I,J))*DBLE(X(J))
 11     CONTINUE
        R(I)=SDP
 12   CONTINUE
      CALL LUBKSB(ALUD,N,NP,INDX,R)
      DO 13 I=1,N
        X(I)=X(I)-R(I)
 13   CONTINUE
      RETURN
      END

      SUBROUTINE DIAGMA(NCONF,IRPA)
ccc   nouveau diagma du 12 Fevrier 1991
      implicit real*8 (a-h,o-z)
      parameter (ncf=680,ncf2=2*ncf)
ccc   method of Ullah-Rowe, Nucl. Phys. 163, p.257
      common/bdiam/aa(ncf,ncf),bb(ncf,ncf),evr(ncf),vecr(ncf2,ncf)  
      common/hsqr1/t1
      common/hsqr2/t2
      common/hsqr3/t3 
c     common/hsqr1,t1,hsqr2,t2,hsqr3,t3 
      dimension x(ncf,ncf),v(ncf),t1(ncf),t2(ncf),t3(ncf),e(ncf)
      dimension cec(ncf),tep(ncf,ncf),alp(ncf),valp(ncf),vec(ncf)
      dimension nord(ncf)
      ndim=nconf
      n=ndim
      eps=(1.0d-015)*n*n
      go to (401,400),irpa
  400 continue
      do1 i=1,n
      do1 j=1,i
      aa(i,j)=aa(i,j)+bb(i,j)
      bb(i,j)=aa(i,j)-2.d0*bb(i,j)
      aa(j,i)=aa(i,j)
      bb(j,i)=bb(i,j)
      x(i,j)=bb(i,j)
      x(j,i)=x(i,j)
    1 continue
      call diasym(x,ncf,n,v,e,vec,nord)
      do2 i=1,n
      cec(i)=v(i)
      do2 j=1,n
      tep(i,j)=x(i,j)
    2 continue
      do1000 i=1,n
      do1001 j=1,n
      xnew=0.d0
      do1002 k=1,n
      xnew=xnew+tep(k,i)*aa(k,j)
 1002 continue
      x(i,j)=xnew
 1001 continue
 1000 continue
      do1003 i=1,n
      do1004 j=1,i
      bbij=0.d0
      do1005 k=1,n
      bbij=bbij+x(i,k)*tep(k,j)
 1005 continue
      bb(i,j)=bbij
      bb(j,i)=bbij
 1004 continue
 1003 continue
      do5 i=1,n
      do5 j=1,i
      arpms=dsqrt(cec(i))*bb(i,j)*dsqrt(cec(j))
      x(i,j)=arpms
      x(j,i)=x(i,j)
    5 continue
      call diasym(x,ncf,n,v,e,vec,nord) 
      do6 i=1,n
      alp(i)=v(i)
      if (alp(i).gt.0.d0) go to 150
      write(2,151) i
  151 format(2x,'state number',i3,'is imaginary')
      go to 6
  150 valp(i)=dsqrt(alp(i))
    6 continue
      do60 i=1,n
      do61 j=1,n
      x(i,j)=dsqrt(valp(j))*x(i,j)
   61 continue
   60 continue
      do70 i=1,n
      do7 j=1,n
      aa(i,j)=0.d0
      bb(i,j)=0.d0
      do8 k=1,n
      yy1=0.5d0*tep(i,k)*dsqrt(cec(k))*x(k,j)
      yy2=(0.5d0*tep(i,k)/dsqrt(cec(k)))*x(k,j)
      aa(i,j)=aa(i,j)+yy1
      bb(i,j)=bb(i,j)+yy2
    8 continue
      aa(i,j)=aa(i,j)/valp(j)
      if(alp(i)) 25,25,24
   24 aa(i,j)=aa(i,j)+bb(i,j)
      bb(i,j)=aa(i,j)-2.d0*bb(i,j)
   25 continue
    7 continue
   70 continue
      nd2=2*n
      do130 j=1,n
      evr(j)=valp(j)
      do131 i=1,nd2
      if(i.gt.n) go to 140
      ev=aa(i,j)
      go to 141
  140 ev=bb(i-n,j)
  141 vecr(i,j)=ev
  131 continue
  130 continue
      return
  401 do402 i=1,n
      do403 j=1,n
      x(i,j)=aa(i,j)
  403 continue
  402 continue
      call diasym(x,ncf,n,v,e,vec,nord) 
      do404 j=1,n
      evr(j)=v(j)
      do405 i=1,n
  405 vecr(i,j)=x(i,j)
  404 continue
      return
      end

      subroutine ccg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end

      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end

      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end

      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end

      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end

      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision s,tr,ti,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end

      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end

      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      function radho(n1,l1,n2,l2,l)
      implicit double precision(a-h,o-z)
      common/bafa0/afa0
      radho=0.d0
      if(n1.eq.n2)then
      ldiff=iabs(l1-l2)
      if(ldiff.ne.1)go to 1
      lmin=min(l1,l2)
      radho=dsqrt(dfloat(n1+lmin)+0.5)/dsqrt(afa0)
      else 
      ndiff=iabs(n1-n2)
      ldiff=iabs(l1-l2)
      if(ldiff.ne.1.or.ndiff.ne.1)go to 1
      if(n1.gt.n2.and.l1.gt.l2)go to 1
      if(n2.gt.n1.and.l2.gt.l1)go to 1
      nmin=min(n1,n2)
      lmin=min(l1,l2)
      radho=-dsqrt(dfloat(nmin))/dsqrt(afa0)
      end if
    1 return
      end 

      function radho2(n1,l1,n2,l2,l)
      implicit double precision(a-h,o-z)
      common/bafa0/afa0
      common/bnri/nri
      common/bnrir/dr
      b=1/dsqrt(afa0)
      radho2=0.d0
      do j=1,nri
      x=float(j)*dr
      radho2=radho2+dr*hwf(n1-1,l1,b,x)*
     1hwf(n2-1,l2,b,x)*x**l
      return
      end do
      end
