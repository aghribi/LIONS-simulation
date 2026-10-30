    program lions

	use lion_m
	implicit none
	integer nreste,npi,ier,k,icalc

    write(6,*)'----------------------------------------------------'
    write(6,*)'LIONS       : Fortran 90 portable'
	write(6,*)'auteur      : Patrick Bertrand - GANIL'
	write(6,*)'version  1  : decembre 1998...decembre 2004'
	write(6,*)'version  2  : mars 2017, modif seed, swap etc...'
	write(6,*)'              voir les modifs avec le mot clé : 2017'
    write(6,*)'----------------------------------------------------'
    write(6,*)' '
    write(6,*)'modif swap de mars 2017 :'
    write(6,*)'on swap les octets ? (y/n) (y si anciens fichiers data unformatted...)'
    read(5,*) repswap

	idebug	= 0
	call init

	write(6,*)'--------------------'
	write(6,*)'test trieur vertical :'
	write(6,*)'--------------------'
	write(6,*)'rmin_ver rmax_ver teta1_ver teta2_ver :'
    write(6,*)'1000 2000 0 0 si pas de tri vertical ...'
	read(5,*)rmin_ver,rmax_ver,teta1_ver,teta2_ver
	teta1_ver = teta1_ver*pi/180._double
	teta2_ver = teta2_ver*pi/180._double
	write(6,*)teta1_ver,teta2_ver
	write(6,*)'pot_ver,h_ver,phi_ver anuz coupnuz dpotspot :'
    write(6,*)'0 0 0 0 0 0 si pas de tri vertical...'
	read(5,*)pot_ver,h_ver,phi_ver,anuz,coupnuz,dpotspot
	phi_ver = phi_ver*pi/180._double
	write(6,*)'-------------------'

	call initgra
	ite = 0
	call boite
        call calcbe
	call info

	do ite=1,nbdt

	  call piege
	  call boite
          call calcbe
	  call pousse

	  t 	= ite*dt+t0


	  call info
	  if (numbid .ne. 0 .and. (mod(ite,ifreb).eq.0   &
     &	      .or.ite.eq.1))call bid

	end do

	if(numbid.ne.0)    call clobid(ier)
	if(ficz.ne.'no')   close(15)
	if(ficxy.ne.'no')  close(16)
	if(ficcor.ne.'no') close(17)
	if(ficmat.ne.'no') close(18)
	if(ficdc .ne.'no') close(19)
	if(nbrpie.ne.0)then
	  nreste	= count(ok)
	  write(6,*)'nombre de particules non piegees : ',nreste
	  write(6,*)'écriture fichier des numero particules piegees ? (o/n)'
	  read(5,*)reponse
	  if(reponse.eq.'o')then
	    write(6,*)'choisissez nom du fichier : '
	    read(5,*)ficpie
	    open(60,file=ficpie)
	    npi = 0
	    do k=1,n
	       if(.not.ok(k))then
	         npi = npi + 1
	         write(60,6000)npi,k
6000	         format(1x,i5,1x,i5)
	       endif
	    enddo
	    close(60)
	  endif
	endif

	end

!	------------------------------------------------------------------

	subroutine init

	use lion_m
	implicit none
	character*80 title
	real(kind=double)tetdfid,tetdfod,tetdmid,tetdmod
	real(kind=double)w0mev,am0k,alpha,r2,rc,co,si,rr
	integer kg,iharm,i1,i2,m,ntetall
	integer nx1,ny1,nx2,ny2,nbr,i,j,k

	pi	= acos(-1._double)
	clum		=  299792458._double
	rac8		=  sqrt(8._double)

	write(6,*)'fichier_de_donnees de départ :'
	read(5,*)ficin
	open(20,file=ficin,status='old')

	read(20,*)corps
	read(20,*)npaspt
	read(20,*)r1th,ntours
	read(20,*)fhf,pothf,teta1,harm,dbsb
        read(20,*)dphid,dphig
        dphid   = dphid*pi/180._double
        dphig   = dphig*pi/180._double
        read(20,*)pothfd,pothfg
        pothfd  = -pothfd
        pothfg  = -pothfg
        read(20,*)ncoh(1)
        if(ncoh(1).ne.0)read(20,*)(coh(i,1),i=1,ncoh(1))
        if(ncoh(1).ne.0)read(20,*)phihar(1)
        read(20,*)ncoh(2)
        if(ncoh(2).ne.0)read(20,*)(coh(i,2),i=1,ncoh(2))
        if(ncoh(2).ne.0)read(20,*)phihar(2)

	yomega	= 0._double
	if(teta1.lt.0._double)then
	  write(6,*)'attention : gap en biais...'
	  write(6,*)'on remet teta1 positif'
	  write(6,*)'et on va lire le decalage...'
	  teta1	= -teta1
	  read(20,*)yomega
	  write(6,*)'decalage yomega : ',yomega
	endif
	read(20,*)zmin,zmax,rmax1
	read(20,*)rex,tet1ex,tet2ex,eex
	ndf = 0
	ndm = 0
	if(rex.lt.0._double)then
	  if(tet1ex.le.0._double)then
	    iteext = 0
	    write(6,*)'extraction : le paquet sera coupe'
	  else
	    iteext = tet1ex + .5_double
	    write(6,*)'extraction : le paquet sera entier'
	  endif
	  read(20,*)ndf
	  if(ndf.gt.0)then
	    do i=1,ndf
	      read(20,*)rdf(i),edf(i),xdfi(i),ydfi(i), &
     &                  xdfo(i),ydfo(i),tetdfid,tetdfod
	      tetdfi(i) = tetdfid*pi/180._double
	      tetdfo(i) = tetdfod*pi/180._double
	      uxdfi(i)  = cos(tetdfi(i))
	      uydfi(i)  = sin(tetdfi(i))
	      uxdfo(i)  = cos(tetdfo(i))
	      uydfo(i)  = sin(tetdfo(i))
	    enddo
	  endif
	  read(20,*)ndm
	  if(ndm.gt.0)then
	    do i=1,ndm
	      read(20,*)rdm(i),bdm(i),xdm(i),ydm(i), &
     &                  tetdmid,tetdmod
	      tetdmi(i) = tetdmid*pi/180._double
	      tetdmo(i) = tetdmod*pi/180._double
	      uxdmi(i)  = cos(tetdmi(i))
	      uydmi(i)  = sin(tetdmi(i))
	      uxdmo(i)  = cos(tetdmo(i))
	      uydmo(i)  = sin(tetdmo(i))
	    enddo
	  endif
	  if(tet2ex.gt.0._double)then
	    write(6,*)'grilles extraction...'
	    itegex	= eex + .5_double
	    read(20,*)ngex
	    do kg=1,ngex
	      read(20,*)tetagex(kg),xgex(kg),ygex(kg),xgexc(kg),ygexc(kg)
	      xgex(kg) = xgex(kg)/1000._double
	      ygex(kg) = ygex(kg)/1000._double
	      xgexc(kg) = xgexc(kg)/1000._double
	      ygexc(kg) = ygexc(kg)/1000._double
	      drgex(kg) = sqrt( (xgex(kg)-xgexc(kg))**2 + &
     &                          (ygex(kg)-ygexc(kg))**2 )
	      uxgex(kg) = (xgex(kg)-xgexc(kg))/drgex(kg)
	      uygex(kg) = (ygex(kg)-ygexc(kg))/drgex(kg)
	    enddo
	    read(20,*)ficdebgex
	    do kg=1,ngex
              if(kg.lt.10)then
                write(ficfingex,'(i1)')kg
              else
                write(ficfingex,'(i2)')kg
	      endif
	      ficgex(kg) = ficdebgex(1:len_trim(ficdebgex)) // '.' &
     &                     // ficfingex(1:len_trim(ficfingex))
	      open(70+kg,file=ficgex(kg))
	      write(70+kg,*)'teta x y xc yc : '
	      write(70+kg,*)tetagex(kg),xgex(kg),ygex(kg) &
     &                       ,xgexc(kg),ygexc(kg)
	    enddo
	  endif
	endif

	dt 		=  harm/fhf/npaspt
	omhf 		=  2._double*pi*fhf
	omr 		=  omhf/harm
	cdt     	=  clum*dt
	teta1 		=  teta1*pi/180._double
	tet1ex		=  tet1ex*pi/180._double
	tet2ex		=  tet2ex*pi/180._double
	iharm		=  int(harm + .5_double)
	heps		=  1._double
	if((iharm/2)*2.eq.iharm)heps = -1._double

	read(20,*)nbrfam
	if(nbrfam.gt.1)read(20,*)n
	index(1) 	=  1
	np 		=  0
	do k=1,nbrfam
	  read(20,*)nbrion(k),amat(k),a(k),qn(k)
	  if(k.eq.1)then
	    if(nbrfam.eq.1)n = nbrion(1)
	    call allo
	    sortant 	= .false.
	    indgex	= 0
	    ok		= .true.
	    nbrok	= n
	    nbrpie	= 0
	    iflpie 	= 1
	  endif
	  read(20,*)r0(k),teta0(k),z0(k),pr0(k),pz0(k),w0mev,phi0(k)
	  if(k.eq.1)w0mev1 = w0mev
	  teta0(k)	=  teta0(k)*pi/180._double
	  pr0(k)	=  pr0(k)*.001_double
	  pz0(k)	=  pz0(k)*.001_double
	  np		=  np + nbrion(k)
	  if(np.gt.n)then
	    write(6,*)'n np : ',n,np
	    write(6,*)'changer n ou np...'
	   stop
	  end if
	  index(k+1)    =  index(k) + nbrion(k)
	  i1		=  index(k)
	  i2		=  index(k) + nbrion(k) - 1
	  x(i1:i2)	=  r0(k)*cos(teta0(k))
	  y(i1:i2)	=  r0(k)*sin(teta0(k))
	  z(i1:i2)	=  z0(k)
	  w0(k)		=  w0mev*1.60217733e-13_double*a(k)
	  am0k		=  1.6605402e-27_double*(amat(k) -             &
     &                      .0005486_double*qn(k))
	  am0(i1:i2)	=  am0k
	  am0c(i1:i2)	=  am0(i1:i2)*clum
	  am0c2(i1:i2)	=  am0(i1:i2)*clum**2
	  am02c2(i1:i2)	=  am0c(i1:i2)**2
	  ga(i1:i2)     =  1._double + w0(k)/am0c2(i1:i2)
	  alpha		=  (pi/2._double + teta0(k)) - pr0(k)
	  p(i1:i2)	=  sqrt((ga(i1:i2)-1._double)*    &
     &                          (ga(i1:i2)+1._double))
	  px(i1:i2)	=  p(i1:i2)*cos(alpha)*cos(pz0(k))
	  py(i1:i2)	=  p(i1:i2)*sin(alpha)*cos(pz0(k))
	  pz(i1:i2)	=  p(i1:i2)*sin(pz0(k))
          q(i1:i2)      =  qn(k)*1.60217733e-19_double
	  qdts2m0(i1:i2)=  q(i1:i2)*dt/2._double/am0(i1:i2)
	  qdtsm0c(i1:i2)=  q(i1:i2)*dt/am0c(i1:i2)
	end do
	if((tet2ex.gt.0._double).and.(ngex.gt.0))then
	  do kg = 1,ngex
	    write(70+kg,*)'caracteristiques particules :'
	    write(70+kg,*)np,a(1),amat(1),qn(1)
	    write(70+kg,*)'caracteristiques HF :'
	    write(70+kg,*)npaspt,fhf,harm
	  enddo
	endif

	read(20,*)iflcor
    ficsig = "pas de fichier correl..." ! Patrick Bertrand : ajout mars 2017
	if(iflcor.ne.0)then
	  m		=  6
	  read(20,*)ficsig
	  read(20,*)espace
	  open(25,file=ficsig,status='old')
	  read(20,*)ficcor,ficmat,ficdc
	  if(ficcor.ne.'no')then
	    open(17,file=ficcor,status='new')
	    write(17,2117)
2117	    format('  tr     r      teta    dr    ds      dphi',   &
     &      '   epsr    dwsw        wmpa     dw          dwdphi',  &
     &      '       dz     epsz')
	  end if
	  if(ficmat.ne.'no')open(18,file=ficmat,status='new')
	  if(ficdc.ne.'no')then
	    open(19,file=ficdc,status='new')
	    write(19,2119)
2119	    format('  tr     r      teta   drdc  dsdc  dphidc',     &
     &      ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc')
	  end if
	end if

	if(np.ne.n)then
	  write(6,*)'n np : ',n,np
	  write(6,*)'changer n ou np...'
	  stop
	end if
	t0		=  phi0(1)*pi/180._double/omhf
	t		=  t0

	read(20,*)nbdt,ifreq
	read(20,*)ficz,ficxy
	read(20,*)ficb,ficv1,ficv2
        ficloa  = 'no'
        ficsto  = 'no'
        read(20,*,end=800)ficloa,ficsto
800     close(20)

	open(21,file=ficb,status='old')
	read(21,*)title
	read(21,*)rmin,rmax,nr,ntet
	if(ntet.ne.91)then
	  write(6,*)'lecture carte tosca.big ou carte mesuree ...'
	  read(21,*)((tbz(i,j),j=1,ntet),i=1,nr)
          read(21,*)((d2bz(i,j),j=1,ntet),i=1,nr)
          read(21,*)((dbzx(i,j),j=1,ntet),i=1,nr)
          read(21,*)((dbzy(i,j),j=1,ntet),i=1,nr)
          close(21)
	  htet  = 2._double*pi/(ntet-1._double)
	  ntetall = ntet
	  goto 900
	endif
	write(6,*)'carte tosca.small (1/8eme) ...'
	read(21,*)((tbz(i,92-j),j=1,ntet),i=1,nr)
	read(21,*)((d2bz(i,92-j),j=1,ntet),i=1,nr)
	read(21,*,end=444)((dbzxlu(i,j),j=1,ntet),i=1,nr)
	read(21,*)((dbzylu(i,j),j=1,ntet),i=1,nr)
	close(21)
	goto 445
444	write(6,*)'attention fichier magnetique incomplet...'
	write(6,*)'bx et by seront nuls...'
	write(6,*)'vous resterez dans le plan median...'
	dbzxlu(1:nr,1:ntet) = 0
	dbzylu(1:nr,1:ntet) = 0
	close(21)
445	continue
        tbz(1:nr,181:92:-1) = tbz(1:nr,1:90)
        tbz(1:nr,182:361)   = tbz(1:nr,2:181)
        tbz(1:nr,362:721)   = tbz(1:nr,2:361)
	d2bz(1:nr,181:92:-1)= d2bz(1:nr,1:90)
	d2bz(1:nr,182:361)  = d2bz(1:nr,2:181)
	d2bz(1:nr,362:721)  = d2bz(1:nr,2:361)
	r2		    = sqrt(2._double)/2._double
        dbzx(1:nr,91:1:-1)  = r2*(dbzxlu(1:nr,1:91) + &
     &                            dbzylu(1:nr,1:91))
        dbzy(1:nr,91:1:-1)  = r2*(dbzxlu(1:nr,1:91) - &
     &                            dbzylu(1:nr,1:91))
        dbzx(1:nr,181:92:-1)=  dbzy(1:nr,1:90)
        dbzy(1:nr,181:92:-1)=  dbzx(1:nr,1:90)
        dbzx(1:nr,182:361)  = -dbzx(1:nr,180:1:-1)
        dbzy(1:nr,182:361)  =  dbzy(1:nr,180:1:-1)
        dbzx(1:nr,362:721)  =  dbzx(1:nr,360:1:-1)
        dbzy(1:nr,362:721)  = -dbzy(1:nr,360:1:-1)
	htet   		= pi/4._double/(ntet-1._double)
	ntetall = 721
900	continue
	if( (ncoh(1).eq.0) .and. (ncoh(2).eq.0) )goto 910
	write(6,*)'carte mag. + rajout harmonique 1 et/ou 2...'
        call addh(rmin,rmax,nr,ntetall,tbz,d2bz,dbzx,dbzy,ncoh, &
     &            coh,phihar)
910	hr   		= (rmax-rmin)/(nr-1._double)
	hrm1		= 1._double/hr
	htetm1 		= 1._double/htet
	hrhtet		= hr*htet
	hrhtetm1	= 1._double/hrhtet
	tbz		=  tbz*(1._double + dbsb)
	d2bz		=  d2bz*(1._double + dbsb)

	if(ficv1.ne.'no')then

!	  Patrick Bertrand : modif swap mars 2017 (pour lire anciens fichiers sun, cray etc...
!     lié à problème d'ordre des octets sur processuer Intel...

      if (repswap.eq.'n')then
	    open(22,file=ficv1,form='unformatted',status='old')
      else if (repswap.eq.'y')then
        open(22,file=ficv1,form='unformatted',convert='BIG_ENDIAN',status='old')
      else
        write(6,*) 'stop : erreur de la réponse sur le swap des octets...'
        stop
      endif

	  read(22)xmin1,xmax1,ymin1,ymax1,nx1,ny1
	  xmin1		=  xmin1/1000._double
	  xmax1		=  xmax1/1000._double
	  ymin1		=  ymin1/1000._double
	  ymax1		=  ymax1/1000._double
	  do i=1,nx1
	    read(22)(pot1(i,j),j=1,ny1)
	  end do
	  do i=1,nx1
	    read(22,end=544)(d2pot1(i,j),j=1,ny1)
	  end do
	  close(22)
	else
	  write(6,*)'pas de fichier centre cyclo...'
	  write(6,*)'xmin1 xmax1 ymin1 ymax1 nx1 ny1 ?'
	  read(5,*)xmin1,xmax1,ymin1,ymax1,nx1,ny1
	  pot1 = 0._double
	  d2pot1 = 0._double
	endif
	goto 545
544	write(6,*)'attention fichier potentiel centre incomplet...'
	write(6,*)'ez sera nul...'
	write(6,*)'vous resterez dans le plan median...'
	d2pot1 = 0._double
545	hx1		=  (xmax1-xmin1)/(nx1-1._double)
	hy1		=  (ymax1-ymin1)/(ny1-1._double)
	hx1hy1		=  hx1*hy1
        hx1m1           = 1._double/hx1
        hy1m1           = 1._double/hy1
        hx1hy1m1        = 1._double/hx1hy1

	open(23,file=ficv2,status='old')
	read(23,*)xmin2,xmax2,ymin2,ymax2,nx2,ny2
	read(23,*)((pot2(i,j),j=1,ny2),i=1,nx2)
	read(23,*)((d2pot2(i,j),j=1,ny2),i=1,nx2)
	close(23)

	if(iflcor.ne.0)then
          call boite
          call calcb
          do k=1,nbrfam
            nbr           =  nbrion(k)
	    call remp(m,nbr,dxs)
            do i=1,6
              read(25,*)(cor(i,j),j=1,6)
            end do
            call choles(m,cor,rcor)
            trcor         =  transpose(rcor)
	    dxtmp(1:6,1:nbr) = 0._double
	    do i=1,6
	      do j=1,6
	        dxtmp(i,1:nbr) = dxtmp(i,1:nbr) + trcor(i,j)*   &
     &                           dxs(j,1:nbr)
	      end do
	    end do
            dxs(1:6,1:nbr)=  dxtmp(1:6,1:nbr)
            i1            =  index(k)
            i2            =  index(k) + nbrion(k) - 1
            rc            =  -am0c(k)*p(i1)/q(i1)/bz(i1)
	    if(r0(k).gt.rmax1)then
	      write(6,*)'prise en compte de sqrt(cos(phi))...'
              vec4        =  sqrt(cos(harm*dxs(5,1:nbr)/rc))
              dxs(1,1:nbr)=  dxs(1,1:nbr)*vec4(1:nbr) - &
     &                       rc*(1._double-vec4(1:nbr))
	    else
	      write(6,*)'pas de sqrt(cos(phi))...'
	    endif
            vec1(1:nbr)   =  cos(dxs(5,1:nbr)/rc)
            vec2(1:nbr)   =  sin(dxs(5,1:nbr)/rc)
            vec3(1:nbr)   =  sin(dxs(5,1:nbr)/rc/2._double)**2
            dxd(1,1:nbr)  =  dxs(1,1:nbr) - 2._double*   &
     &                       (rc + dxs(1,1:nbr))*vec3(1:nbr)
	    dxs(6,1:nbr)  =  sqrt((p(i1)+dxs(6,1:nbr))**2         &
     &	                     - dxs(2,1:nbr)**2 - dxs(4,1:nbr)**2) &
     &                       - p(i1)
            dxd(5,1:nbr)  =  (rc + dxs(1,1:nbr))*vec2(1:nbr)
            dxd(2,1:nbr)  =  dxs(2,1:nbr)*vec1(1:nbr)  -          &
     &                       dxs(6,1:nbr)*vec2(1:nbr)
            dxd(6,1:nbr)  =  dxs(2,1:nbr)*vec2(1:nbr)  +          &
     &                       dxs(6,1:nbr)*vec1(1:nbr)
            dxd(3:4,1:nbr)=  dxs(3:4,1:nbr)
            co      =  cos(teta0(k)-pr0(k))
            si      =  sin(teta0(k)-pr0(k))
	    rot	      =  0._double
            rot(1,1)  =  co
            rot(1,5)  =  -si
            rot(2,2)  =  co
            rot(2,6)  =  -si
            rot(3,3)  =  1._double
            rot(4,4)  =  1._double
            rot(5,1)  =  si
            rot(5,5)  =  co
            rot(6,2)  =  si
            rot(6,6)  =  co
	    dxtmp(1:6,1:nbr) = 0._double
	    do i=1,6
	      do j=1,6
	        dxtmp(i,1:nbr) = dxtmp(i,1:nbr) + rot(i,j)*   &
     &                           dxd(j,1:nbr)
	      end do
	    end do
            dxd(1:6,1:nbr)=  dxtmp(1:6,1:nbr)
            x(i1:i2)      =  x(i1:i2)  + dxd(1,1:nbr)
            y(i1:i2)      =  y(i1:i2)  + dxd(5,1:nbr)
	    vec(1:nbr)	  =  px(i1:i2)*vec1(1:nbr) -              &
     &                       py(i1:i2)*vec2(1:nbr)
            py(i1:i2)     =  px(i1:i2)*vec2(1:nbr) +              &
     &                       py(i1:i2)*vec1(1:nbr) + dxd(6,1:nbr)
            px(i1:i2)     =  vec(1:nbr) + dxd(2,1:nbr)
            if(espace(1:2).eq.'3d')then
              z(i1:i2)    =  z(i1:i2)  + dxd(3,1:nbr)
              pz(i1:i2)   =  pz(i1:i2) + dxd(4,1:nbr)
            end if
	    if(r0(k).gt.rmax1)then
              px(i1:i2) = px(i1:i2)*vec4(1:nbr)
              py(i1:i2) = py(i1:i2)*vec4(1:nbr)
              pz(i1:i2) = pz(i1:i2)*vec4(1:nbr)
	    endif
          end do
          close(25)

          if(ficloa.ne.'no')then
            write(6,*)'attention : chargement paquet d ions...'

!	        Patrick Bertrand : modif swap mars 2017

            if (repswap.eq.'n')then
              open(26,file=ficloa,form='unformatted',status='old')
            else if (repswap.eq.'y')then
              open(26,file=ficloa,form='unformatted',convert='SWAP',status='old')
            else
              write(6,*) 'stop : erreur a l open sur reponse swap octets...'
              stop
            endif
            read(26)(x(i),y(i),z(i),i=1,n)
            read(26)(px(i),py(i),pz(i),i=1,n)
            close(26)
          endif
	end if

	hx2		=  (xmax2-xmin2)/(nx2-1._double)
	hy2		=  (ymax2-ymin2)/(ny2-1._double)
	hx2hy2		=  hx2*hy2
        hx2m1           = 1._double/hx2
        hy2m1           = 1._double/hy2
        hx2hy2m1        = 1._double/hx2hy2
	rr		=  sqrt(x(1)**2+y(1)**2)
	tetnew		=  acos(x(1)/rr)
	if(y(1).lt.0._double)tetnew = 2._double*pi - tetnew
	tetnew		=  180._double/pi*tetnew
	numtou		=  0
	if(npaspt.lt.0)then
	  numtou  =  -nbdt/npaspt
	end if

	if(ficz.ne.'no') open(15,file =ficz,status='new')
	if(ficxy.ne.'no')open(16,file =ficxy,status='new')

	end

!	------------------------------------------------------------------

        subroutine calcb

        use lion_m
	implicit none
	integer ik,ikp1,jk,jkp1,k
	real(kind=double)vk

!       calcul de b

        teta	= acos(x/r)
        where(y.lt.0._double)
          teta = 2._double*pi - teta
        end where
        where (box)
          ivec          =  (r-rmin)/hr + 1._double
          jvec          =  teta/htet   + 1._double
          veca          =  r - rmin - (ivec-1._double)*hr
          vecb          =  teta - (jvec-1._double)*htet
          vecd1         =  veca-hr
          vecd2         =  vecb-htet
          vecp1         =  vecd1 * vecd2
          vecp2         = -veca  * vecd2
          vecp3         =  veca  * vecb
          vecp4         = -vecd1 * vecb
          vecz          =  z*z/2._double
        end where
        do k=1,np
          if(box(k))then
            ik          =  ivec(k)
            ikp1        =  ik + 1
            jk          =  jvec(k)
            jkp1        =  jk + 1
            vk          =  vecz(k)
            vecb1(k)    =  tbz(ik,jk)
            vecb2(k)    =  tbz(ikp1,jk)
            vecb3(k)    =  tbz(ikp1,jkp1)
            vecb4(k)    =  tbz(ik,jkp1)
            vec1(k)     = -vk*d2bz(ik,jk)
            vec2(k)     = -vk*d2bz(ikp1,jk)
            vec3(k)     = -vk*d2bz(ikp1,jkp1)
            vec4(k)     = -vk*d2bz(ik,jkp1)
            vecx1(k)    =  dbzx(ik,jk)
            vecx2(k)    =  dbzx(ikp1,jk)
            vecx3(k)    =  dbzx(ikp1,jkp1)
            vecx4(k)    =  dbzx(ik,jkp1)
            vecy1(k)    =  dbzy(ik,jk)
            vecy2(k)    =  dbzy(ikp1,jk)
            vecy3(k)    =  dbzy(ikp1,jkp1)
            vecy4(k)    =  dbzy(ik,jkp1)
          end if
        end do
        where (box)
          vec1          =  vec1 + vecb1
          vec2          =  vec2 + vecb2
          vec3          =  vec3 + vecb3
          vec4          =  vec4 + vecb4
          bz            =  (vecp1*vec1 + vecp2*vec2 +           &
     &                      vecp3*vec3 + vecp4*vec4)/hrhtet
          bx            =  z*(vecp1*vecx1 + vecp2*vecx2 +       &
     &                      vecp3*vecx3 + vecp4*vecx4)/hrhtet
          by            =  z*(vecp1*vecy1 + vecp2*vecy2 +       &
     &                      vecp3*vecy3 + vecp4*vecy4)/hrhtet
        end where
        where (r.le.rmin)
          bx    = 0._double
          by    = 0._double
          bz    = tbz(1,1)
        end where

	end

!	-------------------------------------------------------------------

        subroutine remp(m,n,ps)
        implicit none
	    integer,parameter :: double = selected_real_kind(P=10)

        integer iseed(12),i,j,m,n   ! Patrick Bertrand : 20/03/2017 : passage de 10 à 12 demandé par gfortran

        real(kind=double)ps(m,n),c(m,m),r(m,m),tr(m,m)
        real(kind=double)vtab(m,n-1),vec(n-1),vec1(n-1)
        real(kind=double)vmax
        logical  vmask(n-1)

        do i=1,12
          iseed = i
        end do
        call random_seed(put=iseed)
        vtab            =  10._double
100     vec             =  sum(vtab*vtab,dim=1)
        vmax            =  maxval(vec)
        if(vmax.le.1.)goto 200
        vmask           =  vec.gt.1._double
        do j=1,m
          call random_number(vec1)
          where(vmask)
            vtab(j,1:n-1) = 2._double*vec1 - 1._double
          end where
        end do
	goto 100
200     ps(1:m,2:n)     =  vtab
        ps(1:m,1)       =  0._double
        c               = 0._double
        do i=1,m
          do j=1,m
            c(i,j) = sum(ps(i,1:n)*ps(j,1:n))/n
          end do
        end do
        call choles(m,c,r)
        tr              =  transpose(r)
        ps(1,1:n)       =  ps(1,1:n)/tr(1,1)
        do i=2,m
          ps(i,1:n)     =  (ps(i,1:n) - sum(spread(tr(i,1:i-1),dim=2,  &
     &                     ncopies=n)*ps(1:i-1,1:n),dim=1))/tr(i,i)
        end do

        end

!       -----------------------------------------------------------------

        subroutine choles(m,a,r)
	implicit none
	integer,parameter :: double = selected_real_kind(P=10)
	integer i,m
        real(kind=double)a(m,m),r(m,m),vec(m)
        r               =  0._double
        r(1,1)          =  sqrt(a(1,1))
        r(1,2:m)        =  a(1,2:m)/r(1,1)
        do i=2,m-1
          r(i,i)        =  sqrt(a(i,i)-sum(r(1:i-1,i)**2))
          vec(i+1:m)    =  sum(spread(r(1:i-1,i),dim=2 &
     &    ,ncopies=m-i)*r(1:i-1,i+1:m),dim=1)
          r(i,i+1:m)    =  (a(i,i+1:m) - vec(i+1:m))/r(i,i)
        end do
        r(m,m)          =  sqrt(a(m,m)-sum(r(1:m-1,m)**2))
        end

!       -----------------------------------------------------------------

        subroutine initgra

	use lion_m
	implicit none
	real(kind=double)drmax,dprsp0max
	real(kind=double)dzmax,dpzsp0max
	real(kind=double)dphimax,dwswmax
	real(kind=double)xbamin,xbamax,ybamin,ybamax
	real(kind=double)xr1,yr1,xr2,yr2,omer
	real(kind=double)xd1,yd1,xd2,yd2
	real(kind=double)sxd,syd,xd,yd
	real(kind=double)xce,yce,rce,tetce1,tetce2,flacer
	real(kind=double)dtetce,tetce,xck,yck,xbary,ybary
	real(kind=double)gai1,wmpa
	integer nbpo,npbid,nxr,nyr,niter,nbrdro,num
	integer nbrcer,npcer,nf1000,imanq,nbr,i1,i2
	integer i,j,k,ier,ix,iy

	character*8 typcou
        character*25 titre(20),labelx(20),labely(20)
        integer nbx(20),nby(20)
        real(kind=double)xbid(20,2)
        real(kind=double)bornes(2,2,50)
	character*80 ficrep

	fenlog		= .true.
	ifreb		= 1
        write(6,*)'numero premier bidim (0 si aucun bidim)  :'
        read(5,*)numbid
        if(numbid.ne.0)then
	  nbx = 13
	  nby = 13
	  write(6,*)'frequence_de_trace par rapport aux pas de temps: '
	  read(5,*)ifreb
	  write(6,*)'tracé des centres ou des barycentres (c ou b) :'
	  read(5,*)baryc
	  ngravm	= npaspt/ifreb
	  if(ngravm.lt.0)ngravm = -ngravm
	  ngrav		= 0
      write(6,*)'boite des bidims (xmin xmax ymin ymax en mètres) :'
	  xbid(1,1)   	= -rmax
	  xbid(1,2)	= rmax
	  xbid(2,1)   	= -rmax
	  xbid(2,2)	= rmax
      read(5,*)((xbid(i,j),j=1,2),i=1,2)
      titre(1)      = 'trajectoire ion'
      labelx(1)     = 'x'
      labely(1)     = 'y'
      if(baryc.eq.'c')titre(2)      = 'centres de courbure'
      if(baryc.eq.'b')titre(2)      = 'barycentres'
      labelx(2)     = 'x'
      labely(2)     = 'y'
	  if(n.gt.1)then
            titre(3)      = 'trajectoires ions'
            labelx(3)     = 'x (m)'
            labely(3)     = 'y (m)'
            titre(4)      = 'correlation dphi-dw/w'
            labelx(4)     = 'dphi (deg)'
            labely(4)     = 'dw/w'
            titre(5)      = 'correlation r-dpr/p0'
            labelx(5)     = 'dr (mm)'
            labely(5)     = 'dpr/p0 (mrad)'
	  else
            titre(3)      = 'z(t) particule 1'
            labelx(3)     = 't'
            labely(3)     = 'z (mm)'
            titre(4)      = 'ez/1000000.*10.'
            labelx(4)     = 't'
            labely(4)     = 'ez'
            titre(5)      = '(vx*by-vy*bx)/1000000*10.'
            labelx(5)     = 't'
            labely(5)     = '(v*b)(z)'
	  endif
          titre(6)      = 'correlation r-dphi'
          labelx(6)     = 'dr (mm)'
          labely(6)     = 'dphi (deg)'
          titre(7)      = 'courbe dr(t)'
          labelx(7)     = 't'
          labely(7)     = 'dr (mm)'
          titre(8)      = 'correlation z-dpz/p0'
          labelx(8)     = 'dz (mm)'
          labely(8)     = 'dpz/p0 (mrad)'
          titre(9)      = 'courbe dz(t)'
          labelx(9)     = 't'
          labely(9)     = 'dz (mm)'
          titre(10)      = 'correlation r-z'
          labelx(10)     = 'dr (mm)'
          labely(10)     = 'dz (mm)'
          titre(11)      = 'correlation dphi-dpr/p0'
          labelx(11)     = 'dphi (deg)'
          labely(11)     = 'dpr/p0 (mrad)'
          titre(12)      = 'correlation dr-dw/w'
          labelx(12)     = 'dr (mm)'
          labely(12)     = 'dw/w'
          titre(13)      = 'correlation dpr/p0-dw/w'
          labelx(13)     = 'dpr/p0 (mrad)'
          labely(13)     = 'dw/w'
          titre(14)      = 'corr-grille r-dpr/p0'
          labelx(14)     = 'dr (mm)'
          labely(14)     = 'dpr/p0 (mrad)'
          titre(15)      = 'corr-grille dphi-dpr/p0'
          labelx(15)     = 'dphi (deg)'
          labely(15)     = 'dpr/p0 (deg)'
          titre(16)      = 'corr-grille r-dphi'
          labelx(16)     = 'dr (deg)'
          labely(16)     = 'dphi (deg)'
          titre(17)      = 'correlation dz-dw/w'
          labelx(17)     = 'dz (mm)'
          labely(17)     = 'dw/w'
	      itebid	 = 1
	  if(n.gt.1)then
	    write(6,*)'bidim dphi-dw/w : dphi_max dw/w_max   :'
	    read(5,*)dphimax,dwswmax
	    write(6,*)'bidim dr-dpr/p0 : dr_max   dpr/p0_max :'
	    read(5,*)drmax,dprsp0max
	    write(6,*)'bidim dz-dpz/p0 : dz_max   dpz/p0_max :'
	    read(5,*)dzmax,dpzsp0max
	    write(6,*)'itemax : ',nbdt
	    write(6,*)'debut_iteration_stockage_bidim (entre 1 et itemax) :'
	    read(5,*)itebid
	  end if
	  write(6,*)'valeurs de la fenetre des barycentres en mètres: '
	  read(5,*)xbamin,xbamax,ybamin,ybamax
	  bornes(1,1,1) = xbid(1,1)
	  bornes(2,1,1) = xbid(1,2)
	  bornes(1,2,1) = xbid(2,1)
	  bornes(2,2,1) = xbid(2,2)
	  bornes(1,1,2) = xbamin
	  bornes(2,1,2) = xbamax
	  bornes(1,2,2) = ybamin
	  bornes(2,2,2) = xbamax
	  if(n.eq.1)then
	    bornes(1,1,3) = 0._double
	    bornes(2,1,3) = (nbdt*1._double)/ifreb
	    bornes(1,2,3) = -dzmax
	    bornes(2,2,3) = dzmax
	    bornes(1,1,4) = 0._double
	    bornes(2,1,4) = (nbdt*1._double)/ifreb
	    bornes(1,2,4) = -dzmax
	    bornes(2,2,4) = dzmax
	    bornes(1,1,5) = 0._double
	    bornes(2,1,5) = (nbdt*1._double)/ifreb
	    bornes(1,2,5) = -dzmax
	    bornes(2,2,5) = dzmax
            call openbi(5,numbid,nbx,nby,titre,labelx,labely, &
     &                  bornes,ier)
	  else
	    bornes(1,1,3) = xbid(1,1)
	    bornes(2,1,3) = xbid(1,2)
	    bornes(1,2,3) = xbid(2,1)
	    bornes(2,2,3) = xbid(2,2)
	    bornes(1,1,4) = -dphimax
	    bornes(2,1,4) =  dphimax
	    bornes(1,2,4) = -dwswmax
	    bornes(2,2,4) =  dwswmax
	    bornes(1,1,5) = -drmax
	    bornes(2,1,5) =  drmax
	    bornes(1,2,5) = -dprsp0max
	    bornes(2,2,5) =  dprsp0max
	    bornes(1,1,6) = -drmax
	    bornes(2,1,6) =  drmax
	    bornes(1,2,6) = -dphimax
	    bornes(2,2,6) =  dphimax
	    bornes(1,1,7) =  0._double
	    bornes(2,1,7) =  (nbdt*1._double)/ifreb
	    bornes(1,2,7) = -drmax
	    bornes(2,2,7) =  drmax
	    bornes(1,1,8) = -dzmax
	    bornes(2,1,8) =  dzmax
	    bornes(1,2,8) = -dpzsp0max
	    bornes(2,2,8) =  dpzsp0max
	    bornes(1,1,9) =  0._double
	    bornes(2,1,9) =  (nbdt*1._double)/ifreb
	    bornes(1,2,9) = -dzmax
	    bornes(2,2,9) =  dzmax
	    bornes(1,1,10) =  -drmax
	    bornes(2,1,10) =   drmax
	    bornes(1,2,10) =  -dzmax
	    bornes(2,2,10) =   dzmax
	    bornes(1,1,11) =  -dphimax
	    bornes(2,1,11) =   dphimax
	    bornes(1,2,11) =  -dprsp0max
	    bornes(2,2,11) =   dprsp0max
	    bornes(1,1,12) =  -drmax
	    bornes(2,1,12) =   drmax
	    bornes(1,2,12) =  -dwswmax
	    bornes(2,2,12) =   dwswmax
	    bornes(1,1,13) =  -dprsp0max
	    bornes(2,1,13) =   dprsp0max
	    bornes(1,2,13) =  -dwswmax
	    bornes(2,2,13) =   dwswmax
	    bornes(1,1,14) =  -drmax
	    bornes(2,1,14) =   drmax
	    bornes(1,2,14) =  -dprsp0max
	    bornes(2,2,14) =   dprsp0max
	    bornes(1,1,15) =  -dphimax
	    bornes(2,1,15) =   dphimax
	    bornes(1,2,15) =  -dprsp0max
	    bornes(2,2,15) =   dprsp0max
	    bornes(1,1,16) =  -drmax
	    bornes(2,1,16) =   drmax
	    bornes(1,2,16) =  -dphimax
	    bornes(2,2,16) =   dphimax
        call openbi(17,numbid,nbx,nby,titre,labelx,labely, &
     &              bornes,ier)
	  endif
          if(ier.ne.0)then
            write(6,*)'bidim : openbi : stop erreur :',ier
            stop
          endif
          ax1(1)        = 8191._double/(xbid(1,2)-xbid(1,1))
          ax0(1)        = -ax1(1)*xbid(1,1)
          ay1(1)        = 8191._double/(xbid(2,2)-xbid(2,1))
          ay0(1)        = -ay1(1)*xbid(2,1)
          ax1(2)        = 8191._double/(xbamax-xbamin)
          ax0(2)        = -ax1(2)*xbamin
          ay1(2)        = 8191._double/(ybamax-ybamin)
          ay0(2)        = -ay1(2)*ybamin
	  if(n.eq.1)then
            ax1(7)        = (8191._double*ifreb)/nbdt
            ax0(7)        = 0._double
	    dzmax	  = 40._double
            ay1(9)        = 8191._double/2._double/dzmax
            ay0(9)        = ay1(9)*dzmax
	  endif
	  if(n.gt.1)then
            ax1(4)        = 8191._double/2._double/dphimax
            ax0(4)        = ax1(4)*dphimax
            ay1(4)        = 8191._double/2._double/dwswmax
            ay0(4)        = ay1(4)*dwswmax
            ax1(5)        = 8191._double/2._double/drmax
            ax0(5)        = ax1(5)*drmax
            ay1(5)        = 8191._double/2._double/dprsp0max
            ay0(5)        = ay1(5)*dprsp0max
            ax1(6)        = 8191._double/2._double/drmax
            ax0(6)        = ax1(6)*drmax
            ay1(6)        = 8191._double/2._double/dphimax
            ay0(6)        = ay1(6)*dphimax
            ax1(7)        = (8191._double*ifreb)/nbdt
            ax0(7)        = 0._double
            ay1(7)        = 8191._double/2._double/drmax
            ay0(7)        = ay1(7)*drmax
            ax1(8)        = 8191._double/2._double/dzmax
            ax0(8)        = ax1(8)*dzmax
            ay1(8)        = 8191._double/2._double/dpzsp0max
            ay0(8)        = ay1(8)*dpzsp0max
            ax1(9)        = (8191._double*ifreb)/nbdt
            ax0(9)        = 0._double
            ay1(9)        = 8191._double/2._double/dzmax
            ay0(9)        = ay1(9)*dzmax
	  end if
	  write(6,*)'lecture du fichier_des_reperes (lignes, piliers etc...) :'
	  read(5,*)ficrep
	  nbpo = 360
	  npbid = 0
	  if(ficrep.ne.'no')then
	    nbrpie = 0
	    open(24,file=ficrep)
	    read(24,*,err=500,end=500)xr1,yr1,xr2,yr2
	    read(24,*,err=500,end=500)nxr,nyr
	    read(24,*,err=500,end=500)niter,omer
10	    read(24,*,err=500,end=500)typcou
	    if(typcou.eq.'droite  ')then
	      read(24,*,err=500,end=500)nbrdro
	      do num=1,nbrdro
	        read(24,*,err=500,end=500)xd1,yd1,xd2,yd2
		    sxd  = (xd2-xd1)/(nbpo-1._double)
		    syd  = (yd2-yd1)/(nbpo-1._double)
	        do k=1,nbpo
	          xd  =  xd1 + (k-1._double)*sxd
	          yd  =  yd1 + (k-1._double)*syd
	          ix  =  ax1(1)*xd/1000._double + ax0(1)
	          iy  =  ay1(1)*yd/1000._double + ay0(1)
	          call incbid(numbid,ix,iy,ier)
	          if(n.ne.1)call incbid(numbid+2,ix,iy,ier)
	        end do
	      end do
	      npbid = npbid + nbrdro*nbpo
	      goto 10
	    elseif(typcou.eq.'cercle  ')then
	      read(24,*,err=500,end=500)nbrcer
	      do num=1,nbrcer
	        read(24,*,err=500,end=500)    &
     &        xce,yce,rce,tetce1,tetce2,npcer,flacer
	 	      tetce1	= tetce1*pi/180._double
	 	      tetce2	= tetce2*pi/180._double
	          dtetce = (tetce2-tetce1)/(nbpo-1._double)
	        if(flacer.lt.-10._double)then
	          nbrpie = nbrpie + 1
	          xpie(nbrpie) = xce/1000._double
	          ypie(nbrpie) = yce/1000._double
	          rpie(nbrpie) = rce/1000._double
	        endif
	        if(flacer.gt.-9000._double)then
	          do k=1,nbpo
	            tetce  = tetce1 + (k-1.)*dtetce
		        xd  =  xce + rce*cos(tetce)
		        yd  =  yce + rce*sin(tetce)
	            ix  =  ax1(1)*xd/1000._double + ax0(1)
	            iy  =  ay1(1)*yd/1000._double + ay0(1)
	            call incbid(numbid,ix,iy,ier)
	            if(n.ne.1)call incbid(numbid+2,ix,iy,ier)
	          end do
	          npbid = npbid + nbpo
	        endif
	      end do
	      goto 10
	    elseif(typcou.eq.'fintyp  ')then
	      goto 400
	    else
	      write(6,*)'erreur lecture fichier reperes...'
	      write(6,*)'typcou : ',typcou
	      stop
	    endif
500	    write(6,*)'erreur en lecture du fichier reperes...'
	    stop
400	    close(24)
	    if(nbrpie.ne.0)then
	      write(6,*)'attention : il y a des pieges a ions...'
	      write(6,*)'nombre de pieges : ',nbrpie
	    endif
	  end if
	  do i=0,8191,8
	    call incbid(numbid,i,4096,ier)
	    call incbid(numbid,4096,i,ier)
	    call incbid(numbid+1,i,4096,ier)
	    call incbid(numbid+1,4096,i,ier)
	    if(n.ne.1)then
	      call incbid(numbid+2,i,4096,ier)
	      call incbid(numbid+2,4096,i,ier)
	    endif
	  end do
	  do i=1,952
	    call incbid(numbid+1,4096,4096,ier)
	  end do
	  npbid 	= npbid + 2048
	  nf1000  	= npbid/1000
	  imanq		= npbid - nf1000*1000
	  do i=1,imanq
	    call incbid(numbid,4096,4096,ier)
	    if(n.ne.1)call incbid(numbid+2,4096,4096,ier)
	  end do

        endif
	if(nbrpie.ne.0)then
	  write(6,*)'nombre de cercles-pieges : ',nbrpie
	  write(6,*)'iteration max pour les pieges: '
      write(6,*)'(pour gagner du temps de calcul quand plus aucun piège...)'
	  read(5,*)itepiem
	  write(6,*)'flag de choix particules piegees: (0 ou 1)'
	  write(6,*)'0   : reste scotchée dans les pieges'
	  write(6,*)'1   : continue, confondue avec la particule centrale'
	  read(5,*)iflpie
	endif

        end

!       -----------------------------------------------------------------

	subroutine bid

	use lion_m
	implicit none
	real(kind=double)xck,yck,xbary,ybary
	real(kind=double)gai1,wmpa
	integer nbr,i1,k,ix,iy,ier,i

	if(ite.lt.itebid)return
	do k=1,np
          if(.not.fenlog(k))goto 310
	  ix	= ax1(1)*x(k) + ax0(1)
	  iy	= ay1(1)*y(k) + ay0(1)
	  ier = 0
	  if( (ix.ge.0 .and. ix .le. 8191)  .and.    &
     &        (iy.ge.0 .and. iy .le. 8191)  )then
	    if(k.eq.1)then
	      call incbid(numbid,ix,iy,ier)
	      if(np.gt.1)call incbid(numbid+2,ix,iy,ier)
	    end if
	    if(k.ne.1)call incbid(numbid+2,ix,iy,ier)
	  end if
	  if(k.eq.1)then
	    xck	= x(k) + py(k)*dt/qdtsm0c(k)/bz(k)
	    yck	= y(k) - px(k)*dt/qdtsm0c(k)/bz(k)
	    if(baryc.eq.'c')then
	      ix	= ax1(1)*xck + ax0(1)
	      iy	= ay1(1)*yck + ay0(1)
	      ier 	= 0
	      if( (ix.ge.0 .and. ix .le. 8191)  .and.    &
     &            (iy.ge.0 .and. iy .le. 8191)  )then
	        call incbid(numbid+1,ix,iy,ier)
	      end if
	    end if
	    if(baryc.eq.'b')then
	      ngrav	   = ngrav + 1
	      xgrav(ngrav) = xck
	      ygrav(ngrav) = yck
	      if(ngrav.eq.ngravm)then
	        ngrav   = 0
	        xbary	= sum(xgrav(1:ngravm))/ngravm
	        ybary	= sum(ygrav(1:ngravm))/ngravm
	        ix	= ax1(2)*xbary + ax0(2)
	        iy	= ay1(2)*ybary + ay0(2)
	        ier 	= 0
	        if( (ix.ge.0 .and. ix .le. 8191)  .and.  &
     &              (iy.ge.0 .and. iy .le. 8191)  )then
	          call incbid(numbid+1,ix,iy,ier)
	        end if
	      end if
	    end if
	  end if
310     continue
	end do

	if(np.eq.1)then
	    ix	=  (ax1(7)*ite)/ifreb + ax0(7)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    iy	=  ay1(9)*z(1)*1000. + ay0(9)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+2,ix,iy,ier)
	    iy	= ay1(9)*ez(1)/100000. + ay0(9)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+3,ix,iy,ier)
	    iy	= ay1(9)*clum/ga(1)* &
     &            (px(1)*by(1)-py(1)*bx(1))/100000. + ay0(9)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+4,ix,iy,ier)
	endif

	if(np.gt.1.and.nbrfam.eq.1)then
	  call caldxs
	  k		=  1
	  nbr		=  nbrion(k)
	  i1		=  index(k)
	  gai1		=  sqrt(1._double + p(i1)**2)
	  vec1(1:nbr)	=  dxs(5,1:nbr)*omhf/clum/p(i1)*180./pi*gai1
	  vec2(1:nbr)	=  dxs(6,1:nbr)*p(i1)/gai1/(gai1-1._double)
	  wmpa    	=  am0c2(k)*(gai1-1._double) &
     &                     /1.60217733e-13_double/a(k)
	  do i=1,nbr
            if(.not.fenlog(i))goto 320
	    ix	=  ax1(4)*vec1(i) + ax0(4)
	    iy	=  ay1(4)*vec2(i) + ay0(4)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+3,ix,iy,ier)
320         continue
	  end do
	  vec3(1:nbr)	=  dxs(1,1:nbr)*1000._double
	  vec4(1:nbr)	=  dxs(2,1:nbr)/p(i1)*1000._double
	  vec5(1:nbr)	=  dxs(3,1:nbr)*1000._double
	  vec6(1:nbr)	=  dxs(4,1:nbr)/p(i1)*1000._double
	  do i=1,nbr
            if(.not.fenlog(i))goto 330
	    ix	=  ax1(5)*vec3(i) + ax0(5)
	    iy	=  ay1(5)*vec4(i) + ay0(5)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+4,ix,iy,ier)
	    iy	=  ay1(4)*vec2(i) + ay0(4)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+11,ix,iy,ier)
	    ix	=  ay1(5)*vec4(i) + ay0(5)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    call incbid(numbid+12,ix,iy,ier)
	    ix	=  ax1(8)*vec5(i) + ax0(8)
	    iy	=  ay1(8)*vec6(i) + ay0(8)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+7,ix,iy,ier)
            iy  =  ay1(4)*vec2(i) + ay0(4)
            if(iy.lt.0)iy = 0
            if(iy.gt.8191)iy = 8191
            call incbid(numbid+16,ix,iy,ier)
330         continue
	  end do
	  do i=1,nbr
            if(.not.fenlog(i))goto 340
	    ix	=  ax1(6)*vec3(i) + ax0(6)
	    iy	=  ay1(6)*vec1(i) + ay0(6)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+5,ix,iy,ier)
340         continue
	  end do
	  ix	=  (ax1(7)*ite)/ifreb + ax0(7)
	  if(ix.lt.0)ix = 0
	  if(ix.gt.8191)ix = 8191
	  do i=1,nbr
            if(.not.fenlog(i))goto 350
	    iy	=  ay1(7)*vec3(i) + ay0(7)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+6,ix,iy,ier)
	    iy	=  ay1(9)*vec5(i) + ay0(9)
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+8,ix,iy,ier)
350         continue
	  end do
	  do i=1,nbr
            if(.not.fenlog(i))goto 360
	    ix	=  ax1(5)*vec3(i) + ax0(5)
	    iy	=  ax1(8)*vec5(i) + ax0(8)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+9,ix,iy,ier)
	    ix	=  ax1(4)*vec1(i) + ax0(4)
	    iy	=  ay1(5)*vec4(i) + ay0(5)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+10,ix,iy,ier)
360         continue
	  end do
	  call calgri
	  vec3(1:nbr) 	=  rgri(1:nbr)*1000._double
	  vec1(1:nbr) 	=  phigri(1:nbr)*180._double/pi
	  vec4(1:nbr) 	=  tprgri(1:nbr)*1000._double
	  do i=1,nbr
            if(.not.fenlog(i))goto 370
	    ix	=  ax1(5)*vec3(i) + ax0(5)
	    iy	=  ay1(5)*vec4(i) + ay0(5)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+13,ix,iy,ier)
	    ix	=  ax1(4)*vec1(i) + ax0(4)
	    iy	=  ay1(5)*vec4(i) + ay0(5)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+14,ix,iy,ier)
	    ix	=  ax1(6)*vec3(i) + ax0(6)
	    iy	=  ay1(6)*vec1(i) + ay0(6)
	    if(ix.lt.0)ix = 0
	    if(ix.gt.8191)ix = 8191
	    if(iy.lt.0)iy = 0
	    if(iy.gt.8191)iy = 8191
	    call incbid(numbid+15,ix,iy,ier)
370         continue
	  enddo
	end if

	end

!       -----------------------------------------------------------------

	subroutine caldxs

	use lion_m
	implicit none
	integer i1,i2,nbr,i,j,k
	real(kind=double)rc,ri1,tetai1,pri1,tpri1,co,si
	real(kind=double)ri1mm

        do k=1,nbrfam
          nbr           =  nbrion(k)
          i1            =  index(k)
          i2            =  index(k) + nbrion(k) - 1
          p(i1)         =  sqrt(px(i1)**2 + py(i1)**2 + pz(i1)**2)
          rc            =  -am0c(k)*p(i1)/q(i1)/bz(i1)
          dxd(1,1:nbr)  =  x(i1:i2) - x(i1)
          dxd(5,1:nbr)  =  y(i1:i2) - y(i1)
          ri1           =  sqrt(x(i1)*x(i1)+y(i1)*y(i1))
          tetai1        =  acos(x(i1)/ri1)
          if(y(i1).lt.0._double)tetai1 = 2._double*pi - tetai1
          pri1          = (x(i1)*px(i1) + y(i1)*py(i1))/ri1
          tpri1         = asin(pri1/p(i1))
          co            =  cos(tetai1 - tpri1)
          si            =  sin(tetai1 - tpri1)
          vec(1:nbr)    =   dxd(1,1:nbr)*co + dxd(5,1:nbr)*si
          dxd(5,1:nbr)  =  -dxd(1,1:nbr)*si + dxd(5,1:nbr)*co
          dxd(1,1:nbr)  =  vec(1:nbr)
          do i=1,nbr
            dxs(5,i)    =  rc*atan(dxd(5,i)/(rc+dxd(1,i)))
          end do
          vec(1:nbr)    =  sin(dxs(5,1:nbr)/2._double/rc)**2
          dxs(1,1:nbr)  =  (dxd(1,1:nbr) + 2._double*rc*vec(1:nbr))   &
     &                     /(1._double - 2._double*vec(1:nbr))
          vec1(1:nbr)   =  cos(dxs(5,1:nbr)/rc)
          vec2(1:nbr)   =  sin(dxs(5,1:nbr)/rc)
          dxd(2,1:nbr)  =  px(i1:i2) - px(i1)*vec1(1:nbr)         &
     &                               + py(i1)*vec2(1:nbr)
          dxd(6,1:nbr)  =  py(i1:i2) - px(i1)*vec2(1:nbr)         &
     &                               - py(i1)*vec1(1:nbr)
          vec(1:nbr)    =   dxd(2,1:nbr)*co + dxd(6,1:nbr)*si
          dxd(6,1:nbr)  =  -dxd(2,1:nbr)*si + dxd(6,1:nbr)*co
          dxd(2,1:nbr)  =  vec(1:nbr)
          dxs(2,1:nbr)  =  dxd(2,1:nbr)*vec1(1:nbr) +             &
     &                     dxd(6,1:nbr)*vec2(1:nbr)
          dxs(6,1:nbr)  =  - dxd(2,1:nbr)*vec2(1:nbr) +             &
     &                     dxd(6,1:nbr)*vec1(1:nbr)
          dxs(3,1:nbr)  =  z(i1:i2) - z(i1)
          dxs(4,1:nbr)  =  pz(i1:i2) - pz(i1)
          dxs(6,1:nbr)  =  sqrt((p(i1)+dxs(6,1:nbr))**2          &
     &                     + dxs(2,1:nbr)**2 + dxs(4,1:nbr)**2)   &
     &                     - p(i1)
	  do i=1,6
	    do j=1,6
	      cor(i,j)	=  sum(dxs(i,1:nbr)*dxs(j,1:nbr))/nbr
	    end do
	  end do
          tetai1        =  tetai1*180._double/pi
          ri1mm         =  ri1*1000._double
          t1            =  0._double
          t1(1,1)       =  -1._double
          t1(2,2)       =  -1._double
          do i=3,6
            t1(i,i)     =  1._double
          end do
          t1(1,6)       =  ri1/p(i1)
          t1(5,2)       =  -t1(1,6)
          cordc         =  matmul(cor,transpose(t1))
          cordc         =  matmul(t1,cordc)
	end do


	end

!       -----------------------------------------------------------------

        subroutine calgri

        use lion_m
	implicit none
	integer i1,i2,nbr,k
	real(kind=double)ux,uy

        do k=1,nbrfam
          nbr		=  nbrion(k)
          i1    	=  index(k)
          i2    	=  index(k) + nbrion(k) - 1
          p(i1:i2) 	=  sqrt(px(i1:i2)**2 + py(i1:i2)**2  &
     &                     + pz(i1:i2)**2)
          vec   	=  -am0c(k)*p(i1:i2)/q(i1:i2)/bz(i1)
          vec1(i1:i2)  	=  (x(i1:i2) + py(i1:i2)*dt/qdtsm0c(i1)/bz(i1))
          vec2(i1:i2)  	=  (y(i1:i2) - px(i1:i2)*dt/qdtsm0c(i1)/bz(i1))
	  ux		=  (x(i1)-vec1(i1))/vec(i1)
	  uy		=  (y(i1)-vec2(i1))/vec(i1)
	  vec3(i1:i2)	=  (x(i1)-vec1(i1:i2))*ux +  &
     &                      (y(i1)-vec2(i1:i2))*uy
	  vec4(i1:i2)	=  (x(i1)-vec1(i1:i2))**2 +  &
     &                     (y(i1)-vec2(i1:i2))**2 - vec(i1:i2)**2
	  rgri(i1:i2)	=  -vec3 + sqrt(abs(vec3**2 - vec4))
	  xgri(i1:i2)	=  x(i1) + rgri(i1:i2)*ux
	  ygri(i1:i2)	=  y(i1) + rgri(i1:i2)*uy
	  vec3(i1:i2)	=  sqrt( (xgri(i1:i2)-vec1(i1:i2))**2 + &
     &                           (ygri(i1:i2)-vec2(i1:i2))**2 )
	  vec4(i1:i2)	=  sqrt( (xgri(i1:i2)-vec1(i1))**2 + &
     &                           (ygri(i1:i2)-vec2(i1))**2 )
	  vec5(i1:i2)	=  ( (xgri(i1:i2)-vec1(i1:i2))*   &
     &                       (ygri(i1:i2)-vec2(i1)) -     &
     &                       (ygri(i1:i2)-vec2(i1:i2))*   &
     &                       (xgri(i1:i2)-vec1(i1)) )/    &
     &                       vec3(i1:i2)/vec4(i1:i2)
	  tprgri(i1:i2) =  asin(vec5(i1:i2))
	  vec4(i1:i2)	=  sqrt( (x(i1:i2)-vec1(i1:i2))**2 + &
     &                           (y(i1:i2)-vec2(i1:i2))**2 )
	  vec6(i1:i2)	=  ( (xgri(i1:i2)-vec1(i1:i2))*   &
     &                       (y(i1:i2)-vec2(i1:i2)) -     &
     &                       (ygri(i1:i2)-vec2(i1:i2))*   &
     &                       (x(i1:i2)-vec1(i1:i2)) )/    &
     &                       vec3(i1:i2)/vec4(i1:i2)
	  phigri(i1:i2) =  harm*asin(vec6(i1:i2))
	end do

	end

!       -----------------------------------------------------------------

	subroutine boite

	use lion_m
	implicit none

	r	=  sqrt(x*x+y*y)
	box     = r.gt.rmin  .and.  r.lt.(rmax-.02_double) .and. &
     &            z.gt.zmin  .and.  z.lt.zmax
	box	= box .and. ok
	nbox	=  count(box)

	end

!       -----------------------------------------------------------------

	subroutine box1

	use lion_m
	implicit none

        centre  = box .and. r.lt.rmax1 .and. abs(y).lt.(ymax1-hy1)
	nbox1	=  count(centre)

	end

!       -----------------------------------------------------------------

        subroutine box2

        use lion_m
	implicit none

	periph  =  vecx.gt.xmin2 .and. vecx.lt.xmax2  .and. &
     &             vecy.gt.ymin2 .and. vecy.lt.ymax2  .and. &
     &             box                                .and. &
     &             .not.centre
	nbox2	=  count(periph)

        end

!       -----------------------------------------------------------------

        subroutine boxdf(i)

        use lion_m
	implicit none
	integer i

	if(i.gt.1)goto 100
	if(iteext.eq.0)then
	  extra  = r.gt.rdf(i)  &
     &    .and.                                               &
     &    (((x-xdfi(i))*uydfi(i)-(y-ydfi(i))*uxdfi(i)).le.0._double) &
     &    .and.                                               &
     &    (((x-xdfo(i))*uydfo(i)-(y-ydfo(i))*uxdfo(i)).ge.0._double)
	  extra = extra .or.                                  &
     &    ( sortant                                           &
     &    .and.                                               &
     &    (((x-xdfi(i))*uydfi(i)-(y-ydfi(i))*uxdfi(i)).le.0._double) &
     &    .and.                                               &
     &    (((x-xdfo(i))*uydfo(i)-(y-ydfo(i))*uxdfo(i)).ge.0._double))
	  sortant  = sortant .or. extra
	else
          extra  = ite.ge.iteext  &
     &    .and.                                               &
     &    (((x-xdfi(i))*uydfi(i)-(y-ydfi(i))*uxdfi(i)).le.0._double) &
     &    .and.                                               &
     &    (((x-xdfo(i))*uydfo(i)-(y-ydfo(i))*uxdfo(i)).ge.0._double)
          extra = extra .or.                                  &
     &    ( sortant                                           &
     &    .and.                                               &
     &    (((x-xdfi(i))*uydfi(i)-(y-ydfi(i))*uxdfi(i)).le.0._double) &
     &    .and.                                               &
     &    (((x-xdfo(i))*uydfo(i)-(y-ydfo(i))*uxdfo(i)).ge.0._double))
          sortant = sortant .or. extra
	endif
	sortant   = sortant .and. ok
	extra 	  = extra .and. ok
        nextra    =  count(extra)

	return
100	extra = sortant                                     &
     &  .and.                                               &
     &  (((x-xdfi(i))*uydfi(i)-(y-ydfi(i))*uxdfi(i)).le.0._double) &
     &  .and.                                               &
     &  (((x-xdfo(i))*uydfo(i)-(y-ydfo(i))*uxdfo(i)).ge.0._double)
	extra    = extra .and. ok
        nextra   =  count(extra)

        end

!       -----------------------------------------------------------------

        subroutine boxdm(i)

        use lion_m
	implicit none
	integer i

	extra  = sortant  &
     &  .and.                                               &
     &  (((x-xdm(i))*uydmi(i)-(y-ydm(i))*uxdmi(i)).le.0._double) &
     &  .and.                                               &
     &  (((x-xdm(i))*uydmo(i)-(y-ydm(i))*uxdmo(i)).ge.0._double)
	extra    = extra .and. ok
        nextra   =  count(extra)

        end

!       -----------------------------------------------------------------

        subroutine boxgex(i)

        use lion_m
	implicit none
	integer i

	if(ndf.eq.0 .and. ndm.eq.0 .and. ngex.eq.1)then
	  if(ite.lt.iteext)then
	    franchi = .false.
	    return
	  endif
          franchi =  (indgex.eq.0) .and. ok .and.  	    &
     &    (((x-xgex(i))*uygex(i)-(y-ygex(i))*uxgex(i)).le.0._double)
	  where( (indgex.eq.1) .and. (((x-xgex(i))*uygex(i) &
     &           -(y-ygex(i))*uxgex(i)).gt.0._double) )
	    indgex = 0
	  end where
	  return
	endif

        franchi =  sortant	.and. ok  	  	    &
     &  .and.                                               &
	(indgex.eq.(i-1))		                    &
     &  .and.                                               &
     &  (((x-xgex(i))*uygex(i)-(y-ygex(i))*uxgex(i)).le.0._double)

        end

!       -----------------------------------------------------------------

        subroutine info

        use lion_m
	implicit none
	integer iprinf,nbr,i1,i2,k,i,j
	real(kind=double)potth,r0mm,pr0m,pz0m,dtnano,phihf,dph
	real(kind=double)exkvcm,eykvcm,ezkvcm,ekvcm,pp,ppp
	real(kind=double)eskvcm,erkvcm,pr,rc,ri1,pri1,co,ri1mm
	real(kind=double)fendphi,gai1,dr,drmm,ds,dsmm,dz,dzmm
	real(kind=double)dphi,dpr,drdpr,epsr,dzdpz,epsz
	real(kind=double)dp,dsdp,dwsw,dwmpa,dwdphi,drdc,drmmdc
	real(kind=double)dsdc,dsmmdc,dphidc,dprdc,drdprdc,epsrdc
	real(kind=double)dpdc,dsdpdc,dwswdc,dwmpadc,dwdphidc
	real(kind=double)pzf,copzf,alphaf,v0th,teta0d,z0mm,rr
	real(kind=double)wmevpa,xmm,ymm,zmm,rmm,xc1mm,yc1mm
	real(kind=double)tpr,tpz,tetai1,tpri1,si,wmpa
	data iprinf/1/

	if(iprinf.eq.1)then
	  iprinf	=0
          ga0th           =  ga(1)
          w0th            =  w0(1)
          w0mevth         =  w0th/1.60217733e-13_double/a(1)
          v0th            =  clum*sqrt((ga0th-1._double)* &
     &                       (ga0th+1._double))/ga0th
          r0th            =  v0th/omr
          b0th            =  -v0th*am0(1)*ga0th/q(1)/r0th
          ga1th           =  1._double/sqrt(1._double -  &
     &                       (ga0th*ga0th-1._double)  &
     &                     *r1th*r1th/ga0th/ga0th/r0th/r0th)
          w1th            =  am0c2(1)*(ga1th-1._double)
          w1mevth         =  w1th/1.6021773e-13_double/a(1)
          b1th            =  b0th/ga0th*ga1th
          potth           =  am0c2(1)*(ga1th-ga0th)/4._double/ntours/q(1)

          if(ficz.ne.'no')then
            write(15,*)'------------------------'
            write(15,*)'Programme LIONS'
	    write(15,*)'version  : decembre 1998'
            write(15,*)'------------------------'
            write(15,*)'-----------------------'
	    write(15,*)'informations selon z   '
            write(15,*)'-----------------------'
	    write(15,*)'fichier associe : ',ficxy(1:60)
	    write(15,4101)
	    write(15,*)' '
	    write(15,4111)
4111	    format('    ite   tr   teta      x       y       z',    &
     &             '       r        bx        by        bz',        &
     &             '    ez      tpz')
	    write(15,*)' '
	    write(15,4101)
	    write(15,*)' '
          end if
          if(ficxy.ne.'no')then
            write(16,*)'------------------------'
            write(16,*)'Programme LIONS'
	    write(16,*)'version  : decembre 1998'
            write(16,*)'------------------------'
            write(16,*)'-----------------------'
            write(16,*)corps(1:40)
            write(16,*)'-----------------------'
            write(16,*)'                        '
	    write(16,6001)amat(1),a(1),qn(1)
6001	    format(' masse atomique        : ',3x,f11.7,/,  &
     &             ' A                     : ',9x,f5.1,/,   &
     &             ' q                     : ',9x,f5.1)
            write(16,*)'                        '
	    r0mm		=  r0(1)*1000._double
	    teta0d		=  teta0(1)*180._double/pi
	    z0mm		=  z0(1)*1000._double
	    pr0m 		=  pr0(1)*1000._double
	    pz0m		=  pz0(1)*1000._double
	    write(16,6002)r0mm,teta0d,z0mm,pr0m,pz0m,w0mev1,phi0(1)
6002	    format(' r0                    : ',8x,f6.1,/,       &
     &             ' teta0                 : ',8x,f6.1,/,       &
     &             ' z0                    : ',8x,f6.1,/,       &
     &             ' pr0                   : ',6x,f8.2,/,       &
     &             ' pz0                   : ',6x,f8.2,/,       &
     &             ' w0                    : ',f14.10,/,        &
     &             ' phi0                  : ',3x,f11.2)
            write(16,*)'                        '
	    write(16,6003)fhf,harm,-pothfd,-pothfg,dbsb
6003	    format(' frequence hf          : ',5x,f9.0,/,       &
     &             ' harmonique            : ',5x,f9.0,/,       &
     &             ' potentiel hf 1        : ',5x,f9.1,/,       &
     &             ' potentiel hf 2        : ',5x,f9.1,/,       &
     &             ' correction db/b       : ',5x,f9.6)
            write(16,*)'                        '
            write(16,*)'fichier donnees  : ',ficin(1:80)
            write(16,*)'carte pot centre : ',ficv1(1:80)
            write(16,*)'carte pot periph : ',ficv2(1:80)
            write(16,*)'carte champ b    : ',ficb(1:100)
            write(16,*)'fichier correls  : ',ficsig(1:80)
            write(16,*)'                        '
            write(16,*)'------------------------------------------'
            write(16,*)'                        '
	    dtnano	= dt*1000000000.
	    write(16,6004)n,nbdt,npaspt,dtnano
6004	    format(' nb particules         : ',i14,/,       &
     &             ' nb iterations         : ',i14,/,       &
     &             ' nb iterations / tour  : ',i14,/,       &
     &             ' pas de temps          : ',5x,f9.6)
            write(16,*)'                        '
            write(16,*)'------------------------------------------'
            write(16,*)' '
	    write(16,*)' '
	    write(16,4101)
4101	    format(' -----------------------------------------',    &
     &             '----------------------------------------------',&
     &             '-----------------------------------')
	    write(16,*)' '
	    write(16,4100)
4100	    format('    ite   tr   teta      x       y      xc',    &
     &             '      yc        bz ',                           &
     &             '    es      er           r      tpr',&
     &             '       w           dphi')
	    write(16,*)' '
	    write(16,4101)
	    write(16,*)' '
	  end if
	end if
	tetold	=  tetnew
        rr 	=  sqrt(x(1)*x(1)+y(1)*y(1))
	tetnew	=  acos(x(1)/rr)
	if(y(1).lt.0._double)tetnew = 2._double*pi - tetnew
	tetnew	=  180._double/pi*tetnew
	if(npaspt.gt.0)then
	  if(tetnew.lt.tetold)numtou = numtou + 1
	else
	  if(tetnew.gt.tetold)numtou = numtou - 1
	end if
	tetall  =  tetnew + 360._double*(numtou-1)
	if(mod(ite,ifreq).ne.0  .and.  ite.ne.1   &
     &     .and. ite.ne.nbdt .and. ite.ne.(nbdt-1) )return

        phihf 	=  180._double/pi*omhf*t
        dph 	=  phihf-harm*tetall - 90._double
        wmevpa  = am0c2(1)*(ga(1)-1.)/1.60217733e-13/a(1)
	xmm	=  x(1)*1000._double
	ymm	=  y(1)*1000._double
	zmm     =  z(1)*1000._double
	exkvcm	=  ex(1)/100000._double
	eykvcm	=  ey(1)/100000._double
	ezkvcm	=  ez(1)/100000._double
	ekvcm	=  sqrt(ex(1)**2+ey(1)**2)/100000._double
	pp	=  sqrt(px(1)**2+py(1)**2)
	ppp	=  sqrt(px(1)**2+py(1)**2+pz(1)**2)
	eskvcm  =  (ex(1)*px(1)+ey(1)*py(1))/pp/100000._double
	erkvcm  =  (ex(1)*py(1)-ey(1)*px(1))/pp/100000._double
	ezkvcm  =  ez(1)/100000._double
	rmm	=  rr*1000._double
	xc1mm	=  (x(1) + py(1)*dt/qdtsm0c(1)/bz(1))*1000._double
	yc1mm	=  (y(1) - px(1)*dt/qdtsm0c(1)/bz(1))*1000._double
	pr	=  (xmm*px(1) + ymm*py(1))/rmm
	tpr	=  asin(pr/ppp)*1000._double
	tpz	=  asin(pz(1)/ppp)*1000._double

	if(mod(ite,ifreq).ne.0  .and.  ite.ne.1 )goto 800

	if(ficz.ne.'no')then
	  if(abs(z(1)).gt.coupnuz)then
	    anuz_cyclo = -1./am0(1)/ga(1)*(harm/omhf)**2/z(1)*  &
     &              q(1)*(ez(1) + clum/ga(1)*(px(1)*by(1)-py(1)*bx(1)))
	    anuz_cyclo = sqrt(abs(anuz_cyclo))
	  else
	    anuz_cyclo = 0.
	  endif
	  write(15,5900)ite,numtou,tetnew,xmm,ymm,zmm,rmm, &
     &                  bx(1),by(1),bz(1),ezkvcm,tpz,anuz_cyclo,pot_ver_vrai
5900	  format(1x,i6,2x,i3,2x,f5.1,4(1x,f8.2),        &
     &           3x,3(1x,f7.4),1x,f7.3,5x,f6.1,1x,4(1x,f12.5))
	end if

	if(ficxy.ne.'no')then
	  write(16,6000)ite,numtou,tetnew,xmm,ymm,xc1mm,yc1mm,  &
     &                 bz(1),eskvcm,erkvcm,rmm,tpr,wmevpa,dph
6000	  format(1x,i7,2x,i3,2x,f5.1,4(1x,f7.1),        &
     &          3x,f7.4,2(1x,f7.3),5x,f7.1,2x,f6.1,1x,    &
     &          f9.5,1x,f15.2)
	end if

	if(iflcor.eq.0.or.(ficcor.eq.'no'.and.ficmat.eq.'no'   &
     &     .and.ficdc.eq.'no'))goto 800
	do k=1,nbrfam
	  nbr		=  nbrion(k)
          i1            =  index(k)
          i2            =  index(k) + nbrion(k) - 1
	  p(i1)		=  sqrt(px(i1)**2 + py(i1)**2 + pz(i1)**2)
          rc            =  -am0c(k)*p(i1)/q(i1)/bz(i1)
	  dxd(1,1:nbr)  =  x(i1:i2) - x(i1)
	  dxd(5,1:nbr)  =  y(i1:i2) - y(i1)
          ri1 		=  sqrt(x(i1)*x(i1)+y(i1)*y(i1))
	  tetai1	=  acos(x(i1)/ri1)
	  if(y(i1).lt.0.)tetai1 = 2._double*pi - tetai1
	  pri1		= (x(i1)*px(i1) + y(i1)*py(i1))/ri1
	  tpri1 	= asin(pri1/p(i1))
	  co		=  cos(tetai1 - tpri1)
	  si		=  sin(tetai1 - tpri1)
	  vec(1:nbr)    =   dxd(1,1:nbr)*co + dxd(5,1:nbr)*si
	  dxd(5,1:nbr)  =  -dxd(1,1:nbr)*si + dxd(5,1:nbr)*co
	  dxd(1,1:nbr)	=  vec(1:nbr)
	  do i=1,nbr
	    dxs(5,i) 	=  rc*atan(dxd(5,i)/(rc+dxd(1,i)))
	  end do
	  vec(1:nbr)  	=  sin(dxs(5,1:nbr)/2._double/rc)**2
	  dxs(1,1:nbr)	=  (dxd(1,1:nbr) + 2._double*rc*vec(1:nbr))   &
     &                     /(1._double - 2._double*vec(1:nbr))
          vec1(1:nbr)   =  cos(dxs(5,1:nbr)/rc)
          vec2(1:nbr)   =  sin(dxs(5,1:nbr)/rc)
	  dxd(2,1:nbr)  =  px(i1:i2) - px(i1)*vec1(1:nbr)         &
     &                               + py(i1)*vec2(1:nbr)
	  dxd(6,1:nbr)  =  py(i1:i2) - px(i1)*vec2(1:nbr)         &
     &                               - py(i1)*vec1(1:nbr)
	  vec(1:nbr)    =   dxd(2,1:nbr)*co + dxd(6,1:nbr)*si
	  dxd(6,1:nbr)  =  -dxd(2,1:nbr)*si + dxd(6,1:nbr)*co
	  dxd(2,1:nbr)	=  vec(1:nbr)
	  dxs(2,1:nbr)	=  dxd(2,1:nbr)*vec1(1:nbr) +             &
     &                     dxd(6,1:nbr)*vec2(1:nbr)
	  dxs(6,1:nbr)	=  - dxd(2,1:nbr)*vec2(1:nbr) +             &
     &                     dxd(6,1:nbr)*vec1(1:nbr)
	  dxs(3,1:nbr)	=  z(i1:i2) - z(i1)
	  dxs(4,1:nbr)  =  pz(i1:i2) - pz(i1)
	  dxs(6,1:nbr)  =  sqrt((p(i1)+dxs(6,1:nbr))**2          &
     &                     + dxs(2,1:nbr)**2 + dxs(4,1:nbr)**2)   &
     &                     - p(i1)
	  cor           =  0._double
	  do i=1,6
	    do j=1,6
	      cor(i,j)  =  sum(dxs(i,1:nbr)*dxs(j,1:nbr),mask=ok)/nbrok
	    end do
	  end do
	  tetai1	=  tetai1*180._double/pi
	  ri1mm		=  ri1*1000._double
	  t1		=  0._double
	  t1(1,1)	=  -1._double
	  t1(2,2)	=  -1._double
	  do i=3,6
	    t1(i,i)	=  1._double
	  end do
	  t1(1,6)	=  ri1/p(i1)
	  t1(5,2)	=  -t1(1,6)
	  cordc		=  matmul(cor,transpose(t1))
	  cordc		=  matmul(t1,cordc)

	  if(ficmat.ne.'no')then
	    write(18,*)' '
	    write(18,*)'ite numtou r tetai1 : ',ite,numtou,ri1mm,tetai1
	    write(18,*)' '
	    do i=1,6
	      write(18,2018)(cor(i,j),j=1,6)
2018	      format(6(1x,e15.8))
	    end do
	    write(18,*)' '
	    call calgric
	    do i=1,6
	      do j=1,6
	        corgri(i,j)  =  sum(dxsgri(i,1:nbr)*dxsgri(j,1:nbr) &
     &                          ,mask=ok)/nbrok
	      end do
	    end do
	    do i=1,6
	      write(18,2018)(corgri(i,j),j=1,6)
	    end do
	    if(ite.eq.nbdt)then
	      write(6,*)'fenetre en dw/w : '
	      read(5,*)fendwsw
	      write(6,*)'fenetre en dphi : '
	      read(5,*)fendphi
	      if(fendwsw.gt.0.)then
	        fenlog(1:nbr) = abs(dxs(6,1:nbr)*2/p(i1)) .le. fendwsw
	        gai1          = sqrt(1. + p(i1)*p(i1))
	        vec(1:nbr)    = dxs(5,1:nbr)*omhf/clum/p(i1)*180.  &
     &                          /pi*gai1
	        fenlog(1:nbr) = (abs(vec(1:nbr)).le.fendphi) .and.  &
     &                          fenlog(1:nbr)
	        fenlog(1:nbr) = fenlog(1:nbr) .and. ok(1:nbr)
	        nfen	= count(fenlog(1:nbr))
	        write(6,*)'nb. particules filtrees par fenètre : ',nfen
	        corfen  =  0._double
	        do i=1,6
	          do j=1,6
	            corfen(i,j)  =  sum(dxs(i,1:nbr)*dxs(j,1:nbr),  &
     &                              mask=fenlog(1:nbr))/nfen
	          end do
	        end do
	        do i=1,6
	          do j=1,6
	            corfeng(i,j) =  sum(dxsgri(i,1:nbr)*  &
     &                              dxsgri(j,1:nbr),      &
     &                              mask=fenlog(1:nbr))/nfen
	          end do
	        end do
	        write(18,*)' '
	        write(18,*)'matrices fitrees en t et grille  : '
	        write(18,*)'fenetre dw/w   : ',fendwsw
	        write(18,*)'nb particules  : ',nfen
	        write(18,*)' '
	        do i=1,6
	          write(18,2018)(corfen(i,j),j=1,6)
	        end do
	        write(18,*)' '
	        do i=1,6
	          write(18,2018)(corfeng(i,j),j=1,6)
	        end do
	      endif
	      call outinf
	    endif
	  end if

	  if(ficcor.ne.'no')then
	    dr		=  sqrt(cor(1,1))*rac8
	    drmm	=  dr*1000._double
	    ds		=  sqrt(cor(5,5))*rac8
	    dsmm	=  ds*1000._double
	    dz		=  sqrt(cor(3,3))*rac8
	    dzmm	=  dz*1000._double
	    gai1	=  sqrt(1._double + p(i1)*p(i1))
	    dphi	=  ds*omhf/clum/p(i1)*180._double/pi*gai1
	    dpr		=  sqrt(cor(2,2))*rac8
	    drdpr	=  8._double*sqrt(cor(1,1)*cor(2,2) -     &
     &                                cor(1,2)**2)
	    epsr	=  drdpr*1000000._double/p(i1)
	    dzdpz	=  8.*sqrt(cor(3,3)*cor(4,4) -     &
     &                                cor(3,4)**2)
	    epsz	=  dzdpz*1000000._double/p(i1)
	    dp		=  sqrt(cor(6,6))*rac8
	    dsdp	=  8._double*sqrt(cor(5,5)*cor(6,6) -     &
     &                                cor(5,6)**2)
	    dwsw	=  dp*p(i1)/gai1/(gai1-1._double)
	    wmpa	=  am0c2(k)*(gai1-1._double)/ &
     &                     1.60217733e-13_double/a(k)
	    dwmpa	=  dwsw*wmpa
	    dwdphi	=  omhf/clum/p(i1)*180._double/pi*gai1 * &
     &                     p(i1)/gai1/(gai1-1.)*wmpa *     &
     &                     dsdp*1000._double
	    write(17,2017)numtou,ri1mm,tetai1,drmm,dsmm,dphi,    &
     &                    epsr,dwsw,wmpa,dwmpa,dwdphi,dzmm,epsz
2017	    format(1x,i3,2x,f7.2,2x,f5.1,1x,f6.2,1x,f6.2,1x,f7.2, &
     &             1x,f6.2,1x,e12.7,1x,f6.3,1x,e12.7,1x,e12.7,    &
     &             2x,f6.2,1x,f6.2)
	  end if

          if(ficdc.ne.'no')then
            drdc        =  sqrt(cordc(1,1))*rac8
            drmmdc      =  drdc*1000._double
            dsdc        =  sqrt(cordc(5,5))*rac8
            dsmmdc      =  dsdc*1000._double
            gai1        =  sqrt(1. + p(i1)*p(i1))
            dphidc      =  dsdc*omhf/clum/p(i1)*180./pi*gai1
            dprdc       =  sqrt(cordc(2,2))*rac8
            drdprdc     =  8._double*sqrt(cordc(1,1)*cordc(2,2) -     &
     &                                cordc(1,2)**2)
            epsrdc      =  drdprdc*1000000._double/p(i1)
            dpdc        =  sqrt(cordc(6,6))*rac8
            dsdpdc      =  8._double*sqrt(cordc(5,5)*cordc(6,6) -     &
     &                                cordc(5,6)**2)
            dwswdc      =  dpdc*p(i1)/gai1/(gai1-1._double)
            wmpa        =  am0c2(k)*(gai1-1._double)/    &
     &                     1.60217733e-13_double/a(k)
            dwmpadc     =  dwswdc*wmpa
            dwdphidc    =  omhf/clum/p(i1)*180._double/pi*gai1 *       &
     &                     p(i1)/gai1/(gai1-1.)*wmpa * dsdpdc*1000.
            write(19,2019)numtou,ri1mm,tetai1,drmmdc,dsmmdc,dphidc,    &
     &                    epsrdc,dwswdc,wmpa,dwmpadc,dwdphidc
2019        format(1x,i3,2x,f7.2,2x,f5.1,1x,f6.2,1x,f6.2,1x,f5.2, &
     &             1x,f6.2,1x,e12.7,1x,f6.3,1x,e12.7,1x,e12.7)
          end if

	end do

800	if(ite.eq.(nbdt-1))then
	  xf    =  x(1)
	  yf    =  y(1)
	  zf	=  z(1)
	  rf    =  rr
	  tetf  =  tetnew
 	  phif  =  phihf
          if(ficsto.ne.'no')then
            write(6,*)'attention : sauvegarde paquet d ions...'
            open(26,file=ficsto,form='unformatted',convert='BIG_ENDIAN',status='new')
            write(26)(x(i),y(i),z(i),i=1,n)
          endif
	endif
	if(ite.eq.nbdt)then
	  pxf   = px(1)
	  pyf   = py(1)
	  pzf	= pz(1)
	  pf    = sqrt(pxf**2 + pyf**2 + pzf**2)
	  copzf = cos(pzf)
	  alphaf = acos(pxf/pf/copzf)
	  if(pyf.lt.0._double)alphaf = 2._double*pi - alphaf
	  tprf	= ((pi/2.+tetf*pi/180.)-alphaf)*1000._double
	  wf    = wmevpa
	  tpzf  = asin(pzf/pf)*1000._double
          if(ficsto.ne.'no')then
            write(26)(px(i),py(i),pz(i),i=1,n)
            close(26)
          endif
          if(ficxy.ne.'no')then
	    write(16,*)'---------------------------------'
	    write(16,*)'info pour le mouvement inverse...'
	    write(16,*)'---------------------------------'
	    write(16,*)'r0 teta0 z0 pr0 pz0 w0mev phi0 :'
	    write(16,*)rf,tetf,zf,tprf
	    write(16,*)tpzf,wf,phif
	    write(16,*)'---------------------------------'
	    write(16,*)'------------------------------------'
	    write(16,*)'info pour continuer le mouvement ...'
	    write(16,*)'------------------------------------'
	    write(16,*)'r0 teta0 z0 pr0 pz0 w0mev phi0 :'
	    write(16,*)rr,tetnew,z(1),tprf
	    write(16,*)tpzf,wf,phihf
	    write(16,*)'---------------------------------'
	    write(16,*)'------------------------------------'
	    write(16,*)'info utile pour DEDUCE_SIGMA ...'
	    write(16,*)'------------------------------------'
	    write(16,*)'amat a qn w0mev : '
	    write(16,*)amat(1),a(1),qn(1),wf
	    write(16,*)'champ magnetique (val. absolue) :'
	    write(16,*)-bz(1)
	    write(16,*)'fhf harmonique :'
	    write(16,*)fhf,harm
	    write(16,*)'------------------------------------'
	  endif
	endif

        end

!       -----------------------------------------------------------------

	subroutine outinf
        use lion_m
	implicit none
	integer i,j,ii,npfil
	real(kind=double)r0mm,teta0d,z0mm,pr0m,pz0m
	real(kind=double)dtnano,b00,p00,am0c00,r00

	write(6,*)'nom fichier pour inflecteur : (no si aucun)'
	read(5,*)ficinf
	if(ficinf.eq.'no')return
	call bid
	open(31,file=ficinf,status='new')
        write(31,*)'-----------------------'
        write(31,*)corps(1:40)
        write(31,*)'-----------------------'
        write(31,*)'                        '
        write(31,6001)amat(1),a(1),qn(1)
6001    format(' masse atomique        : ',3x,f11.7,/,  &
     &         ' A                     : ',9x,f5.1,/,   &
     &         ' q                     : ',9x,f5.1)
        write(31,*)'                        '
        r0mm                =  r0(1)*1000.
        teta0d              =  teta0(1)*180./pi
        z0mm                =  z0(1)*1000.
        pr0m                =  pr0(1)*1000.
        pz0m                =  pz0(1)*1000.
        write(31,6002)r0mm,teta0d,z0mm,pr0m,pz0m,w0mev1,phi0(1)
6002    format(' r0                    : ',8x,f6.1,/,       &
     &         ' teta0                 : ',8x,f6.1,/,       &
     &         ' z0                    : ',8x,f6.1,/,       &
     &         ' pr0                   : ',6x,f8.2,/,       &
     &         ' pz0                   : ',6x,f8.2,/,       &
     &         ' w0                    : ',f14.10,/,        &
     &         ' phi0                  : ',3x,f11.2)
        write(31,*)'                        '
        write(31,6003)fhf,harm,pothf,dbsb
6003    format(' frequence hf          : ',5x,f9.0,/,       &
     &         ' harmonique            : ',5x,f9.0,/,       &
     &         ' potentiel hf          : ',5x,f9.1,/,       &
     &         ' correction db/b       : ',5x,f9.6)
        write(31,*)'                        '
        write(31,*)'fichier donnees  : ',ficin(1:80)
        write(31,*)'carte pot centre : ',ficv1(1:80)
        write(31,*)'carte pot periph : ',ficv2(1:80)
        write(31,*)'carte champ b    : ',ficb(1:100)
        write(31,*)'                        '
        write(31,*)'------------------------------------------'
        write(31,*)'                        '
        dtnano      = dt*1000000000._double
        write(31,6004)nbdt,npaspt,dtnano
6004    format(' nb iterations               ',i14,/,       &
     &         ' nb iterations par tour      ',i14,/,       &
     &         ' pas de temps                ',5x,f9.6)
        write(31,*)'                        '
        write(31,*)'------------------------------------------'
        write(31,*)' '
        write(31,*)' '
        write(31,4101)
4101    format(' -----------------------------------------',    &
     &         '----------------------------------------------',&
     &         '-----------------------------------')
        write(31,*)' '
	write(31,*)'-------------------------'
	write(31,*)'fenetrage en dw/w      : '
	write(31,*)'-------------------------'
	write(31,*)'fenetre                 :  ',fendwsw
	write(31,*)'nb. particules filtrees :  ',nfen
	write(31,*)' '
        p00	=  sqrt(px(1)**2 + py(1)**2  &
     &                     + pz(1)**2)
	b00	=  bz(1)
	am0c00	=  am0c(1)
        r00   	=  -am0c(1)*p00/q(1)/b00
	write(31,*)'bz local : ',b00
	write(31,*)'r0       : ',r00
	write(31,*)'p0       : ',p00
	write(31,*)'m0*c     : ',am0c00
	write(31,*)' '
	write(31,*)'matrice faisceau correspondante en t : '
	write(31,*)' '
	do i=1,6
	  write(31,2018)(corfen(i,j),j=1,6)
2018	  format(6(1x,e15.8))
	end do
	write(31,*)' '
	write(31,*)'matrice faisceau correspondante en grillt : '
	write(31,*)' '
	do i=1,6
	  write(31,2018)(corfeng(i,j),j=1,6)
	end do
	write(31,*)' '
	write(31,*)'liste des particules :'
	write(31,*)' '
	ii	= 1
	do npfil=1,nfen
555	  continue
	  if(.not.fenlog(ii))then
	    ii = ii + 1
	    if(ii.gt.np)then
	      write(6,*)'pb fenlog...'
	      stop
	    endif
	    goto 555
	  endif
	  write(31,2555)npfil,x(ii),y(ii),z(ii),px(ii),py(ii), &
     &                  pz(ii)
2555	  format(1x,i4,3(1x,e11.4),3(1x,e12.5))
	  ii = ii + 1
	enddo

	end

!       -----------------------------------------------------------------

        subroutine calgric

        use lion_m
	implicit none
	integer i1,i2,nbr,k
	real(kind=double)ux,uy

        do k=1,nbrfam
          nbr		=  nbrion(k)
          i1    	=  index(k)
          i2    	=  index(k) + nbrion(k) - 1
          p(i1:i2) 	=  sqrt(px(i1:i2)**2 + py(i1:i2)**2  &
     &                     + pz(i1:i2)**2)
          vec   	=  -am0c(k)*p(i1:i2)/q(i1:i2)/bz(i1)
          vec1(i1:i2)  	=  (x(i1:i2) + py(i1:i2)*dt/qdtsm0c(i1)/bz(i1))
          vec2(i1:i2)  	=  (y(i1:i2) - px(i1:i2)*dt/qdtsm0c(i1)/bz(i1))
	  ux		=  (x(i1)-vec1(i1))/vec(i1)
	  uy		=  (y(i1)-vec2(i1))/vec(i1)
	  vec3(i1:i2)	=  (x(i1)-vec1(i1:i2))*ux +  &
     &                      (y(i1)-vec2(i1:i2))*uy
	  vec4(i1:i2)	=  (x(i1)-vec1(i1:i2))**2 +  &
     &                     (y(i1)-vec2(i1:i2))**2 - vec(i1:i2)**2
	  rgri(i1:i2)	=  -vec3 + sqrt(abs(vec3**2 - vec4))
	  xgri(i1:i2)	=  x(i1) + rgri(i1:i2)*ux
	  ygri(i1:i2)	=  y(i1) + rgri(i1:i2)*uy
	  vec3(i1:i2)	=  sqrt( (xgri(i1:i2)-vec1(i1:i2))**2 + &
     &                           (ygri(i1:i2)-vec2(i1:i2))**2 )
	  vec4(i1:i2)	=  sqrt( (xgri(i1:i2)-vec1(i1))**2 + &
     &                           (ygri(i1:i2)-vec2(i1))**2 )
	  vec5(i1:i2)	=  ( (xgri(i1:i2)-vec1(i1:i2))*   &
     &                       (ygri(i1:i2)-vec2(i1)) -     &
     &                       (ygri(i1:i2)-vec2(i1:i2))*   &
     &                       (xgri(i1:i2)-vec1(i1)) )/    &
     &                       vec3(i1:i2)/vec4(i1:i2)
	  tprgri(i1:i2) =  asin(vec5(i1:i2))
	  vec4(i1:i2)	=  sqrt( (x(i1:i2)-vec1(i1:i2))**2 + &
     &                           (y(i1:i2)-vec2(i1:i2))**2 )
	  vec6(i1:i2)	=  ( (xgri(i1:i2)-vec1(i1:i2))*   &
     &                       (y(i1:i2)-vec2(i1:i2)) -     &
     &                       (ygri(i1:i2)-vec2(i1:i2))*   &
     &                       (x(i1:i2)-vec1(i1:i2)) )/    &
     &                       vec3(i1:i2)/vec4(i1:i2)
	  phigri(i1:i2) =  harm*asin(vec6(i1:i2))
	end do

	dxsgri(1,1:nbr) = rgri(1:nbr)
	dxsgri(2,1:nbr) = tprgri(1:nbr)*p(1)
	dxsgri(3,1:nbr) = dxs(3,1:nbr)
	dxsgri(4,1:nbr) = dxs(4,1:nbr)
	dxsgri(5,1:nbr) = phigri(1:nbr)*vec(1:nbr)/harm
	dxsgri(6,1:nbr) = dxs(6,1:nbr)

	end

	subroutine openbi(nbr,n1,nbx,nby,titre,labelx,labely, &
     &                    bornes,ier)
	integer,parameter :: double = selected_real_kind(P=10)
	character*10 type,lieu,auteur,forma
	character*25 titre(*),labelx(*),labely(*)
	character*120 vide
	character*80 filnam
	character*20 cxmi(50),cxma(50),cymi(50),cyma(50)
	integer idate(3),itime(3)
	integer nbx(*),nby(*)
	real(kind=double)bornes(2,2,50)
	common /carbid/nbrbid,n1bid,n2bid
	common /carbid/nbxbid(50),nbybid(50),lmxbid(50),lmybid(50)
	common /carbid/ind(50)
	ier=0
	type='bidim'
	lieu='ganil'
	auteur='inconnu'
	forma='Binaire'
	nbrbid=0
	idate = 0
	itime = 0
	do i=1,200
	  vide(i:i) = ' '
	end do
	if(nbr.le.0.or.nbr.gt.50)then
	  write(6,*)'openbi : nombre de bidims incorrect'
	  ier=1
	  return
	endif
	if(n1.le.0)then
	  write(6,*)'openbi : numero du premier bidim incorrect'
	  ier=2
	  return
	endif
	n1bid=n1
	n2bid=n1+nbr-1
	do 10 i=1,nbr
	mx=2**nbx(i)-1
	if(mx.gt.8191)then
	  write(6,1001)i,nbx(i)
1001	  format(' openbi : pour bidim ',i10,/,   &
     &    ' nombre de bits en x incorrect : ',i10)
	  ier=3
	  return
	endif
	nbxbid(i)=nbx(i)
	lmxbid(i)=mx
	my=2**nby(i)-1
	if(my.gt.8191)then
	  write(6,1002)i,nby(i)
1002	  format(' openbi : pour bidim ',i10,/,  &
     &    ' nombre de bits en y incorrect : ',i10)
	  ier=4
	  return
	endif
	nbybid(i)=nby(i)
	lmybid(i)=my
10	continue
	do 20 i=1,nbr
	num=n1+i-1
	if(num.lt.10)then
	write(filnam,'(i1,4h.bid)')num
	else if(num.lt.100)then
	write(filnam,'(i2,4h.bid)')num
	else if(num.lt.1000)then
	write(filnam,'(i3,4h.bid)')num
	else if(num.lt.10000)then
	write(filnam,'(i4,4h.bid)')num
	else if(num.lt.100000)then
	write(filnam,'(i5,4h.bid)')num
	else if(num.lt.1000000)then
	write(filnam,'(i6,4h.bid)')num
	endif
	write(6,*)'filnam : ',filnam
	open(40+i,file=filnam,form='UNFORMATTED',status='NEW', &
     &  err=800)
	minx=0
	maxx=lmxbid(i)
	miny=0
	maxy=lmybid(i)
	write(cxmi(i),'(d20.10)')bornes(1,1,i)
	write(cxma(i),'(d20.10)')bornes(2,1,i)
	write(cymi(i),'(d20.10)')bornes(1,2,i)
	write(cyma(i),'(d20.10)')bornes(2,2,i)
	write(40+i)num,type,idate,itime,lieu,auteur,titre(i),labelx(i), &
     &  labely(i),maxx,maxy,forma,cxmi(i),cxma(i),cymi(i),cyma(i),vide
	ind(i)=1
20	continue
	write(6,*)'openbi : open des bidims sur disque effectue correctement'
	nbrbid=nbr
	return
800	write(6,1003)i
1003	format('openbi : probleme open du bidim : ',i10)
	ier=5
	if(i.eq.1)return
	do 30 k=1,i-1
	close(40+i)
30	continue
	return
	end
!
!	---------------------------------------------------------------
!
	subroutine incbid(n,ix,iy,ier)
	common /carbid/nbrbid,n1bid,n2bid
	common /carbid/nbxbid(50),nbybid(50),lmxbid(50),lmybid(50)
	common /carbid/ind(50)
	common /bufbid/ibubid(4096,50)
	ier=0
	if(n.lt.n1bid.or.n.gt.n2bid)then
	  write(6,*)'incbid : numero bidim incorrect :',n
	  write(6,*)'...en dehors de : ',n1bid,n2bid
	  ier=1
	  return
	endif
	num=n-n1bid+1
	if(ix.lt.0)ix = 0
	if(ix.gt.lmxbid(num))ix = lmxbid(num)
	if(iy.lt.0)iy = 0
	if(iy.gt.lmybid(num))iy = lmybid(num)
	if(ind(num).eq.4097)then
	  call wribid(40+num,ibubid(1,num))
	  ind(num)=1
	endif
	ibubid(ind(num),num)=ix*65536 + iy
	ind(num)=ind(num)+1
	end
!
!	---------------------------------------------------------------
!
	subroutine clobid(ier)
	common /carbid/nbrbid,n1bid,n2bid
	common /carbid/nbxbid(50),nbybid(50),lmxbid(50),lmybid(50)
	common /carbid/ind(50)
	common /bufbid/ibubid(4096,50)
	ier=0
	do 10 num=1,nbrbid
	if(ind(num).eq.4097)then
	  call wribid(40+num,ibubid(1,num))
	  ind(num)=1
	else if(ind(num).ne.1)then
	  ibubid(ind(num),num)=9999*65536 + 9999
	  call wribid(40+num,ibubid(1,num))
	endif
	close(40+num)
10	continue
	end
!
!	---------------------------------------------------------------
!
        subroutine wribid(iu,ibuf)
	integer iu,ibuf(4096)
        write(iu)ibuf(1:4096)
        end
!
!	---------------------------------------------------------------
!
	subroutine addh(rmin,rmax,nr,ntet,tbz,d2bz,dbzx,dbzy,ncoh, &
     &            coh,phihar)

        implicit none
	integer,parameter :: double = selected_real_kind(P=10)
	real(kind=double)tbz(250,721)
	real(kind=double)d2bz(250,721)
	real(kind=double)dbzx(250,721)
	real(kind=double)dbzy(250,721)
	real(kind=double)phihar(10),coh(10,10)
	real(kind=double)rmin,rmax,pi,dr,dteta,phih1,phih2,r
	real(kind=double)x,y,r2,co1,si1,s,pb,px,py,p2
	real(kind=double)r2pk1,r2pk2,co2,si2,s2,s3,teta
	integer ncoh(10),nr,ntet,i,j,k,k1,k2

	pi	= acos(-1._double)
	dr	= (rmax-rmin)/(nr-1._double)
	dteta	= 360._double/(ntet-1._double)*pi/180._double
	phih1	= phihar(1)*pi/180._double
	phih2	= phihar(2)*pi/180._double

	if(ncoh(1).eq.0)goto 200
	do i=1,nr
	  r	= rmin + (i-1._double)*dr
	  do j=1,ntet-1
	    teta	= (j-1._double)*dteta
	    x		= r*cos(teta)
	    y		= r*sin(teta)
	    r2		= x**2+y**2
	    co1		= cos(-phih1)
	    si1		= sin(-phih1)
	    s		= x*co1 - y*si1
	    pb	=  coh(1,1)*s
	    px	=  coh(1,1)*co1
	    py	= -coh(1,1)*si1
  	    p2	=  0.
	    do k=2,ncoh(1)
	      k1	= k - 1
	      k2	= k - 2
	      r2pk1	= r2**k1
	      if(k.eq.2 .and. r2.le.0._double)then
	        r2pk2	= 1._double
	      else
	        r2pk2	= r2**k2
	      endif
	      pb  = pb + coh(k,1)*r2pk1*s
	      px  = px + coh(k,1)*r2pk2*(2._double*x*k1*s + r2*co1)
	      py  = py + coh(k,1)*r2pk2*(2._double*y*k1*s - r2*si1)
! 3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
	      p2	= p2 + coh(k,1)*r2pk2*4._double*s*k1*k
	    enddo
	    tbz(i,j)	= tbz(i,j)  + pb
	    dbzx(i,j)	= dbzx(i,j) + px
	    dbzy(i,j)	= dbzy(i,j) + py
	    d2bz(i,j)	= d2bz(i,j) + p2
	  enddo
	enddo

200	if(ncoh(2).eq.0)goto 300
        do i=1,nr
          r     = rmin + (i-1._double)*dr
          do j=1,ntet-1
            teta        = (j-1._double)*dteta
            x           = r*cos(teta)
            y           = r*sin(teta)
            r2          = x**2+y**2
            co2         = cos(-phih2)
            si2         = sin(-phih2)
            s           = (x**2-y**2)*co2 - 2._double*x*y*si2
	    s2		= 2._double*x*co2 - 2._double*y*si2
	    s3		= -2._double*y*co2 - 2._double*x*si2
            pb  =  coh(1,2)*s
            px  =  coh(1,2)*s2
            py  = -coh(1,2)*s3
            p2  =  0._double
            do k=2,ncoh(2)
              k1        = k - 1
              k2        = k - 2
	      r2pk1	= r2**k1
	      if(k.eq.2 .and. r2.le.0._double)then
	        r2pk2	= 1._double
	      else
	        r2pk2	= r2**k2
	      endif
              pb   = pb + coh(k,2)*r2pk1*s
              px   = px + coh(k,2)*r2pk2*(2._double*x*k1*s + r2*s2)
              py   = py + coh(k,2)*r2pk2*(2._double*y*k1*s + r2*s3)
              p2   = p2 + coh(k,2)*r2pk2*4._double*s*k1*(k+1)
            enddo
            tbz(i,j)    = tbz(i,j)  + pb
            dbzx(i,j)   = dbzx(i,j) + px
            dbzy(i,j)   = dbzy(i,j) + py
            d2bz(i,j)   = d2bz(i,j) + p2
          enddo
        enddo

300	continue
	tbz(:,ntet)	= tbz(:,1)
	dbzx(:,ntet)	= dbzx(:,1)
	dbzy(:,ntet)	= dbzy(:,1)
	d2bz(:,ntet)	= d2bz(:,1)

	end

	subroutine piege

	use lion_m
	implicit none
	integer k
	real(kind=double)rpie2

	if(nbrpie.eq.0)return
    where(z.lt.zmin  .or.  z.gt.zmax )
	      ok	= .false.
	end where
	if(ite.le.itepiem)then
	  do k=1,nbrpie
	    rpie2	= rpie(k)**2
	    where(((x-xpie(k))**2+(y-ypie(k))**2).le.rpie2)
	      ok	= .false.
	    end where
	  enddo
	endif
	nbrok	= count(ok)
        if(iflpie.eq.1)then
          where (.not. ok)
           x  = x(1)
           y  = y(1)
           z  = z(1)
           px = px(1)
           py = py(1)
           pz = pz(1)
           ga = ga(1)
          end where
        endif

	end

	subroutine pousse

        use lion_m
	    implicit none

	    if(iflpie.eq.0)then

          where (ok)

            vec  = qdts2m0/ga
            a1   = vec*bx
            a2   = vec*by
            a3   = vec*bz
            vec1 = a1*a1
            vec2 = a2*a2
            vec3 = a3*a3
            vec  = 1._double/(1._double + vec1 + vec2 + vec3)
            s12  =    a1*a2
            s21  =    s12 - a3
            s12  =    s12 + a3
            s13  =    a1*a3
            s31  =    s13 + a2
            s13  =    s13 - a2
            s23  =    a2*a3
            s32  =    s23 - a1
            s23  =    s23 + a1

            fx   = qdtsm0c*ex
            fy   = qdtsm0c*ey
            fz   = qdtsm0c*ez


            vec4 = vec*(fx + 2._double*px)
            vec5 = vec*(fy + 2._double*py)
            vec6 = vec*(fz + 2._double*pz)
            px   = px - (vec2+vec3)*vec4 + s12*vec5 + s13*vec6 + fx
            py   = py + s21*vec4 - (vec1+vec3)*vec5 + s23*vec6 + fy
            pz   = pz + s31*vec4 + s32*vec5 - (vec1+vec2)*vec6 + fz


            ga  = sqrt(1._double + px*px + py*py + pz*pz)
            vec = cdt/ga
            x   = x + vec*px
            y   = y + vec*py
            z   = z + vec*pz

          end where

	    else

            vec  = qdts2m0/ga
            a1   = vec*bx
            a2   = vec*by
            a3   = vec*bz
            vec1 = a1*a1
            vec2 = a2*a2
            vec3 = a3*a3
            vec  = 1._double/(1._double + vec1 + vec2 + vec3)
            s12  =    a1*a2
            s21  =    s12 - a3
            s12  =    s12 + a3
            s13  =    a1*a3
            s31  =    s13 + a2
            s13  =    s13 - a2
            s23  =    a2*a3
            s32  =    s23 - a1
            s23  =    s23 + a1

            fx   = qdtsm0c*ex
            fy   = qdtsm0c*ey
            fz   = qdtsm0c*ez

            vec4 = vec*(fx + 2._double*px)
            vec5 = vec*(fy + 2._double*py)
            vec6 = vec*(fz + 2._double*pz)
            px   = px - (vec2+vec3)*vec4 + s12*vec5 + s13*vec6 + fx
            py   = py + s21*vec4 - (vec1+vec3)*vec5 + s23*vec6 + fy
            pz   = pz + s31*vec4 + s32*vec5 - (vec1+vec2)*vec6 + fz

            ga  = sqrt(1._double + px*px + py*py + pz*pz)
            vec = cdt/ga
            x   = x + vec*px
            y   = y + vec*py
            z   = z + vec*pz

	    endif

	end

        subroutine calcbe

        use lion_m
	    implicit none
	    integer nbiz,ik,ikp1,jk,jkp1,nto,kg,i,k
	    integer ikp2,ikm1,jkp2,jkm1
	    real(kind=double)pcright,pcleft
	    real(kind=double)rap,vk,omhft,pc,dc1,ds1

!       calcul de b

        nbiz    = count( (x/r .gt.1._double) .or. (x/r .lt.-1._double) )
        if (nbiz.ne.0)then
          write(6,*)'nbiz : ',nbiz
          do i=1,np
            rap = x(i)/r(i)
            if ( (rap .gt. 1._double) .or. (rap .lt. -1._double) )then
             write(6,*)i,x(i),y(i),r(i),rap
             write(6,4444)x(i),y(i),r(i),rap
4444         format(4(1x,z16))
            endif
          enddo
        endif

        veca    = x/r
        where ((veca.gt.-1._double).and.(veca.lt.1._double))
         teta = acos(veca)
        end where
        where (veca.ge.1._double)
         teta = acos(1._double)
        end where
        where (veca.le.-1._double)
         teta = acos(-1._double)
        end where
        where(y.lt.0._double)
          teta = 2._double*pi - teta
        end where

        ivec 	= 1
        jvec 	= 1
        where (box)
          ivec          =  (r-rmin)*hrm1 + 1._double
          jvec          =  teta*htetm1   + 1._double
        end where
        veca          =  r - rmin - (ivec-1._double)*hr
        vecb          =  teta - (jvec-1._double)*htet
        vecd1         =  veca-hr
        vecd2         =  vecb-htet
        vecz          =  z*z*.5_double
        do k=1,np
          bz(k)       =  (vecd1(k)*vecd2(k)*(tbz(ivec(k),jvec(k)) -   &
     &                    vecz(k)*d2bz(ivec(k),jvec(k))) -            &
     &                    veca(k)*vecd2(k)*(tbz(ivec(k)+1,jvec(k)) -  &
     &                    vecz(k)*d2bz(ivec(k)+1,jvec(k))) +          &
     &                    veca(k)*vecb(k)*(tbz(ivec(k)+1,jvec(k)+1) - &
     &                    vecz(k)*d2bz(ivec(k)+1,jvec(k)+1))  -       &
     &                    vecd1(k)*vecb(k)*(tbz(ivec(k),jvec(k)+1) -  &
     &                    vecz(k)*d2bz(ivec(k),jvec(k)+1)))           &
     &                    *hrhtetm1
          bx(k)       =   ( vecd1(k)*vecd2(k)*dbzx(ivec(k),jvec(k)) -  &
     &                      veca(k)*vecd2(k)*dbzx(ivec(k)+1,jvec(k))+  &
     &                      veca(k)*vecb(k)*dbzx(ivec(k)+1,jvec(k)+1) -&
     &                      vecd1(k)*vecb(k)*dbzx(ivec(k),jvec(k)+1) ) &
     &                     *z(k)*hrhtetm1
          by(k)       =   ( vecd1(k)*vecd2(k)*dbzy(ivec(k),jvec(k)) -  &
     &                      veca(k)*vecd2(k)*dbzy(ivec(k)+1,jvec(k))+  &
     &                      veca(k)*vecb(k)*dbzy(ivec(k)+1,jvec(k)+1) -&
     &                      vecd1(k)*vecb(k)*dbzy(ivec(k),jvec(k)+1) ) &
     &                     *z(k)*hrhtetm1
        end do
        where ((box .and. bz.eq.0._double).or.(.not.box))
          bx    = 0._double
          by    = 0._double
          bz    = -.00001_double
        end where
        where (r.le.rmin)
          bx    = 0._double
          by    = 0._double
          bz    = tbz(1,1)
        end where

!       calcul de e au centre

        ex              =  0._double
        ey              =  0._double
        ez              =  0._double
        call box1
        if(nbox1.eq.0)goto 100
        where( centre )
          ivec          =  (x-xmin1)*hx1m1 + 1._double
          jvec          =  (y-ymin1)*hy1m1 + 1._double
          veca          =  x-xmin1 - (ivec-1._double)*hx1
          vecb          =  y-ymin1 - (jvec-1._double)*hy1
          vecd1         =  veca-hx1
          vecd2         =  vecb-hy1
          vecp1         =  vecd1 * vecd2
          vecp2         = -veca  * vecd2
          vecp3         =  veca  * vecb
          vecp4         = -vecd1 * vecb
          vecz          =  z*z*.5_double
        end where
        do k=1,np
          if(centre(k))then
            ik          =  ivec(k)
            ikp1        =  ik + 1
            ikp2        =  ik + 2
            ikm1        =  ik - 1
            jk          =  jvec(k)
            jkp1        =  jk + 1
            jkp2        =  jk + 2
            jkm1        =  jk - 1
            vecb1(k)    =  d2pot1(ik,jk)
            vecb2(k)    =  d2pot1(ikp1,jk)
            vecb3(k)    =  d2pot1(ikp1,jkp1)
            vecb4(k)    =  d2pot1(ik,jkp1)
            vex1(k)     =  (-(pot1(ikp1,jk)-pot1(ikm1,jk)) + vecz(k)  &
     &                      *(d2pot1(ikp1,jk)-d2pot1(ikm1,jk)) )      &
     &                      *.5_double/hx1
            vex2(k)     =  (-(pot1(ikp2,jk)-pot1(ik,jk)) + vecz(k)    &
     &                      *(d2pot1(ikp2,jk)-d2pot1(ik,jk))   )      &
     &                      *.5_double/hx1
            vex3(k)     =  (-(pot1(ikp2,jkp1)-pot1(ik,jkp1)) + vecz(k)&
     &                      *(d2pot1(ikp2,jkp1)-d2pot1(ik,jkp1)) )    &
     &                      *.5_double/hx1
            vex4(k)     =  (-(pot1(ikp1,jkp1)-pot1(ikm1,jkp1))+vecz(k)&
     &                      *(d2pot1(ikp1,jkp1)-d2pot1(ikm1,jkp1)) )  &
     &                      *.5_double/hx1
            vey1(k)     =  (-(pot1(ik,jkp1)-pot1(ik,jkm1)) + vecz(k)  &
     &                      *(d2pot1(ik,jkp1)-d2pot1(ik,jkm1)) )      &
     &                      *.5_double/hy1
            vey2(k)     =  (-(pot1(ikp1,jkp1)-pot1(ikp1,jkm1))+vecz(k)&
     &                      *(d2pot1(ikp1,jkp1)-d2pot1(ikp1,jkm1))  ) &
     &                      *.5_double/hy1
            vey3(k)     =  (-(pot1(ikp1,jkp2)-pot1(ikp1,jk)) + vecz(k)&
     &                      *(d2pot1(ikp1,jkp2)-d2pot1(ikp1,jk)) )    &
     &                      *.5_double/hy1
            vey4(k)     =  (-(pot1(ik,jkp2)-pot1(ik,jk)) + vecz(k)    &
     &                      *(d2pot1(ik,jkp2)-d2pot1(ik,jk)) )  &
     &                      *.5_double/hy1
          end if
        end do
        nto             =  omhf*t/2._double/pi
        omhft           =  omhf*t - nto*2._double*pi
        pcright	 	= pothfd*cos(omhft+dphid)
        pcleft		= pothfg*cos(omhft+dphig)
        where(x.ge.0._double)
          vecpc     =  pcright
        end where
        where(x.lt.0._double)
          vecpc     =  pcleft
        end where
        where( centre )
          ex   =  vecpc*(vex1*vecp1+vex2*vecp2+vex3*vecp3+ &
     &                         vex4*vecp4)*hx1hy1m1
          ey   =  vecpc*(vey1*vecp1+vey2*vecp2+vey3*vecp3+ &
     &                         vey4*vecp4)*hx1hy1m1
          ez   =   z*vecpc*( vecp1*vecb1 + vecp2*vecb2         &
     &                           + vecp3*vecb3 + vecp4*vecb4 )*hx1hy1m1
        end where

!       calcul de e peripherique

100     if((nbox.eq.nbox1) .and. (ndm.eq.0) .and.                       &
     &     (ndf.eq.0) .and. (ngex.eq.1))goto 300
        if(nbox.eq.nbox1)return
        dc1             =  cos(teta1)
        ds1             =  sin(teta1)
        where(.not.centre)
          vecx          =   abs(x)*dc1 + (abs(y)-yomega)*ds1
          vecy          =  -abs(x)*ds1 + (abs(y)-yomega)*dc1
        end where
        call box2
          nto           =  omhf*t/2._double/pi
          omhft         =  omhf*t - nto*2._double*pi
        if(nbox2.eq.0)goto 200
        if(nbox1.eq.0)then
          nto           =  omhf*t/2._double/pi
          omhft         =  omhf*t - nto*2._double*pi
          pcright	=  pothfd*cos(omhft+dphid)
          pcleft	=  pothfg*cos(omhft+dphig)
          where(x.ge.0)
            vecpc         =  pcright
          end where
          where(x.lt.0)
            vecpc         =  pcleft
          end where
        end if
        where( periph )
          ivec          =  (vecx-xmin2)*hx2m1 + 1._double
          jvec          =  (vecy-ymin2)*hy2m1 + 1._double
          veca          =  vecx-xmin2 - (ivec-1._double)*hx2
          vecb          =  vecy-ymin2 - (jvec-1._double)*hy2
          vecd1         =  veca-hx2
          vecd2         =  vecb-hy2
          vecp1         =  vecd1 * vecd2
          vecp2         = -veca  * vecd2
          vecp3         =  veca  * vecb
          vecp4         = -vecd1 * vecb
          vecz          =  z*z*.5_double
        end where
        do k=1,np
          if(periph(k))then
            vec1(k)     =  pot2(ivec(k),jvec(k))
            vec2(k)     =  pot2(ivec(k)+1,jvec(k))
            vec3(k)     =  pot2(ivec(k)+1,jvec(k)+1)
            vec4(k)     =  pot2(ivec(k),jvec(k)+1)
            vecb1(k)    =  d2pot2(ivec(k),jvec(k))
            vecb2(k)    =  d2pot2(ivec(k)+1,jvec(k))
            vecb3(k)    =  d2pot2(ivec(k)+1,jvec(k)+1)
            vecb4(k)    =  d2pot2(ivec(k),jvec(k)+1)
          end if
        end do
        where( periph )
          ex            =  -vecpc*(vecd2*(vec1-vec2-vecz*(vecb1-vecb2)) &
     &                     +vecb*(vec3-vec4+vecz*(vecb3-vecb4)))        &
     &                     *hx2hy2m1
          ey            =  -vecpc*(vecd1*(vec1-vec4-vecz*(vecb1-vecb4)) &
     &                     -veca*(vec2-vec3+vecz*(vecb2-vecb3)))        &
     &                     *hx2hy2m1
          vec           =   ex*dc1 - ey*ds1
          ey            =   ex*ds1 + ey*dc1
          ex            =   vec
          ez            =   z*vecpc*( vecp1*vecb1 + vecp2*vecb2         &
     &                      + vecp3*vecb3 + vecp4*vecb4 )               &
     &                      *hx2hy2m1
        end where
        where(periph .and. x.lt.0._double .and. y.gt.0._double)
          ex            =   heps*ex
          ey            =  -heps*ey
        end where
        where(periph .and. x.gt.0._double .and. y.lt.0._double)
          ey            =  -ey
        end where
        where(periph .and. x.lt.0._double .and. y.lt.0._double)
          ex            =   heps*ex
          ey            =   heps*ey
        end where
        where(periph .and. x.lt.0._double)
          ez            =  -heps*ez
        end where

!	ajout test tri vertical

200	continue
	where(r.ge.rmin_ver .and. r.le.rmax_ver .and. &
     &        teta.ge.teta1_ver .and. teta.le.teta2_ver )
	  ez   = ez - pot_ver*(1.-dpotspot)/h_ver/2._double*cos(2.*omhf*t*(1.+anuz/2./harm) + &
     &           2.*phi_ver-pi/2.) &
     &              - pot_ver*(1.+dpotspot)/h_ver/2._double*cos(2.*omhf*t*(1.-anuz/2./harm) + &
     &           2.*phi_ver-pi/2.)
	end where
	pot_ver_vrai = 0.
	if(r(1).ge.rmin_ver .and. r(1).le.rmax_ver .and. &
     &        teta(1).ge.teta1_ver .and. teta(1).le.teta2_ver )pot_ver_vrai = -ez(1)*h_ver

!	fin ajout test tri vertical

!       calcul de l'extraction

        do i=1,ndf
          call boxdf(i)
          if(nextra.ne.0)then
            where( extra )
              vec       = sqrt(px**2+py**2)
              ex        = edf(i)*py/vec
              ey        = -edf(i)*px/vec
            end where
          endif
        enddo
        do i=1,ndm
          call boxdm(i)
          if(nextra.ne.0)then
            where(extra)
              vec = sqrt((x-xdm(i))**2+(y-ydm(i))**2)
              vec1 = vec-rdm(i)
              bx = bx + z*bdm(i)*(x-xdm(i))/vec
              by = by + z*bdm(i)*(y-ydm(i))/vec
              bz = bz + vec1*bdm(i)
            end where
          endif
        enddo

!       mise a jour pour grilles d'extraction ou grille sonde

300     if(rex.lt.0._double .and. ngex.gt.0)then
          if(ite.lt.itegex)goto 888
          do kg=1,ngex
            call boxgex(kg)
            where(franchi)indgex = kg
            do i=1,n
              if(franchi(i))then
                write(70+kg,2777)i,ite,x(i),y(i),z(i),px(i),py(i),      &
     &                           pz(i),bz(i)
2777            format(1x,i5,1x,i6,6(1x,e15.8),1x,f8.5)
              endif
            enddo
          enddo
        endif
888     continue

        end
