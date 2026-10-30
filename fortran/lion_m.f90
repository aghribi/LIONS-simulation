module lion_m
! 3456789 123456789 123456789 123456789 123456789 123456789 123456789 12 

        implicit none 

       	integer,parameter :: double = selected_real_kind (P=10)
       	integer,parameter :: nf=100,lx1=500,ly1=300,lx2=50,ly2=1000
    	integer,parameter :: ndfm=10,ndmm=10

	integer n,nbrok,nbrpie,itepiem,iflpie

	character*80 corps
	character*100 ficb
	character*80 ficin,ficz,ficxy,ficv1,ficv2,ficsig
	character*80 ficcor,ficmat,ficdc
	character*80 ficloa,ficsto,ficinf
	character*80 ficdebgex,ficfingex,ficgex(10),ficpie
	character*4  espace
	character*4  baryc
	character*1  reponse,repswap  ! Patrick Bertrand, ajour de repswap, mars 2017

	integer nbrion(nf),index(nf)

	real(kind=double) rdeb1,dbdeb1
	real(kind=double) amat(nf),a(nf),qn(nf)
	real(kind=double) r0(nf),teta0(nf),z0(nf),pr0(nf),pz0(nf),w0(nf),phi0(nf)
        real(kind=double),dimension (:),allocatable :: x,y,z,px,py,pz,ga,p
        real(kind=double),dimension (:),allocatable :: r,teta
	real(kind=double),dimension (:),allocatable :: am0,q
	real(kind=double),dimension (:),allocatable :: am0c,am0c2,am02c2
        real(kind=double),dimension (:),allocatable :: bx,by,bz
        real(kind=double),dimension (:),allocatable :: ex,ey,ez
        real(kind=double),dimension (:),allocatable :: fx,fy,fz
        real(kind=double),dimension (:),allocatable :: a1,a2,a3
        real(kind=double),dimension (:),allocatable :: s12,s13,s21,s23,s31,s32
	real(kind=double),dimension (:),allocatable :: qdts2m0,qdtsm0c
	real(kind=double),dimension (:),allocatable :: vec,vec1,vec2,vec3
	real(kind=double),dimension (:),allocatable :: vec4,vec5,vec6
	real(kind=double),dimension (:),allocatable :: vecp1,vecp2,vecp3,vecp4
	real(kind=double),dimension (:),allocatable :: vecb1,vecb2,vecb3,vecb4
	real(kind=double),dimension (:),allocatable :: vecd1,vecd2,vecd3,vecd4
	real(kind=double),dimension (:),allocatable :: vecx,vecy,vecz
	real(kind=double),dimension (:),allocatable :: vecx1,vecx2,vecx3,vecx4
	real(kind=double),dimension (:),allocatable :: vecy1,vecy2,vecy3,vecy4
	real(kind=double),dimension (:),allocatable :: veca,vecb,vecpc
	real(kind=double),dimension (:),allocatable :: vex1,vex2,vex3,vex4
	real(kind=double),dimension (:),allocatable :: vey1,vey2,vey3,vey4
	real(kind=double),dimension (:),allocatable :: xgri,ygri,rgri,tprgri,phigri
	real(kind=double),dimension (:,:),allocatable :: dxs,dxd,dxtmp,dxsgri
        logical,dimension (:),allocatable :: box,centre,periph
        logical,dimension (:),allocatable :: extra,fenlog,sortant,franchi
        logical,dimension (:),allocatable :: ok
	integer,dimension (:),allocatable :: ivec,jvec,iflag,indgex

	real(kind=double) pot1(lx1,ly1),d2pot1(lx1,ly1)
	real(kind=double) pot2(lx2,ly2),d2pot2(lx2,ly2)
        real(kind=double) rdf(ndfm),edf(ndfm),xdfi(ndfm),xdfo(ndfm)
        real(kind=double) ydfi(ndfm),ydfo(ndfm)
        real(kind=double) tetdfi(ndfm),tetdfo(ndfm)
        real(kind=double) uxdfi(ndfm),uydfi(ndfm),uxdfo(ndfm),uydfo(ndfm)
        real(kind=double) rdm(ndmm),bdm(ndmm),xdm(ndmm),ydm(ndmm)
        real(kind=double) tetdmi(ndmm),tetdmo(ndmm)
        real(kind=double) uxdmi(ndmm),uydmi(ndmm),uxdmo(ndmm),uydmo(ndmm)
	real(kind=double) tetagex(10),xgex(10),ygex(10),xgexc(10),ygexc(10)
	real(kind=double) drgex(10),uxgex(10),uygex(10)
	real(kind=double) xpie(100),ypie(100),rpie(100)

	integer idebug
	integer ifirst

	integer nbrfam
	real(kind=double)w0th,w0mevth,ga0th,r0th,b0th
	real(kind=double)w1th,w1mevth,ga1th,r1th,b1th
	real(kind=double)w0mev1

	integer np

	integer iflcor
	real(kind=double)fendwsw
	real(kind=double)cor(6,6),rcor(6,6),trcor(6,6),rot(6,6)
	real(kind=double)t1(6,6),cordc(6,6),corgri(6,6),corfen(6,6)
	real(kind=double)corfeng(6,6)


	integer nbox,nbox1,nbox2,nextra,nfen
	real(kind=double)zmin,zmax
	real(kind=double)xmin1,xmax1,ymin1,ymax1,rmax1
	real(kind=double)hx1,hy1,hx1hy1
	real(kind=double)hx1m1,hy1m1,hx1hy1m1
	real(kind=double)xmin2,xmax2,ymin2,ymax2
	real(kind=double)hx2,hy2,hx2hy2
	real(kind=double)hx2m1,hy2m1,hx2hy2m1

	real(kind=double)dbsb

	real(kind=double)t0,t,dt,cdt

	real(kind=double)clum,pi,rac8
	
	integer nr,ntet
	real(kind=double)rmin,rmax,hr,htet,hrhtet
	real(kind=double)hrm1,htetm1,hrhtetm1
	real(kind=double)tbz(250,721),d2bz(250,721)
	real(kind=double)dbzx(250,721),dbzy(250,721)
	real(kind=double)dbzxlu(250,91),dbzylu(250,91)
	real(kind=double)fhf,omhf,omr,pothf,psi(100,2),dpsi,teta1
	real(kind=double)harm,heps,yomega
	integer iteext
	real(kind=double)rex,tet1ex,tet2ex,eex
	integer ndf
	integer ndm

	integer nbdt,ite,ntours,npaspt
	integer numtou
	real(kind=double)tetold,tetnew,tetall
	integer ifreq
	integer numbid,nbrbid,ifreb,itebid
	real(kind=double)ax0(10),ax1(10),ay0(10),ay1(10)
	real(kind=double)xgrav(10800),ygrav(10800)
	integer ngravm,ngrav

	real(kind=double)xf,yf,rf,tetf,phif
	real(kind=double)pxf,pyf,pf,prf,tprf
	real(kind=double)wf,zf,tpzf

	real(kind=double)rpcmax

	real(kind=double)dphid,dphig,pothfd,pothfg

	integer ngex,itegex

	integer ncoh(10)
	real(kind=double)phihar(10),coh(10,10)
	real(kind=double)rbib1,rbib2,tetabib1,tetabib2
	real(kind=double)rbibtest,dbbib
	real(kind=double)xrec,yrec,courec,rinflu
	real(kind=double)covg,covd

	real(kind=double)rmin_ver,rmax_ver,teta1_ver,teta2_ver
	real(kind=double)pot_ver,h_ver,phi_ver,anuz,pot_ver_vrai,dpotspot
	real(kind=double)coef_ver,coef_cos,anuz_cyclo,coupnuz

contains

	subroutine allo
         
        implicit none

        allocate(x(n),y(n),z(n),px(n),py(n),pz(n),ga(n),p(n))
        allocate(r(n),teta(n))
        allocate(am0(n),q(n))
        allocate(am0c(n),am0c2(n),am02c2(n))
        allocate(bx(n),by(n),bz(n))
        allocate(ex(n),ey(n),ez(n))
        allocate(fx(n),fy(n),fz(n))
        allocate(a1(n),a2(n),a3(n))
        allocate(s12(n),s13(n),s21(n),s23(n),s31(n),s32(n))
        allocate(qdts2m0(n),qdtsm0c(n))
        allocate(vec(n),vec1(n),vec2(n),vec3(n))
        allocate(vec4(n),vec5(n),vec6(n))
        allocate(vecp1(n),vecp2(n),vecp3(n),vecp4(n))
        allocate(vecb1(n),vecb2(n),vecb3(n),vecb4(n))
        allocate(vecd1(n),vecd2(n),vecd3(n),vecd4(n))
        allocate(vecx(n),vecy(n),vecz(n))
        allocate(vecx1(n),vecx2(n),vecx3(n),vecx4(n))
        allocate(vecy1(n),vecy2(n),vecy3(n),vecy4(n))
        allocate(veca(n),vecb(n),vecpc(n))
        allocate(vex1(n),vex2(n),vex3(n),vex4(n))
        allocate(vey1(n),vey2(n),vey3(n),vey4(n))
        allocate(xgri(n),ygri(n),rgri(n),tprgri(n),phigri(n))
        allocate(dxs(6,n),dxd(6,n),dxtmp(6,n),dxsgri(6,n))
        allocate(box(n),centre(n),periph(n))
        allocate(extra(n),fenlog(n),sortant(n),franchi(n))
        allocate(ivec(n),jvec(n),iflag(n),indgex(n))
        allocate(ok(n))

        write(6,*)'allocations effectuees...'
        write(6,*)'nombre de particules : ',n

	end subroutine allo

end module lion_m
