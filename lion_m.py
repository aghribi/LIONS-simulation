"""
lion_m.py - Module de données pour LIONS (converti depuis lion_m.f90)

Programme LIONS : Simulation de dynamique d'ions dans un cyclotron
Auteur original : Patrick Bertrand - GANIL
Conversion Python : 2026
"""

import numpy as np


# Constantes de dimensionnement
NF = 100
LX1 = 500
LY1 = 300
LX2 = 50
LY2 = 1000
NDFM = 10
NDMM = 10


class LionData:
    """Classe contenant toutes les données du module lion_m."""

    def __init__(self):
        # --- Scalaires entiers ---
        self.n = 0
        self.nbrok = 0
        self.nbrpie = 0
        self.itepiem = 0
        self.iflpie = 0

        # --- Chaînes de caractères ---
        self.corps = ''
        self.ficb = ''
        self.ficin = ''
        self.ficz = ''
        self.ficxy = ''
        self.ficv1 = ''
        self.ficv2 = ''
        self.ficsig = ''
        self.ficcor = ''
        self.ficmat = ''
        self.ficdc = ''
        self.ficloa = 'no'
        self.ficsto = 'no'
        self.ficinf = ''
        self.ficdebgex = ''
        self.ficfingex = ''
        self.ficgex = [''] * 10
        self.ficpie = ''
        self.espace = ''
        self.baryc = ''
        self.reponse = ''
        self.repswap = ''

        # --- Tableaux de taille NF ---
        self.nbrion = np.zeros(NF, dtype=np.int32)
        self.index = np.zeros(NF, dtype=np.int32)
        self.amat = np.zeros(NF, dtype=np.float64)
        self.a = np.zeros(NF, dtype=np.float64)
        self.qn = np.zeros(NF, dtype=np.float64)
        self.r0 = np.zeros(NF, dtype=np.float64)
        self.teta0 = np.zeros(NF, dtype=np.float64)
        self.z0 = np.zeros(NF, dtype=np.float64)
        self.pr0 = np.zeros(NF, dtype=np.float64)
        self.pz0 = np.zeros(NF, dtype=np.float64)
        self.w0 = np.zeros(NF, dtype=np.float64)
        self.phi0 = np.zeros(NF, dtype=np.float64)

        self.rdeb1 = 0.0
        self.dbdeb1 = 0.0

        # --- Tableaux allocatables (initialisés dans allo()) ---
        self.x = None
        self.y = None
        self.z = None
        self.px = None
        self.py = None
        self.pz = None
        self.ga = None
        self.p = None
        self.r = None
        self.teta = None
        self.am0 = None
        self.q = None
        self.am0c = None
        self.am0c2 = None
        self.am02c2 = None
        self.bx = None
        self.by = None
        self.bz = None
        self.ex = None
        self.ey = None
        self.ez = None
        self.fx = None
        self.fy = None
        self.fz = None
        self.a1 = None
        self.a2 = None
        self.a3 = None
        self.s12 = None
        self.s13 = None
        self.s21 = None
        self.s23 = None
        self.s31 = None
        self.s32 = None
        self.qdts2m0 = None
        self.qdtsm0c = None
        self.vec = None
        self.vec1 = None
        self.vec2 = None
        self.vec3 = None
        self.vec4 = None
        self.vec5 = None
        self.vec6 = None
        self.vecp1 = None
        self.vecp2 = None
        self.vecp3 = None
        self.vecp4 = None
        self.vecb1 = None
        self.vecb2 = None
        self.vecb3 = None
        self.vecb4 = None
        self.vecd1 = None
        self.vecd2 = None
        self.vecd3 = None
        self.vecd4 = None
        self.vecx = None
        self.vecy = None
        self.vecz = None
        self.vecx1 = None
        self.vecx2 = None
        self.vecx3 = None
        self.vecx4 = None
        self.vecy1 = None
        self.vecy2 = None
        self.vecy3 = None
        self.vecy4 = None
        self.veca = None
        self.vecb_arr = None  # renamed to avoid conflict with vecb
        self.vecpc = None
        self.vex1 = None
        self.vex2 = None
        self.vex3 = None
        self.vex4 = None
        self.vey1 = None
        self.vey2 = None
        self.vey3 = None
        self.vey4 = None
        self.xgri = None
        self.ygri = None
        self.rgri = None
        self.tprgri = None
        self.phigri = None
        self.dxs = None
        self.dxd = None
        self.dxtmp = None
        self.dxsgri = None
        self.box = None
        self.centre = None
        self.periph = None
        self.extra = None
        self.fenlog = None
        self.sortant = None
        self.franchi = None
        self.ok = None
        self.ivec = None
        self.jvec = None
        self.iflag = None
        self.indgex = None

        # --- Tableaux de taille fixe (cartes de champ et potentiel) ---
        self.pot1 = np.zeros((LX1, LY1), dtype=np.float64)
        self.d2pot1 = np.zeros((LX1, LY1), dtype=np.float64)
        self.pot2 = np.zeros((LX2, LY2), dtype=np.float64)
        self.d2pot2 = np.zeros((LX2, LY2), dtype=np.float64)

        # Déflecteurs
        self.rdf = np.zeros(NDFM, dtype=np.float64)
        self.edf = np.zeros(NDFM, dtype=np.float64)
        self.xdfi = np.zeros(NDFM, dtype=np.float64)
        self.xdfo = np.zeros(NDFM, dtype=np.float64)
        self.ydfi = np.zeros(NDFM, dtype=np.float64)
        self.ydfo = np.zeros(NDFM, dtype=np.float64)
        self.tetdfi = np.zeros(NDFM, dtype=np.float64)
        self.tetdfo = np.zeros(NDFM, dtype=np.float64)
        self.uxdfi = np.zeros(NDFM, dtype=np.float64)
        self.uydfi = np.zeros(NDFM, dtype=np.float64)
        self.uxdfo = np.zeros(NDFM, dtype=np.float64)
        self.uydfo = np.zeros(NDFM, dtype=np.float64)

        # Dipôles magnétiques
        self.rdm = np.zeros(NDMM, dtype=np.float64)
        self.bdm = np.zeros(NDMM, dtype=np.float64)
        self.xdm = np.zeros(NDMM, dtype=np.float64)
        self.ydm = np.zeros(NDMM, dtype=np.float64)
        self.tetdmi = np.zeros(NDMM, dtype=np.float64)
        self.tetdmo = np.zeros(NDMM, dtype=np.float64)
        self.uxdmi = np.zeros(NDMM, dtype=np.float64)
        self.uydmi = np.zeros(NDMM, dtype=np.float64)
        self.uxdmo = np.zeros(NDMM, dtype=np.float64)
        self.uydmo = np.zeros(NDMM, dtype=np.float64)

        # Grilles d'extraction
        self.tetagex = np.zeros(10, dtype=np.float64)
        self.xgex = np.zeros(10, dtype=np.float64)
        self.ygex = np.zeros(10, dtype=np.float64)
        self.xgexc = np.zeros(10, dtype=np.float64)
        self.ygexc = np.zeros(10, dtype=np.float64)
        self.drgex = np.zeros(10, dtype=np.float64)
        self.uxgex = np.zeros(10, dtype=np.float64)
        self.uygex = np.zeros(10, dtype=np.float64)

        # Pièges à ions
        self.xpie = np.zeros(100, dtype=np.float64)
        self.ypie = np.zeros(100, dtype=np.float64)
        self.rpie = np.zeros(100, dtype=np.float64)

        # --- Scalaires divers ---
        self.idebug = 0
        self.ifirst = 0
        self.nbrfam = 0
        self.w0th = 0.0
        self.w0mevth = 0.0
        self.ga0th = 0.0
        self.r0th = 0.0
        self.b0th = 0.0
        self.w1th = 0.0
        self.w1mevth = 0.0
        self.ga1th = 0.0
        self.r1th = 0.0
        self.b1th = 0.0
        self.w0mev1 = 0.0

        self.np_ = 0  # 'np' est réservé pour numpy

        self.iflcor = 0
        self.fendwsw = 0.0
        self.cor = np.zeros((6, 6), dtype=np.float64)
        self.rcor = np.zeros((6, 6), dtype=np.float64)
        self.trcor = np.zeros((6, 6), dtype=np.float64)
        self.rot = np.zeros((6, 6), dtype=np.float64)
        self.t1 = np.zeros((6, 6), dtype=np.float64)
        self.cordc = np.zeros((6, 6), dtype=np.float64)
        self.corgri = np.zeros((6, 6), dtype=np.float64)
        self.corfen = np.zeros((6, 6), dtype=np.float64)
        self.corfeng = np.zeros((6, 6), dtype=np.float64)

        self.nbox = 0
        self.nbox1 = 0
        self.nbox2 = 0
        self.nextra = 0
        self.nfen = 0
        self.zmin = 0.0
        self.zmax = 0.0
        self.xmin1 = 0.0
        self.xmax1 = 0.0
        self.ymin1 = 0.0
        self.ymax1 = 0.0
        self.rmax1 = 0.0
        self.hx1 = 0.0
        self.hy1 = 0.0
        self.hx1hy1 = 0.0
        self.hx1m1 = 0.0
        self.hy1m1 = 0.0
        self.hx1hy1m1 = 0.0
        self.xmin2 = 0.0
        self.xmax2 = 0.0
        self.ymin2 = 0.0
        self.ymax2 = 0.0
        self.hx2 = 0.0
        self.hy2 = 0.0
        self.hx2hy2 = 0.0
        self.hx2m1 = 0.0
        self.hy2m1 = 0.0
        self.hx2hy2m1 = 0.0

        self.dbsb = 0.0
        self.b_factor = None  # Facteur B (T) : max|Bz| impose a la carte ; None = carte telle quelle
        self.t0 = 0.0
        self.t = 0.0
        self.dt = 0.0
        self.cdt = 0.0
        self.clum = 0.0
        self.pi = 0.0
        self.rac8 = 0.0

        self.nr = 0
        self.ntet = 0
        self.rmin = 0.0
        self.rmax = 0.0
        self.hr = 0.0
        self.htet = 0.0
        self.hrhtet = 0.0
        self.hrm1 = 0.0
        self.htetm1 = 0.0
        self.hrhtetm1 = 0.0
        self.tbz = np.zeros((250, 721), dtype=np.float64)
        self.d2bz = np.zeros((250, 721), dtype=np.float64)
        self.dbzx = np.zeros((250, 721), dtype=np.float64)
        self.dbzy = np.zeros((250, 721), dtype=np.float64)
        self.dbzxlu = np.zeros((250, 91), dtype=np.float64)
        self.dbzylu = np.zeros((250, 91), dtype=np.float64)

        self.fhf = 0.0
        self.omhf = 0.0
        self.omr = 0.0
        self.pothf = 0.0
        self.psi = np.zeros((100, 2), dtype=np.float64)
        self.dpsi = 0.0
        self.teta1 = 0.0
        self.harm = 0.0
        self.heps = 0.0
        self.yomega = 0.0

        self.iteext = 0
        self.rex = 0.0
        self.tet1ex = 0.0
        self.tet2ex = 0.0
        self.eex = 0.0
        self.ndf = 0
        self.ndm = 0

        self.nbdt = 0
        self.ite = 0
        self.ntours = 0
        self.npaspt = 0
        self.numtou = 0
        self.tetold = 0.0
        self.tetnew = 0.0
        self.tetall = 0.0
        self.ifreq = 0
        self.numbid = 0
        self.nbrbid = 0
        self.ifreb = 0
        self.itebid = 0

        self.ax0 = np.zeros(10, dtype=np.float64)
        self.ax1 = np.zeros(10, dtype=np.float64)
        self.ay0 = np.zeros(10, dtype=np.float64)
        self.ay1 = np.zeros(10, dtype=np.float64)
        self.xgrav = np.zeros(10800, dtype=np.float64)
        self.ygrav = np.zeros(10800, dtype=np.float64)
        self.ngravm = 0
        self.ngrav = 0

        self.xf = 0.0
        self.yf = 0.0
        self.rf = 0.0
        self.tetf = 0.0
        self.phif = 0.0
        self.pxf = 0.0
        self.pyf = 0.0
        self.pf = 0.0
        self.prf = 0.0
        self.tprf = 0.0
        self.wf = 0.0
        self.zf = 0.0
        self.tpzf = 0.0

        self.rpcmax = 0.0
        self.dphid = 0.0
        self.dphig = 0.0
        self.pothfd = 0.0
        self.pothfg = 0.0

        self.ngex = 0
        self.itegex = 0

        self.ncoh = np.zeros(10, dtype=np.int32)
        self.phihar = np.zeros(10, dtype=np.float64)
        self.coh = np.zeros((10, 10), dtype=np.float64)
        self.rbib1 = 0.0
        self.rbib2 = 0.0
        self.tetabib1 = 0.0
        self.tetabib2 = 0.0
        self.rbibtest = 0.0
        self.dbbib = 0.0
        self.xrec = 0.0
        self.yrec = 0.0
        self.courec = 0.0
        self.rinflu = 0.0
        self.covg = 0.0
        self.covd = 0.0

        self.rmin_ver = 0.0
        self.rmax_ver = 0.0
        self.teta1_ver = 0.0
        self.teta2_ver = 0.0
        self.pot_ver = 0.0
        self.h_ver = 0.0
        self.phi_ver = 0.0
        self.anuz = 0.0
        self.pot_ver_vrai = 0.0
        self.dpotspot = 0.0
        self.coef_ver = 0.0
        self.coef_cos = 0.0
        self.anuz_cyclo = 0.0
        self.coupnuz = 0.0

    def allo(self):
        """Allocation des tableaux dynamiques (équivalent du subroutine allo)."""
        n = self.n
        self.x = np.zeros(n, dtype=np.float64)
        self.y = np.zeros(n, dtype=np.float64)
        self.z = np.zeros(n, dtype=np.float64)
        self.px = np.zeros(n, dtype=np.float64)
        self.py = np.zeros(n, dtype=np.float64)
        self.pz = np.zeros(n, dtype=np.float64)
        self.ga = np.zeros(n, dtype=np.float64)
        self.p = np.zeros(n, dtype=np.float64)
        self.r = np.zeros(n, dtype=np.float64)
        self.teta = np.zeros(n, dtype=np.float64)
        self.am0 = np.zeros(n, dtype=np.float64)
        self.q = np.zeros(n, dtype=np.float64)
        self.am0c = np.zeros(n, dtype=np.float64)
        self.am0c2 = np.zeros(n, dtype=np.float64)
        self.am02c2 = np.zeros(n, dtype=np.float64)
        self.bx = np.zeros(n, dtype=np.float64)
        self.by = np.zeros(n, dtype=np.float64)
        self.bz = np.zeros(n, dtype=np.float64)
        self.ex = np.zeros(n, dtype=np.float64)
        self.ey = np.zeros(n, dtype=np.float64)
        self.ez = np.zeros(n, dtype=np.float64)
        self.fx = np.zeros(n, dtype=np.float64)
        self.fy = np.zeros(n, dtype=np.float64)
        self.fz = np.zeros(n, dtype=np.float64)
        self.a1 = np.zeros(n, dtype=np.float64)
        self.a2 = np.zeros(n, dtype=np.float64)
        self.a3 = np.zeros(n, dtype=np.float64)
        self.s12 = np.zeros(n, dtype=np.float64)
        self.s13 = np.zeros(n, dtype=np.float64)
        self.s21 = np.zeros(n, dtype=np.float64)
        self.s23 = np.zeros(n, dtype=np.float64)
        self.s31 = np.zeros(n, dtype=np.float64)
        self.s32 = np.zeros(n, dtype=np.float64)
        self.qdts2m0 = np.zeros(n, dtype=np.float64)
        self.qdtsm0c = np.zeros(n, dtype=np.float64)
        self.vec = np.zeros(n, dtype=np.float64)
        self.vec1 = np.zeros(n, dtype=np.float64)
        self.vec2 = np.zeros(n, dtype=np.float64)
        self.vec3 = np.zeros(n, dtype=np.float64)
        self.vec4 = np.zeros(n, dtype=np.float64)
        self.vec5 = np.zeros(n, dtype=np.float64)
        self.vec6 = np.zeros(n, dtype=np.float64)
        self.vecp1 = np.zeros(n, dtype=np.float64)
        self.vecp2 = np.zeros(n, dtype=np.float64)
        self.vecp3 = np.zeros(n, dtype=np.float64)
        self.vecp4 = np.zeros(n, dtype=np.float64)
        self.vecb1 = np.zeros(n, dtype=np.float64)
        self.vecb2 = np.zeros(n, dtype=np.float64)
        self.vecb3 = np.zeros(n, dtype=np.float64)
        self.vecb4 = np.zeros(n, dtype=np.float64)
        self.vecd1 = np.zeros(n, dtype=np.float64)
        self.vecd2 = np.zeros(n, dtype=np.float64)
        self.vecd3 = np.zeros(n, dtype=np.float64)
        self.vecd4 = np.zeros(n, dtype=np.float64)
        self.vecx = np.zeros(n, dtype=np.float64)
        self.vecy = np.zeros(n, dtype=np.float64)
        self.vecz = np.zeros(n, dtype=np.float64)
        self.vecx1 = np.zeros(n, dtype=np.float64)
        self.vecx2 = np.zeros(n, dtype=np.float64)
        self.vecx3 = np.zeros(n, dtype=np.float64)
        self.vecx4 = np.zeros(n, dtype=np.float64)
        self.vecy1 = np.zeros(n, dtype=np.float64)
        self.vecy2 = np.zeros(n, dtype=np.float64)
        self.vecy3 = np.zeros(n, dtype=np.float64)
        self.vecy4 = np.zeros(n, dtype=np.float64)
        self.veca = np.zeros(n, dtype=np.float64)
        self.vecb_arr = np.zeros(n, dtype=np.float64)
        self.vecpc = np.zeros(n, dtype=np.float64)
        self.vex1 = np.zeros(n, dtype=np.float64)
        self.vex2 = np.zeros(n, dtype=np.float64)
        self.vex3 = np.zeros(n, dtype=np.float64)
        self.vex4 = np.zeros(n, dtype=np.float64)
        self.vey1 = np.zeros(n, dtype=np.float64)
        self.vey2 = np.zeros(n, dtype=np.float64)
        self.vey3 = np.zeros(n, dtype=np.float64)
        self.vey4 = np.zeros(n, dtype=np.float64)
        self.xgri = np.zeros(n, dtype=np.float64)
        self.ygri = np.zeros(n, dtype=np.float64)
        self.rgri = np.zeros(n, dtype=np.float64)
        self.tprgri = np.zeros(n, dtype=np.float64)
        self.phigri = np.zeros(n, dtype=np.float64)
        self.dxs = np.zeros((6, n), dtype=np.float64)
        self.dxd = np.zeros((6, n), dtype=np.float64)
        self.dxtmp = np.zeros((6, n), dtype=np.float64)
        self.dxsgri = np.zeros((6, n), dtype=np.float64)
        self.box = np.zeros(n, dtype=bool)
        self.centre = np.zeros(n, dtype=bool)
        self.periph = np.zeros(n, dtype=bool)
        self.extra = np.zeros(n, dtype=bool)
        self.fenlog = np.zeros(n, dtype=bool)
        self.sortant = np.zeros(n, dtype=bool)
        self.franchi = np.zeros(n, dtype=bool)
        self.ok = np.zeros(n, dtype=bool)
        self.ivec = np.zeros(n, dtype=np.int32)
        self.jvec = np.zeros(n, dtype=np.int32)
        self.iflag = np.zeros(n, dtype=np.int32)
        self.indgex = np.zeros(n, dtype=np.int32)

        print('allocations effectuees...')
        print(f'nombre de particules : {n}')
