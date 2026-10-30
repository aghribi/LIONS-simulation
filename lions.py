"""
lions.py - Programme principal LIONS converti de Fortran 90 en Python
Programme de simulation de dynamique d'ions dans un cyclotron
Auteur original : Patrick Bertrand - GANIL
Version 1 : décembre 1998...décembre 2004
Version 2 : mars 2017, modif seed, swap etc...
Conversion Python : 2026
"""

import numpy as np
import struct
import sys
import os
from lion_m import LionData


# =============================================================================
# Classe BidimManager - gestion des histogrammes 2D (bidims)
# =============================================================================
class BidimManager:
    """Gestion des bidims (histogrammes 2D) - remplace les common blocks Fortran."""

    def __init__(self):
        self.nbrbid = 0
        self.n1bid = 0
        self.n2bid = 0
        self.nbxbid = np.zeros(50, dtype=np.int32)
        self.nbybid = np.zeros(50, dtype=np.int32)
        self.lmxbid = np.zeros(50, dtype=np.int32)
        self.lmybid = np.zeros(50, dtype=np.int32)
        self.ind = np.zeros(50, dtype=np.int32)
        self.ibubid = np.zeros((4096, 50), dtype=np.int32)
        self.files = {}

    def openbi(self, nbr, n1, nbx, nby, titre, labelx, labely, bornes):
        """Ouvre les fichiers bidim."""
        self.nbrbid = 0
        idate = [0, 0, 0]
        itime = [0, 0, 0]
        type_str = 'bidim'
        lieu = 'ganil'
        auteur = 'inconnu'
        forma = 'Binaire'
        vide = ' ' * 120

        if nbr <= 0 or nbr > 50:
            print('openbi : nombre de bidims incorrect')
            return 1
        if n1 <= 0:
            print('openbi : numero du premier bidim incorrect')
            return 2

        self.n1bid = n1
        self.n2bid = n1 + nbr - 1

        for i in range(nbr):
            mx = 2**nbx[i] - 1
            if mx > 8191:
                print(f'openbi : pour bidim {i+1}, nombre de bits en x incorrect : {nbx[i]}')
                return 3
            self.nbxbid[i] = nbx[i]
            self.lmxbid[i] = mx
            my = 2**nby[i] - 1
            if my > 8191:
                print(f'openbi : pour bidim {i+1}, nombre de bits en y incorrect : {nby[i]}')
                return 4
            self.nbybid[i] = nby[i]
            self.lmybid[i] = my

        for i in range(nbr):
            num = n1 + i
            filnam = f'{num}.bid'
            print(f'filnam : {filnam}')
            try:
                f = open(filnam, 'wb')
                self.files[40 + i + 1] = f
                maxx = self.lmxbid[i]
                maxy = self.lmybid[i]
                cxmi = f'{bornes[0, 0, i]:.10e}'
                cxma = f'{bornes[1, 0, i]:.10e}'
                cymi = f'{bornes[0, 1, i]:.10e}'
                cyma = f'{bornes[1, 1, i]:.10e}'
                header = struct.pack('i', num)
                header += type_str.encode('ascii').ljust(10)
                header += struct.pack('3i', *idate)
                header += struct.pack('3i', *itime)
                header += lieu.encode('ascii').ljust(10)
                header += auteur.encode('ascii').ljust(10)
                header += titre[i].encode('ascii').ljust(25)
                header += labelx[i].encode('ascii').ljust(25)
                header += labely[i].encode('ascii').ljust(25)
                header += struct.pack('2i', maxx, maxy)
                header += forma.encode('ascii').ljust(10)
                header += cxmi.encode('ascii').ljust(20)
                header += cxma.encode('ascii').ljust(20)
                header += cymi.encode('ascii').ljust(20)
                header += cyma.encode('ascii').ljust(20)
                header += vide.encode('ascii')[:120]
                f.write(header)
                self.ind[i] = 1
            except IOError:
                print(f'openbi : probleme open du bidim : {i+1}')
                for k in range(i):
                    self.files[40 + k + 1].close()
                return 5

        print('openbi : open des bidims sur disque effectue correctement')
        self.nbrbid = nbr
        return 0

    def incbid(self, n, ix, iy):
        """Incrémente un bidim."""
        if n < self.n1bid or n > self.n2bid:
            print(f'incbid : numero bidim incorrect : {n}')
            print(f'...en dehors de : {self.n1bid} {self.n2bid}')
            return 1
        num = n - self.n1bid  # 0-based
        ix = max(0, min(ix, self.lmxbid[num]))
        iy = max(0, min(iy, self.lmybid[num]))
        if self.ind[num] == 4097:
            self._wribid(40 + num + 1, num)
            self.ind[num] = 1
        idx = self.ind[num] - 1  # 0-based
        self.ibubid[idx, num] = ix * 65536 + iy
        self.ind[num] += 1
        return 0

    def clobid(self):
        """Ferme les bidims."""
        for num in range(self.nbrbid):
            if self.ind[num] == 4097:
                self._wribid(40 + num + 1, num)
                self.ind[num] = 1
            elif self.ind[num] != 1:
                idx = self.ind[num] - 1
                self.ibubid[idx, num] = 9999 * 65536 + 9999
                self._wribid(40 + num + 1, num)
            unit = 40 + num + 1
            if unit in self.files:
                self.files[unit].close()
        return 0

    def _wribid(self, iu, num):
        """Écrit un buffer bidim."""
        if iu in self.files:
            data = self.ibubid[:4096, num]
            self.files[iu].write(data.tobytes())


# =============================================================================
# Fonctions utilitaires
# =============================================================================

def read_input(prompt=''):
    """Lit une ligne d'entrée utilisateur."""
    if prompt:
        print(prompt)
    return input().strip()


def fortran_float(s):
    """Convertit un nombre Fortran (notation d/D) en float Python."""
    return float(s.replace('d', 'e').replace('D', 'e'))


def parse_floats(s):
    """Parse une chaîne en liste de floats (supporte la notation Fortran d0)."""
    return [fortran_float(x) for x in s.split()]


def parse_ints(s):
    """Parse une chaîne en liste d'ints."""
    return [int(fortran_float(x)) for x in s.split()]


class FortranFile:
    """Émule la lecture Fortran list-directed (read(unit,*)).

    - Saute les lignes vides
    - Peut lire N tokens sur plusieurs lignes
    - Gère les chaînes entre quotes (simple ou double)
    """

    def __init__(self, f):
        self.f = f
        self._buffer = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    @staticmethod
    def _tokenize(line):
        """Découpe une ligne en tokens, en gérant les chaînes quotées."""
        tokens = []
        i = 0
        while i < len(line):
            while i < len(line) and line[i] in ' \t,':
                i += 1
            if i >= len(line):
                break
            if line[i] in "'\"":
                quote = line[i]
                j = i + 1
                while j < len(line) and line[j] != quote:
                    j += 1
                tokens.append(line[i + 1:j])
                i = j + 1
            else:
                j = i
                while j < len(line) and line[j] not in ' \t,':
                    j += 1
                tokens.append(line[i:j])
                i = j
        return tokens

    def readline(self):
        """Lit la prochaine ligne non vide (comme Fortran)."""
        while True:
            line = self.f.readline()
            if not line:
                return ''
            if line.strip():
                return line

    def read_tokens(self, n):
        """Lit exactement n tokens, en continuant sur les lignes suivantes si nécessaire."""
        result = list(self._buffer)
        self._buffer = []
        while len(result) < n:
            line = self.readline().strip()
            if not line:
                break
            result.extend(self._tokenize(line))
        if len(result) > n:
            self._buffer = result[n:]
            result = result[:n]
        return result


# =============================================================================
# choles - Décomposition de Cholesky
# =============================================================================
def choles(m, a):
    """Décomposition de Cholesky. Retourne R tel que A = R^T * R."""
    r = np.zeros((m, m), dtype=np.float64)
    r[0, 0] = np.sqrt(a[0, 0])
    r[0, 1:m] = a[0, 1:m] / r[0, 0]
    for i in range(1, m - 1):
        r[i, i] = np.sqrt(a[i, i] - np.sum(r[:i, i]**2))
        vec_im = np.sum(r[:i, i].reshape(-1, 1) * r[:i, i+1:m], axis=0)
        r[i, i+1:m] = (a[i, i+1:m] - vec_im) / r[i, i]
    r[m-1, m-1] = np.sqrt(a[m-1, m-1] - np.sum(r[:m-1, m-1]**2))
    return r


# =============================================================================
# remp - Remplissage aléatoire de la matrice de phase space
# =============================================================================
def remp(m, n_part):
    """Génère une distribution aléatoire normalisée dans l'hyper-sphère."""
    np.random.seed(42)
    ps = np.zeros((m, n_part), dtype=np.float64)

    vtab = 10.0 * np.ones((m, n_part - 1), dtype=np.float64)
    while True:
        vec_sum = np.sum(vtab * vtab, axis=0)
        vmax = np.max(vec_sum)
        if vmax <= 1.0:
            break
        vmask = vec_sum > 1.0
        for j in range(m):
            rand_vals = np.random.random(n_part - 1)
            vtab[j, vmask] = 2.0 * rand_vals[vmask] - 1.0

    ps[:, 1:] = vtab
    ps[:, 0] = 0.0

    # Calcul de la matrice de corrélation
    c = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(m):
            c[i, j] = np.sum(ps[i, :] * ps[j, :]) / n_part

    r = choles(m, c)
    tr = r.T

    ps[0, :] = ps[0, :] / tr[0, 0]
    for i in range(1, m):
        ps[i, :] = (ps[i, :] - np.sum(
            tr[i, :i].reshape(-1, 1) * ps[:i, :], axis=0
        )) / tr[i, i]

    return ps


# =============================================================================
# addh - Ajout d'harmoniques à la carte de champ magnétique
# =============================================================================
def addh(rmin, rmax, nr, ntet, tbz, d2bz, dbzx, dbzy, ncoh, coh, phihar):
    """Ajoute des harmoniques à la carte de champ."""
    pi = np.pi
    dr = (rmax - rmin) / (nr - 1.0)
    dteta = 360.0 / (ntet - 1.0) * pi / 180.0
    phih1 = phihar[0] * pi / 180.0
    phih2 = phihar[1] * pi / 180.0

    if ncoh[0] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntet - 1):
                teta = j * dteta
                x = r * np.cos(teta)
                y = r * np.sin(teta)
                r2 = x**2 + y**2
                co1 = np.cos(-phih1)
                si1 = np.sin(-phih1)
                s = x * co1 - y * si1
                pb = coh[0, 0] * s
                px_val = coh[0, 0] * co1
                py_val = -coh[0, 0] * si1
                p2 = 0.0
                for k in range(1, ncoh[0]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coh[k, 0] * r2pk1 * s
                    px_val += coh[k, 0] * r2pk2 * (2.0 * x * k1 * s + r2 * co1)
                    py_val += coh[k, 0] * r2pk2 * (2.0 * y * k1 * s - r2 * si1)
                    p2 += coh[k, 0] * r2pk2 * 4.0 * s * k1 * (k + 1)
                tbz[i, j] += pb
                dbzx[i, j] += px_val
                dbzy[i, j] += py_val
                d2bz[i, j] += p2

    if ncoh[1] != 0:
        for i in range(nr):
            r = rmin + i * dr
            for j in range(ntet - 1):
                teta = j * dteta
                x = r * np.cos(teta)
                y = r * np.sin(teta)
                r2 = x**2 + y**2
                co2 = np.cos(-phih2)
                si2 = np.sin(-phih2)
                s = (x**2 - y**2) * co2 - 2.0 * x * y * si2
                s2 = 2.0 * x * co2 - 2.0 * y * si2
                s3 = -2.0 * y * co2 - 2.0 * x * si2
                pb = coh[0, 1] * s
                px_val = coh[0, 1] * s2
                py_val = -coh[0, 1] * s3
                p2 = 0.0
                for k in range(1, ncoh[1]):
                    k1 = k
                    k2 = k - 1
                    r2pk1 = r2**k1
                    if k == 1 and r2 <= 0.0:
                        r2pk2 = 1.0
                    else:
                        r2pk2 = r2**k2
                    pb += coh[k, 1] * r2pk1 * s
                    px_val += coh[k, 1] * r2pk2 * (2.0 * x * k1 * s + r2 * s2)
                    py_val += coh[k, 1] * r2pk2 * (2.0 * y * k1 * s + r2 * s3)
                    p2 += coh[k, 1] * r2pk2 * 4.0 * s * k1 * (k + 2)
                tbz[i, j] += pb
                dbzx[i, j] += px_val
                dbzy[i, j] += py_val
                d2bz[i, j] += p2

    # Périodicité
    tbz[:, ntet - 1] = tbz[:, 0]
    dbzx[:, ntet - 1] = dbzx[:, 0]
    dbzy[:, ntet - 1] = dbzy[:, 0]
    d2bz[:, ntet - 1] = d2bz[:, 0]


# =============================================================================
# boite - Vérification que les particules sont dans la boîte
# =============================================================================
def boite(d):
    """Vérifie quelles particules sont dans la boîte de simulation."""
    d.r = np.sqrt(d.x * d.x + d.y * d.y)
    d.box = ((d.r > d.rmin) & (d.r < (d.rmax - 0.02)) &
             (d.z > d.zmin) & (d.z < d.zmax))
    d.box = d.box & d.ok
    d.nbox = np.count_nonzero(d.box)


def box1(d):
    """Vérifie les particules dans la zone centrale."""
    d.centre = d.box & (d.r < d.rmax1) & (np.abs(d.y) < (d.ymax1 - d.hy1))
    d.nbox1 = np.count_nonzero(d.centre)


def box2(d):
    """Vérifie les particules dans la zone périphérique."""
    d.periph = ((d.vecx > d.xmin2) & (d.vecx < d.xmax2) &
                (d.vecy > d.ymin2) & (d.vecy < d.ymax2) &
                d.box & ~d.centre)
    d.nbox2 = np.count_nonzero(d.periph)


def boxdf(d, i):
    """Vérifie les particules dans un déflecteur électrique."""
    # i est 0-based en Python
    cond_in = (((d.x - d.xdfi[i]) * d.uydfi[i] - (d.y - d.ydfi[i]) * d.uxdfi[i]) <= 0.0)
    cond_out = (((d.x - d.xdfo[i]) * d.uydfo[i] - (d.y - d.ydfo[i]) * d.uxdfo[i]) >= 0.0)

    if i > 0:
        d.extra = d.sortant & cond_in & cond_out & d.ok
        d.nextra = np.count_nonzero(d.extra)
        return

    if d.iteext == 0:
        d.extra = (d.r > d.rdf[i]) & cond_in & cond_out
        d.extra = d.extra | (d.sortant & cond_in & cond_out)
        d.sortant = d.sortant | d.extra
    else:
        d.extra = (d.ite >= d.iteext) & cond_in & cond_out
        d.extra = d.extra | (d.sortant & cond_in & cond_out)
        d.sortant = d.sortant | d.extra

    d.sortant = d.sortant & d.ok
    d.extra = d.extra & d.ok
    d.nextra = np.count_nonzero(d.extra)


def boxdm(d, i):
    """Vérifie les particules dans un dipôle magnétique."""
    cond_in = (((d.x - d.xdm[i]) * d.uydmi[i] - (d.y - d.ydm[i]) * d.uxdmi[i]) <= 0.0)
    cond_out = (((d.x - d.xdm[i]) * d.uydmo[i] - (d.y - d.ydm[i]) * d.uxdmo[i]) >= 0.0)
    d.extra = d.sortant & cond_in & cond_out & d.ok
    d.nextra = np.count_nonzero(d.extra)


def boxgex(d, i):
    """Vérifie les particules aux grilles d'extraction."""
    cond_cross = (((d.x - d.xgex[i]) * d.uygex[i] -
                   (d.y - d.ygex[i]) * d.uxgex[i]) <= 0.0)

    if d.ndf == 0 and d.ndm == 0 and d.ngex == 1:
        if d.ite < d.iteext:
            d.franchi[:] = False
            return
        d.franchi = (d.indgex == 0) & d.ok & cond_cross
        reset_mask = (d.indgex == 1) & (
            ((d.x - d.xgex[i]) * d.uygex[i] -
             (d.y - d.ygex[i]) * d.uxgex[i]) > 0.0)
        d.indgex[reset_mask] = 0
        return

    d.franchi = d.sortant & d.ok & (d.indgex == i) & cond_cross


# =============================================================================
# piege - Piégeage des ions
# =============================================================================
def piege(d):
    """Gère le piégeage des ions dans les obstacles."""
    if d.nbrpie == 0:
        return
    mask_z = (d.z < d.zmin) | (d.z > d.zmax)
    d.ok[mask_z] = False

    if d.ite <= d.itepiem:
        for k in range(d.nbrpie):
            rpie2 = d.rpie[k]**2
            mask_pie = ((d.x - d.xpie[k])**2 + (d.y - d.ypie[k])**2) <= rpie2
            d.ok[mask_pie] = False

    d.nbrok = np.count_nonzero(d.ok)
    if d.iflpie == 1:
        not_ok = ~d.ok
        d.x[not_ok] = d.x[0]
        d.y[not_ok] = d.y[0]
        d.z[not_ok] = d.z[0]
        d.px[not_ok] = d.px[0]
        d.py[not_ok] = d.py[0]
        d.pz[not_ok] = d.pz[0]
        d.ga[not_ok] = d.ga[0]


# =============================================================================
# pousse - Boris pusher (avance des particules)
# =============================================================================
def pousse(d):
    """Avance les particules avec le Boris pusher."""
    if d.iflpie == 0:
        m = d.ok  # masque
        vec_m = d.qdts2m0[m] / d.ga[m]
        a1_m = vec_m * d.bx[m]
        a2_m = vec_m * d.by[m]
        a3_m = vec_m * d.bz[m]
        v1 = a1_m * a1_m
        v2 = a2_m * a2_m
        v3 = a3_m * a3_m
        vec_m = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1_m * a2_m
        s21 = s12 - a3_m
        s12 = s12 + a3_m
        s13 = a1_m * a3_m
        s31 = s13 + a2_m
        s13 = s13 - a2_m
        s23 = a2_m * a3_m
        s32 = s23 - a1_m
        s23 = s23 + a1_m

        fx_m = d.qdtsm0c[m] * d.ex[m]
        fy_m = d.qdtsm0c[m] * d.ey[m]
        fz_m = d.qdtsm0c[m] * d.ez[m]

        v4 = vec_m * (fx_m + 2.0 * d.px[m])
        v5 = vec_m * (fy_m + 2.0 * d.py[m])
        v6 = vec_m * (fz_m + 2.0 * d.pz[m])
        d.px[m] = d.px[m] - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx_m
        d.py[m] = d.py[m] + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy_m
        d.pz[m] = d.pz[m] + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz_m

        d.ga[m] = np.sqrt(1.0 + d.px[m]**2 + d.py[m]**2 + d.pz[m]**2)
        vec_m = d.cdt / d.ga[m]
        d.x[m] = d.x[m] + vec_m * d.px[m]
        d.y[m] = d.y[m] + vec_m * d.py[m]
        d.z[m] = d.z[m] + vec_m * d.pz[m]
    else:
        # Pas de masque, toutes les particules
        vec = d.qdts2m0 / d.ga
        a1 = vec * d.bx
        a2 = vec * d.by
        a3 = vec * d.bz
        v1 = a1 * a1
        v2 = a2 * a2
        v3 = a3 * a3
        vec = 1.0 / (1.0 + v1 + v2 + v3)
        s12 = a1 * a2
        s21 = s12 - a3
        s12 = s12 + a3
        s13 = a1 * a3
        s31 = s13 + a2
        s13 = s13 - a2
        s23 = a2 * a3
        s32 = s23 - a1
        s23 = s23 + a1

        fx = d.qdtsm0c * d.ex
        fy = d.qdtsm0c * d.ey
        fz = d.qdtsm0c * d.ez

        v4 = vec * (fx + 2.0 * d.px)
        v5 = vec * (fy + 2.0 * d.py)
        v6 = vec * (fz + 2.0 * d.pz)
        d.px = d.px - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx
        d.py = d.py + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy
        d.pz = d.pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz

        d.ga = np.sqrt(1.0 + d.px**2 + d.py**2 + d.pz**2)
        vec = d.cdt / d.ga
        d.x = d.x + vec * d.px
        d.y = d.y + vec * d.py
        d.z = d.z + vec * d.pz


# =============================================================================
# calcbe - Calcul complet du champ B et E
# =============================================================================
def calcbe(d):
    """Calcul des champs magnétique et électrique pour toutes les particules."""
    np_ = d.np_

    # Calcul de B
    nbiz = np.count_nonzero((d.x / d.r > 1.0) | (d.x / d.r < -1.0))
    if nbiz != 0:
        print(f'nbiz : {nbiz}')
        for i in range(np_):
            rap = d.x[i] / d.r[i]
            if rap > 1.0 or rap < -1.0:
                print(f'{i} {d.x[i]} {d.y[i]} {d.r[i]} {rap}')

    d.veca = d.x / d.r
    mask_normal = (d.veca > -1.0) & (d.veca < 1.0)
    mask_ge1 = d.veca >= 1.0
    mask_lem1 = d.veca <= -1.0
    d.teta[mask_normal] = np.arccos(d.veca[mask_normal])
    d.teta[mask_ge1] = np.arccos(1.0)
    d.teta[mask_lem1] = np.arccos(-1.0)
    mask_yneg = d.y < 0.0
    d.teta[mask_yneg] = 2.0 * d.pi - d.teta[mask_yneg]

    d.ivec[:] = 0
    d.jvec[:] = 0
    mask = d.box
    d.ivec[mask] = ((d.r[mask] - d.rmin) * d.hrm1).astype(np.int32)
    d.jvec[mask] = (d.teta[mask] * d.htetm1).astype(np.int32)
    d.veca[:] = d.r - d.rmin - d.ivec * d.hr
    d.vecb_arr[:] = d.teta - d.jvec * d.htet
    d.vecd1[:] = d.veca - d.hr
    d.vecd2[:] = d.vecb_arr - d.htet
    d.vecz[:] = d.z * d.z * 0.5

    for k in range(np_):
        ik = d.ivec[k]
        jk = d.jvec[k]
        d.bz[k] = (d.vecd1[k] * d.vecd2[k] * (d.tbz[ik, jk] -
                   d.vecz[k] * d.d2bz[ik, jk]) -
                   d.veca[k] * d.vecd2[k] * (d.tbz[ik+1, jk] -
                   d.vecz[k] * d.d2bz[ik+1, jk]) +
                   d.veca[k] * d.vecb_arr[k] * (d.tbz[ik+1, jk+1] -
                   d.vecz[k] * d.d2bz[ik+1, jk+1]) -
                   d.vecd1[k] * d.vecb_arr[k] * (d.tbz[ik, jk+1] -
                   d.vecz[k] * d.d2bz[ik, jk+1])) * d.hrhtetm1
        d.bx[k] = (d.vecd1[k] * d.vecd2[k] * d.dbzx[ik, jk] -
                   d.veca[k] * d.vecd2[k] * d.dbzx[ik+1, jk] +
                   d.veca[k] * d.vecb_arr[k] * d.dbzx[ik+1, jk+1] -
                   d.vecd1[k] * d.vecb_arr[k] * d.dbzx[ik, jk+1]
                   ) * d.z[k] * d.hrhtetm1
        d.by[k] = (d.vecd1[k] * d.vecd2[k] * d.dbzy[ik, jk] -
                   d.veca[k] * d.vecd2[k] * d.dbzy[ik+1, jk] +
                   d.veca[k] * d.vecb_arr[k] * d.dbzy[ik+1, jk+1] -
                   d.vecd1[k] * d.vecb_arr[k] * d.dbzy[ik, jk+1]
                   ) * d.z[k] * d.hrhtetm1

    mask_bz0 = (d.box & (d.bz == 0.0)) | (~d.box)
    d.bx[mask_bz0] = 0.0
    d.by[mask_bz0] = 0.0
    d.bz[mask_bz0] = -0.00001

    mask_rle = d.r <= d.rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.tbz[0, 0]

    # Calcul de E au centre
    d.ex[:] = 0.0
    d.ey[:] = 0.0
    d.ez[:] = 0.0
    box1(d)

    if d.nbox1 > 0:
        m_c = d.centre
        d.ivec[m_c] = ((d.x[m_c] - d.xmin1) * d.hx1m1).astype(np.int32)
        d.jvec[m_c] = ((d.y[m_c] - d.ymin1) * d.hy1m1).astype(np.int32)
        d.veca[m_c] = d.x[m_c] - d.xmin1 - d.ivec[m_c] * d.hx1
        d.vecb_arr[m_c] = d.y[m_c] - d.ymin1 - d.jvec[m_c] * d.hy1
        d.vecd1[m_c] = d.veca[m_c] - d.hx1
        d.vecd2[m_c] = d.vecb_arr[m_c] - d.hy1
        d.vecp1[m_c] = d.vecd1[m_c] * d.vecd2[m_c]
        d.vecp2[m_c] = -d.veca[m_c] * d.vecd2[m_c]
        d.vecp3[m_c] = d.veca[m_c] * d.vecb_arr[m_c]
        d.vecp4[m_c] = -d.vecd1[m_c] * d.vecb_arr[m_c]
        d.vecz[m_c] = d.z[m_c]**2 * 0.5

        for k in range(np_):
            if d.centre[k]:
                ik = d.ivec[k]
                jk = d.jvec[k]
                d.vecb1[k] = d.d2pot1[ik, jk]
                d.vecb2[k] = d.d2pot1[ik+1, jk]
                d.vecb3[k] = d.d2pot1[ik+1, jk+1]
                d.vecb4[k] = d.d2pot1[ik, jk+1]
                vz = d.vecz[k]
                d.vex1[k] = (-(d.pot1[ik+1, jk] - d.pot1[ik-1, jk]) +
                             vz * (d.d2pot1[ik+1, jk] - d.d2pot1[ik-1, jk])
                             ) * 0.5 / d.hx1
                d.vex2[k] = (-(d.pot1[ik+2, jk] - d.pot1[ik, jk]) +
                             vz * (d.d2pot1[ik+2, jk] - d.d2pot1[ik, jk])
                             ) * 0.5 / d.hx1
                d.vex3[k] = (-(d.pot1[ik+2, jk+1] - d.pot1[ik, jk+1]) +
                             vz * (d.d2pot1[ik+2, jk+1] - d.d2pot1[ik, jk+1])
                             ) * 0.5 / d.hx1
                d.vex4[k] = (-(d.pot1[ik+1, jk+1] - d.pot1[ik-1, jk+1]) +
                             vz * (d.d2pot1[ik+1, jk+1] - d.d2pot1[ik-1, jk+1])
                             ) * 0.5 / d.hx1
                d.vey1[k] = (-(d.pot1[ik, jk+1] - d.pot1[ik, jk-1]) +
                             vz * (d.d2pot1[ik, jk+1] - d.d2pot1[ik, jk-1])
                             ) * 0.5 / d.hy1
                d.vey2[k] = (-(d.pot1[ik+1, jk+1] - d.pot1[ik+1, jk-1]) +
                             vz * (d.d2pot1[ik+1, jk+1] - d.d2pot1[ik+1, jk-1])
                             ) * 0.5 / d.hy1
                d.vey3[k] = (-(d.pot1[ik+1, jk+2] - d.pot1[ik+1, jk]) +
                             vz * (d.d2pot1[ik+1, jk+2] - d.d2pot1[ik+1, jk])
                             ) * 0.5 / d.hy1
                d.vey4[k] = (-(d.pot1[ik, jk+2] - d.pot1[ik, jk]) +
                             vz * (d.d2pot1[ik, jk+2] - d.d2pot1[ik, jk])
                             ) * 0.5 / d.hy1

        nto = int(d.omhf * d.t / (2.0 * d.pi))
        omhft = d.omhf * d.t - nto * 2.0 * d.pi
        pcright = d.pothfd * np.cos(omhft + d.dphid)
        pcleft = d.pothfg * np.cos(omhft + d.dphig)
        d.vecpc[d.x >= 0.0] = pcright
        d.vecpc[d.x < 0.0] = pcleft

        m_c = d.centre
        d.ex[m_c] = d.vecpc[m_c] * (
            d.vex1[m_c] * d.vecp1[m_c] + d.vex2[m_c] * d.vecp2[m_c] +
            d.vex3[m_c] * d.vecp3[m_c] + d.vex4[m_c] * d.vecp4[m_c]
        ) * d.hx1hy1m1
        d.ey[m_c] = d.vecpc[m_c] * (
            d.vey1[m_c] * d.vecp1[m_c] + d.vey2[m_c] * d.vecp2[m_c] +
            d.vey3[m_c] * d.vecp3[m_c] + d.vey4[m_c] * d.vecp4[m_c]
        ) * d.hx1hy1m1
        d.ez[m_c] = d.z[m_c] * d.vecpc[m_c] * (
            d.vecp1[m_c] * d.vecb1[m_c] + d.vecp2[m_c] * d.vecb2[m_c] +
            d.vecp3[m_c] * d.vecb3[m_c] + d.vecp4[m_c] * d.vecb4[m_c]
        ) * d.hx1hy1m1

    # Calcul de E périphérique
    if (d.nbox == d.nbox1) and (d.ndm == 0) and (d.ndf == 0) and (d.ngex == 1):
        pass  # goto 300
    elif d.nbox == d.nbox1:
        _calcbe_vertical(d)
        _calcbe_extraction(d)
        return
    else:
        dc1 = np.cos(d.teta1)
        ds1 = np.sin(d.teta1)
        mask_nc = ~d.centre
        d.vecx[mask_nc] = np.abs(d.x[mask_nc]) * dc1 + (np.abs(d.y[mask_nc]) - d.yomega) * ds1
        d.vecy[mask_nc] = -np.abs(d.x[mask_nc]) * ds1 + (np.abs(d.y[mask_nc]) - d.yomega) * dc1

        box2(d)
        nto = int(d.omhf * d.t / (2.0 * d.pi))
        omhft = d.omhf * d.t - nto * 2.0 * d.pi

        if d.nbox2 > 0:
            if d.nbox1 == 0:
                pcright = d.pothfd * np.cos(omhft + d.dphid)
                pcleft = d.pothfg * np.cos(omhft + d.dphig)
                d.vecpc[d.x >= 0.0] = pcright
                d.vecpc[d.x < 0.0] = pcleft

            m_p = d.periph
            d.ivec[m_p] = ((d.vecx[m_p] - d.xmin2) * d.hx2m1).astype(np.int32)
            d.jvec[m_p] = ((d.vecy[m_p] - d.ymin2) * d.hy2m1).astype(np.int32)
            d.veca[m_p] = d.vecx[m_p] - d.xmin2 - d.ivec[m_p] * d.hx2
            d.vecb_arr[m_p] = d.vecy[m_p] - d.ymin2 - d.jvec[m_p] * d.hy2
            d.vecd1[m_p] = d.veca[m_p] - d.hx2
            d.vecd2[m_p] = d.vecb_arr[m_p] - d.hy2
            d.vecp1[m_p] = d.vecd1[m_p] * d.vecd2[m_p]
            d.vecp2[m_p] = -d.veca[m_p] * d.vecd2[m_p]
            d.vecp3[m_p] = d.veca[m_p] * d.vecb_arr[m_p]
            d.vecp4[m_p] = -d.vecd1[m_p] * d.vecb_arr[m_p]
            d.vecz[m_p] = d.z[m_p]**2 * 0.5

            for k in range(np_):
                if d.periph[k]:
                    ik = d.ivec[k]
                    jk = d.jvec[k]
                    d.vec1[k] = d.pot2[ik, jk]
                    d.vec2[k] = d.pot2[ik+1, jk]
                    d.vec3[k] = d.pot2[ik+1, jk+1]
                    d.vec4[k] = d.pot2[ik, jk+1]
                    d.vecb1[k] = d.d2pot2[ik, jk]
                    d.vecb2[k] = d.d2pot2[ik+1, jk]
                    d.vecb3[k] = d.d2pot2[ik+1, jk+1]
                    d.vecb4[k] = d.d2pot2[ik, jk+1]

            d.ex[m_p] = -d.vecpc[m_p] * (
                d.vecd2[m_p] * (d.vec1[m_p] - d.vec2[m_p] -
                                d.vecz[m_p] * (d.vecb1[m_p] - d.vecb2[m_p])) +
                d.vecb_arr[m_p] * (d.vec3[m_p] - d.vec4[m_p] +
                                   d.vecz[m_p] * (d.vecb3[m_p] - d.vecb4[m_p]))
            ) * d.hx2hy2m1
            d.ey[m_p] = -d.vecpc[m_p] * (
                d.vecd1[m_p] * (d.vec1[m_p] - d.vec4[m_p] -
                                d.vecz[m_p] * (d.vecb1[m_p] - d.vecb4[m_p])) -
                d.veca[m_p] * (d.vec2[m_p] - d.vec3[m_p] +
                               d.vecz[m_p] * (d.vecb2[m_p] - d.vecb3[m_p]))
            ) * d.hx2hy2m1

            # Rotation inverse
            vec_tmp = d.ex[m_p] * dc1 - d.ey[m_p] * ds1
            d.ey[m_p] = d.ex[m_p] * ds1 + d.ey[m_p] * dc1
            d.ex[m_p] = vec_tmp

            d.ez[m_p] = d.z[m_p] * d.vecpc[m_p] * (
                d.vecp1[m_p] * d.vecb1[m_p] + d.vecp2[m_p] * d.vecb2[m_p] +
                d.vecp3[m_p] * d.vecb3[m_p] + d.vecp4[m_p] * d.vecb4[m_p]
            ) * d.hx2hy2m1

            # Symétries
            m1 = m_p & (d.x < 0.0) & (d.y > 0.0)
            d.ex[m1] = d.heps * d.ex[m1]
            d.ey[m1] = -d.heps * d.ey[m1]

            m2 = m_p & (d.x > 0.0) & (d.y < 0.0)
            d.ey[m2] = -d.ey[m2]

            m3 = m_p & (d.x < 0.0) & (d.y < 0.0)
            d.ex[m3] = d.heps * d.ex[m3]
            d.ey[m3] = d.heps * d.ey[m3]

            m4 = m_p & (d.x < 0.0)
            d.ez[m4] = -d.heps * d.ez[m4]

    # Tri vertical
    _calcbe_vertical(d)

    # Extraction
    _calcbe_extraction(d)


def _calcbe_vertical(d):
    """Ajout du test tri vertical."""
    d.pot_ver_vrai = 0.0
    if d.h_ver == 0.0:
        return
    mask_ver = ((d.r >= d.rmin_ver) & (d.r <= d.rmax_ver) &
                (d.teta >= d.teta1_ver) & (d.teta <= d.teta2_ver))
    d.ez[mask_ver] = (d.ez[mask_ver] -
        d.pot_ver * (1.0 - d.dpotspot) / d.h_ver / 2.0 *
        np.cos(2.0 * d.omhf * d.t * (1.0 + d.anuz / 2.0 / d.harm) +
               2.0 * d.phi_ver - d.pi / 2.0) -
        d.pot_ver * (1.0 + d.dpotspot) / d.h_ver / 2.0 *
        np.cos(2.0 * d.omhf * d.t * (1.0 - d.anuz / 2.0 / d.harm) +
               2.0 * d.phi_ver - d.pi / 2.0))
    if (d.r[0] >= d.rmin_ver and d.r[0] <= d.rmax_ver and
            d.teta[0] >= d.teta1_ver and d.teta[0] <= d.teta2_ver):
        d.pot_ver_vrai = -d.ez[0] * d.h_ver


def _calcbe_extraction(d):
    """Calcul de l'extraction."""
    for i in range(d.ndf):
        boxdf(d, i)
        if d.nextra != 0:
            mask = d.extra
            vec_p = np.sqrt(d.px[mask]**2 + d.py[mask]**2)
            d.ex[mask] = d.edf[i] * d.py[mask] / vec_p
            d.ey[mask] = -d.edf[i] * d.px[mask] / vec_p

    for i in range(d.ndm):
        boxdm(d, i)
        if d.nextra != 0:
            mask = d.extra
            vec_dm = np.sqrt((d.x[mask] - d.xdm[i])**2 + (d.y[mask] - d.ydm[i])**2)
            vec1_dm = vec_dm - d.rdm[i]
            d.bx[mask] += d.z[mask] * d.bdm[i] * (d.x[mask] - d.xdm[i]) / vec_dm
            d.by[mask] += d.z[mask] * d.bdm[i] * (d.y[mask] - d.ydm[i]) / vec_dm
            d.bz[mask] += vec1_dm * d.bdm[i]

    # Grilles d'extraction
    if d.rex < 0.0 and d.ngex > 0:
        if d.ite >= d.itegex:
            for kg in range(d.ngex):
                boxgex(d, kg)
                d.indgex[d.franchi] = kg + 1  # 1-based in data
                for i in range(d.n):
                    if d.franchi[i]:
                        if (70 + kg + 1) in files_dict:
                            files_dict[70 + kg + 1].write(
                                f' {i+1:5d} {d.ite:6d}'
                                f' {d.x[i]:15.8e} {d.y[i]:15.8e} {d.z[i]:15.8e}'
                                f' {d.px[i]:15.8e} {d.py[i]:15.8e} {d.pz[i]:15.8e}'
                                f' {d.bz[i]:8.5f}\n')


# =============================================================================
# calcb - Calcul de B simplifié (utilisé pendant init)
# =============================================================================
def calcb(d):
    """Calcul simplifié du champ magnétique (pendant l'initialisation)."""
    d.teta = np.arccos(np.clip(d.x / d.r, -1.0, 1.0))
    d.teta[d.y < 0.0] = 2.0 * d.pi - d.teta[d.y < 0.0]

    mask = d.box
    d.ivec[mask] = ((d.r[mask] - d.rmin) / d.hr).astype(np.int32)
    d.jvec[mask] = (d.teta[mask] / d.htet).astype(np.int32)
    d.veca[mask] = d.r[mask] - d.rmin - d.ivec[mask] * d.hr
    d.vecb_arr[mask] = d.teta[mask] - d.jvec[mask] * d.htet
    d.vecd1[mask] = d.veca[mask] - d.hr
    d.vecd2[mask] = d.vecb_arr[mask] - d.htet
    d.vecp1[mask] = d.vecd1[mask] * d.vecd2[mask]
    d.vecp2[mask] = -d.veca[mask] * d.vecd2[mask]
    d.vecp3[mask] = d.veca[mask] * d.vecb_arr[mask]
    d.vecp4[mask] = -d.vecd1[mask] * d.vecb_arr[mask]
    d.vecz[mask] = d.z[mask]**2 / 2.0

    for k in range(d.np_):
        if d.box[k]:
            ik = d.ivec[k]
            jk = d.jvec[k]
            vk = d.vecz[k]
            d.vecb1[k] = d.tbz[ik, jk]
            d.vecb2[k] = d.tbz[ik+1, jk]
            d.vecb3[k] = d.tbz[ik+1, jk+1]
            d.vecb4[k] = d.tbz[ik, jk+1]
            d.vec1[k] = -vk * d.d2bz[ik, jk]
            d.vec2[k] = -vk * d.d2bz[ik+1, jk]
            d.vec3[k] = -vk * d.d2bz[ik+1, jk+1]
            d.vec4[k] = -vk * d.d2bz[ik, jk+1]
            d.vecx1[k] = d.dbzx[ik, jk]
            d.vecx2[k] = d.dbzx[ik+1, jk]
            d.vecx3[k] = d.dbzx[ik+1, jk+1]
            d.vecx4[k] = d.dbzx[ik, jk+1]
            d.vecy1[k] = d.dbzy[ik, jk]
            d.vecy2[k] = d.dbzy[ik+1, jk]
            d.vecy3[k] = d.dbzy[ik+1, jk+1]
            d.vecy4[k] = d.dbzy[ik, jk+1]

    d.vec1[mask] += d.vecb1[mask]
    d.vec2[mask] += d.vecb2[mask]
    d.vec3[mask] += d.vecb3[mask]
    d.vec4[mask] += d.vecb4[mask]
    d.bz[mask] = (d.vecp1[mask] * d.vec1[mask] + d.vecp2[mask] * d.vec2[mask] +
                  d.vecp3[mask] * d.vec3[mask] + d.vecp4[mask] * d.vec4[mask]) / d.hrhtet
    d.bx[mask] = d.z[mask] * (
        d.vecp1[mask] * d.vecx1[mask] + d.vecp2[mask] * d.vecx2[mask] +
        d.vecp3[mask] * d.vecx3[mask] + d.vecp4[mask] * d.vecx4[mask]) / d.hrhtet
    d.by[mask] = d.z[mask] * (
        d.vecp1[mask] * d.vecy1[mask] + d.vecp2[mask] * d.vecy2[mask] +
        d.vecp3[mask] * d.vecy3[mask] + d.vecp4[mask] * d.vecy4[mask]) / d.hrhtet

    mask_rle = d.r <= d.rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.tbz[0, 0]


# =============================================================================
# caldxs - Calcul de l'espace de phase
# =============================================================================
def caldxs(d):
    """Calcul des coordonnées dans l'espace de phase."""
    for k in range(d.nbrfam):
        nbr = d.nbrion[k]
        i1 = d.index[k]
        i2 = i1 + nbr
        d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
        rc = -d.am0c[k] * d.p[i1] / d.q[i1] / d.bz[i1]
        d.dxd[0, :nbr] = d.x[i1:i2] - d.x[i1]
        d.dxd[4, :nbr] = d.y[i1:i2] - d.y[i1]
        ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
        tetai1 = np.arccos(d.x[i1] / ri1)
        if d.y[i1] < 0.0:
            tetai1 = 2.0 * d.pi - tetai1
        pri1 = (d.x[i1] * d.px[i1] + d.y[i1] * d.py[i1]) / ri1
        tpri1 = np.arcsin(pri1 / d.p[i1])
        co = np.cos(tetai1 - tpri1)
        si = np.sin(tetai1 - tpri1)
        d.vec[:nbr] = d.dxd[0, :nbr] * co + d.dxd[4, :nbr] * si
        d.dxd[4, :nbr] = -d.dxd[0, :nbr] * si + d.dxd[4, :nbr] * co
        d.dxd[0, :nbr] = d.vec[:nbr]
        for ii in range(nbr):
            d.dxs[4, ii] = rc * np.arctan(d.dxd[4, ii] / (rc + d.dxd[0, ii]))
        d.vec[:nbr] = np.sin(d.dxs[4, :nbr] / 2.0 / rc)**2
        d.dxs[0, :nbr] = (d.dxd[0, :nbr] + 2.0 * rc * d.vec[:nbr]) / (1.0 - 2.0 * d.vec[:nbr])
        d.vec1[:nbr] = np.cos(d.dxs[4, :nbr] / rc)
        d.vec2[:nbr] = np.sin(d.dxs[4, :nbr] / rc)
        d.dxd[1, :nbr] = d.px[i1:i2] - d.px[i1] * d.vec1[:nbr] + d.py[i1] * d.vec2[:nbr]
        d.dxd[5, :nbr] = d.py[i1:i2] - d.px[i1] * d.vec2[:nbr] - d.py[i1] * d.vec1[:nbr]
        d.vec[:nbr] = d.dxd[1, :nbr] * co + d.dxd[5, :nbr] * si
        d.dxd[5, :nbr] = -d.dxd[1, :nbr] * si + d.dxd[5, :nbr] * co
        d.dxd[1, :nbr] = d.vec[:nbr]
        d.dxs[1, :nbr] = d.dxd[1, :nbr] * d.vec1[:nbr] + d.dxd[5, :nbr] * d.vec2[:nbr]
        d.dxs[5, :nbr] = -d.dxd[1, :nbr] * d.vec2[:nbr] + d.dxd[5, :nbr] * d.vec1[:nbr]
        d.dxs[2, :nbr] = d.z[i1:i2] - d.z[i1]
        d.dxs[3, :nbr] = d.pz[i1:i2] - d.pz[i1]
        d.dxs[5, :nbr] = np.sqrt(
            (d.p[i1] + d.dxs[5, :nbr])**2 +
            d.dxs[1, :nbr]**2 + d.dxs[3, :nbr]**2) - d.p[i1]

        for i in range(6):
            for j in range(6):
                d.cor[i, j] = np.sum(d.dxs[i, :nbr] * d.dxs[j, :nbr]) / nbr

        tetai1 = tetai1 * 180.0 / d.pi
        ri1mm = ri1 * 1000.0
        d.t1[:] = 0.0
        d.t1[0, 0] = -1.0
        d.t1[1, 1] = -1.0
        for i in range(2, 6):
            d.t1[i, i] = 1.0
        d.t1[0, 5] = ri1 / d.p[i1]
        d.t1[4, 1] = -d.t1[0, 5]
        d.cordc = d.t1 @ d.cor @ d.t1.T


# =============================================================================
# calgri / calgric - Calcul sur grille
# =============================================================================
def calgri(d):
    """Calcul des coordonnées sur grille."""
    for k in range(d.nbrfam):
        nbr = d.nbrion[k]
        i1 = d.index[k]
        i2 = i1 + nbr
        d.p[i1:i2] = np.sqrt(d.px[i1:i2]**2 + d.py[i1:i2]**2 + d.pz[i1:i2]**2)
        d.vec[:nbr] = -d.am0c[k] * d.p[i1:i2] / d.q[i1:i2] / d.bz[i1]
        d.vec1[i1:i2] = d.x[i1:i2] + d.py[i1:i2] * d.dt / d.qdtsm0c[i1] / d.bz[i1]
        d.vec2[i1:i2] = d.y[i1:i2] - d.px[i1:i2] * d.dt / d.qdtsm0c[i1] / d.bz[i1]
        ux = (d.x[i1] - d.vec1[i1]) / d.vec[0]
        uy = (d.y[i1] - d.vec2[i1]) / d.vec[0]
        d.vec3[i1:i2] = ((d.x[i1] - d.vec1[i1:i2]) * ux +
                         (d.y[i1] - d.vec2[i1:i2]) * uy)
        d.vec4[i1:i2] = ((d.x[i1] - d.vec1[i1:i2])**2 +
                         (d.y[i1] - d.vec2[i1:i2])**2 - d.vec[:nbr]**2)
        d.rgri[i1:i2] = -d.vec3[i1:i2] + np.sqrt(np.abs(d.vec3[i1:i2]**2 - d.vec4[i1:i2]))
        d.xgri[i1:i2] = d.x[i1] + d.rgri[i1:i2] * ux
        d.ygri[i1:i2] = d.y[i1] + d.rgri[i1:i2] * uy
        d.vec3[i1:i2] = np.sqrt((d.xgri[i1:i2] - d.vec1[i1:i2])**2 +
                                (d.ygri[i1:i2] - d.vec2[i1:i2])**2)
        d.vec4[i1:i2] = np.sqrt((d.xgri[i1:i2] - d.vec1[i1])**2 +
                                (d.ygri[i1:i2] - d.vec2[i1])**2)
        d.vec5[i1:i2] = (
            (d.xgri[i1:i2] - d.vec1[i1:i2]) * (d.ygri[i1:i2] - d.vec2[i1]) -
            (d.ygri[i1:i2] - d.vec2[i1:i2]) * (d.xgri[i1:i2] - d.vec1[i1])
        ) / d.vec3[i1:i2] / d.vec4[i1:i2]
        d.tprgri[i1:i2] = np.arcsin(d.vec5[i1:i2])
        d.vec4[i1:i2] = np.sqrt((d.x[i1:i2] - d.vec1[i1:i2])**2 +
                                (d.y[i1:i2] - d.vec2[i1:i2])**2)
        d.vec6[i1:i2] = (
            (d.xgri[i1:i2] - d.vec1[i1:i2]) * (d.y[i1:i2] - d.vec2[i1:i2]) -
            (d.ygri[i1:i2] - d.vec2[i1:i2]) * (d.x[i1:i2] - d.vec1[i1:i2])
        ) / d.vec3[i1:i2] / d.vec4[i1:i2]
        d.phigri[i1:i2] = d.harm * np.arcsin(d.vec6[i1:i2])


def calgric(d):
    """Calcul sur grille pour corrélations."""
    nbr = d.nbrion[0]
    calgri(d)
    d.dxsgri[0, :nbr] = d.rgri[:nbr]
    d.dxsgri[1, :nbr] = d.tprgri[:nbr] * d.p[0]
    d.dxsgri[2, :nbr] = d.dxs[2, :nbr]
    d.dxsgri[3, :nbr] = d.dxs[3, :nbr]
    d.dxsgri[4, :nbr] = d.phigri[:nbr] * d.vec[:nbr] / d.harm
    d.dxsgri[5, :nbr] = d.dxs[5, :nbr]


# =============================================================================
# Variable globale pour les fichiers ouverts
# =============================================================================
files_dict = {}


# =============================================================================
# init - Initialisation du programme
# =============================================================================
def init(d, ficin):
    """Initialisation : lecture du fichier de données ficin, cartes de champ, etc."""
    files_dict.clear()

    d.pi = np.pi
    d.clum = 299792458.0
    d.rac8 = np.sqrt(8.0)

    d.ficin = ficin
    with open(d.ficin, 'r') as f20_raw:
      f20 = FortranFile(f20_raw)
      with f20:  # dummy context for indentation
        d.corps = f20.readline().strip().strip('"').strip("'")
        vals = f20.readline().split()
        d.npaspt = int(vals[0])
        vals = f20.readline().split()
        d.r1th = fortran_float(vals[0])
        d.ntours = int(vals[1])
        vals = f20.readline().split()
        d.fhf = fortran_float(vals[0])
        d.pothf = fortran_float(vals[1])
        d.teta1 = fortran_float(vals[2])
        d.harm = fortran_float(vals[3])
        d.dbsb = fortran_float(vals[4])
        vals = f20.readline().split()
        d.dphid = fortran_float(vals[0]) * d.pi / 180.0
        d.dphig = fortran_float(vals[1]) * d.pi / 180.0
        vals = f20.readline().split()
        d.pothfd = -fortran_float(vals[0])
        d.pothfg = -fortran_float(vals[1])
        vals = f20.readline().split()
        d.ncoh[0] = int(vals[0])
        if d.ncoh[0] != 0:
            vals = f20.readline().split()
            for i in range(d.ncoh[0]):
                d.coh[i, 0] = fortran_float(vals[i])
            d.phihar[0] = fortran_float(f20.readline().split()[0])
        vals = f20.readline().split()
        d.ncoh[1] = int(vals[0])
        if d.ncoh[1] != 0:
            vals = f20.readline().split()
            for i in range(d.ncoh[1]):
                d.coh[i, 1] = fortran_float(vals[i])
            d.phihar[1] = fortran_float(f20.readline().split()[0])

        d.yomega = 0.0
        if d.teta1 < 0.0:
            print('attention : gap en biais...')
            print('on remet teta1 positif')
            print('et on va lire le decalage...')
            d.teta1 = -d.teta1
            d.yomega = fortran_float(f20.readline().split()[0])
            print(f'decalage yomega : {d.yomega}')

        vals = f20.readline().split()
        d.zmin = fortran_float(vals[0])
        d.zmax = fortran_float(vals[1])
        d.rmax1 = fortran_float(vals[2])
        vals = f20.readline().split()
        d.rex = fortran_float(vals[0])
        d.tet1ex = fortran_float(vals[1])
        d.tet2ex = fortran_float(vals[2])
        d.eex = fortran_float(vals[3])
        d.ndf = 0
        d.ndm = 0

        if d.rex < 0.0:
            if d.tet1ex <= 0.0:
                d.iteext = 0
                print('extraction : le paquet sera coupe')
            else:
                d.iteext = int(d.tet1ex + 0.5)
                print('extraction : le paquet sera entier')
            d.ndf = int(f20.readline().split()[0])
            if d.ndf > 0:
                for i in range(d.ndf):
                    vals = f20.readline().split()
                    d.rdf[i] = fortran_float(vals[0])
                    d.edf[i] = fortran_float(vals[1])
                    d.xdfi[i] = fortran_float(vals[2])
                    d.ydfi[i] = fortran_float(vals[3])
                    d.xdfo[i] = fortran_float(vals[4])
                    d.ydfo[i] = fortran_float(vals[5])
                    tetdfid = fortran_float(vals[6])
                    tetdfod = fortran_float(vals[7])
                    d.tetdfi[i] = tetdfid * d.pi / 180.0
                    d.tetdfo[i] = tetdfod * d.pi / 180.0
                    d.uxdfi[i] = np.cos(d.tetdfi[i])
                    d.uydfi[i] = np.sin(d.tetdfi[i])
                    d.uxdfo[i] = np.cos(d.tetdfo[i])
                    d.uydfo[i] = np.sin(d.tetdfo[i])
            d.ndm = int(f20.readline().split()[0])
            if d.ndm > 0:
                for i in range(d.ndm):
                    vals = f20.readline().split()
                    d.rdm[i] = fortran_float(vals[0])
                    d.bdm[i] = fortran_float(vals[1])
                    d.xdm[i] = fortran_float(vals[2])
                    d.ydm[i] = fortran_float(vals[3])
                    tetdmid = fortran_float(vals[4])
                    tetdmod = fortran_float(vals[5])
                    d.tetdmi[i] = tetdmid * d.pi / 180.0
                    d.tetdmo[i] = tetdmod * d.pi / 180.0
                    d.uxdmi[i] = np.cos(d.tetdmi[i])
                    d.uydmi[i] = np.sin(d.tetdmi[i])
                    d.uxdmo[i] = np.cos(d.tetdmo[i])
                    d.uydmo[i] = np.sin(d.tetdmo[i])
            if d.tet2ex > 0.0:
                print('grilles extraction...')
                d.itegex = int(d.eex + 0.5)
                d.ngex = int(f20.readline().split()[0])
                for kg in range(d.ngex):
                    vals = f20.readline().split()
                    d.tetagex[kg] = fortran_float(vals[0])
                    d.xgex[kg] = fortran_float(vals[1]) / 1000.0
                    d.ygex[kg] = fortran_float(vals[2]) / 1000.0
                    d.xgexc[kg] = fortran_float(vals[3]) / 1000.0
                    d.ygexc[kg] = fortran_float(vals[4]) / 1000.0
                    d.drgex[kg] = np.sqrt((d.xgex[kg] - d.xgexc[kg])**2 +
                                          (d.ygex[kg] - d.ygexc[kg])**2)
                    d.uxgex[kg] = (d.xgex[kg] - d.xgexc[kg]) / d.drgex[kg]
                    d.uygex[kg] = (d.ygex[kg] - d.ygexc[kg]) / d.drgex[kg]
                d.ficdebgex = f20.readline().strip()
                for kg in range(d.ngex):
                    ficfingex = str(kg + 1)
                    d.ficgex[kg] = d.ficdebgex + '.' + ficfingex
                    fgex = open(d.ficgex[kg], 'w')
                    files_dict[70 + kg + 1] = fgex
                    fgex.write('teta x y xc yc : \n')
                    fgex.write(f'{d.tetagex[kg]} {d.xgex[kg]} {d.ygex[kg]}'
                               f' {d.xgexc[kg]} {d.ygexc[kg]}\n')

        d.dt = d.harm / d.fhf / d.npaspt
        d.omhf = 2.0 * d.pi * d.fhf
        d.omr = d.omhf / d.harm
        d.cdt = d.clum * d.dt
        d.teta1 = d.teta1 * d.pi / 180.0
        d.tet1ex = d.tet1ex * d.pi / 180.0
        d.tet2ex = d.tet2ex * d.pi / 180.0
        iharm = int(d.harm + 0.5)
        d.heps = 1.0
        if (iharm // 2) * 2 == iharm:
            d.heps = -1.0

        d.nbrfam = int(f20.readline().split()[0])
        if d.nbrfam > 1:
            d.n = int(f20.readline().split()[0])
        d.index[0] = 0  # 0-based
        d.np_ = 0

        for k in range(d.nbrfam):
            vals = f20.readline().split()
            d.nbrion[k] = int(vals[0])
            d.amat[k] = fortran_float(vals[1])
            d.a[k] = fortran_float(vals[2])
            d.qn[k] = fortran_float(vals[3])
            if k == 0:
                if d.nbrfam == 1:
                    d.n = d.nbrion[0]
                d.allo()
                d.sortant[:] = False
                d.indgex[:] = 0
                d.ok[:] = True
                d.nbrok = d.n
                d.nbrpie = 0
                d.iflpie = 1

            vals = f20.readline().split()
            d.r0[k] = fortran_float(vals[0])
            d.teta0[k] = fortran_float(vals[1])
            d.z0[k] = fortran_float(vals[2])
            d.pr0[k] = fortran_float(vals[3])
            d.pz0[k] = fortran_float(vals[4])
            w0mev = fortran_float(vals[5])
            d.phi0[k] = fortran_float(vals[6])
            if k == 0:
                d.w0mev1 = w0mev

            d.teta0[k] *= d.pi / 180.0
            d.pr0[k] *= 0.001
            d.pz0[k] *= 0.001
            d.np_ += d.nbrion[k]
            if d.np_ > d.n:
                print(f'n np : {d.n} {d.np_}')
                print('changer n ou np...')
                sys.exit(1)
            d.index[k + 1] = d.index[k] + d.nbrion[k]
            i1 = d.index[k]
            i2 = d.index[k] + d.nbrion[k]
            d.x[i1:i2] = d.r0[k] * np.cos(d.teta0[k])
            d.y[i1:i2] = d.r0[k] * np.sin(d.teta0[k])
            d.z[i1:i2] = d.z0[k]
            d.w0[k] = w0mev * 1.60217733e-13 * d.a[k]
            am0k = 1.6605402e-27 * (d.amat[k] - 0.0005486 * d.qn[k])
            d.am0[i1:i2] = am0k
            d.am0c[i1:i2] = d.am0[i1:i2] * d.clum
            d.am0c2[i1:i2] = d.am0[i1:i2] * d.clum**2
            d.am02c2[i1:i2] = d.am0c[i1:i2]**2
            d.ga[i1:i2] = 1.0 + d.w0[k] / d.am0c2[i1:i2]
            alpha = (d.pi / 2.0 + d.teta0[k]) - d.pr0[k]
            d.p[i1:i2] = np.sqrt((d.ga[i1:i2] - 1.0) * (d.ga[i1:i2] + 1.0))
            d.px[i1:i2] = d.p[i1:i2] * np.cos(alpha) * np.cos(d.pz0[k])
            d.py[i1:i2] = d.p[i1:i2] * np.sin(alpha) * np.cos(d.pz0[k])
            d.pz[i1:i2] = d.p[i1:i2] * np.sin(d.pz0[k])
            d.q[i1:i2] = d.qn[k] * 1.60217733e-19
            d.qdts2m0[i1:i2] = d.q[i1:i2] * d.dt / 2.0 / d.am0[i1:i2]
            d.qdtsm0c[i1:i2] = d.q[i1:i2] * d.dt / d.am0c[i1:i2]

        if d.tet2ex > 0.0 and d.ngex > 0:
            for kg in range(d.ngex):
                f = files_dict[70 + kg + 1]
                f.write('caracteristiques particules :\n')
                f.write(f'{d.np_} {d.a[0]} {d.amat[0]} {d.qn[0]}\n')
                f.write('caracteristiques HF :\n')
                f.write(f'{d.npaspt} {d.fhf} {d.harm}\n')

        d.iflcor = int(f20.readline().split()[0])
        d.ficsig = "pas de fichier correl..."
        sigma_matrices = []
        if d.iflcor != 0:
            m = 6
            d.ficsig = f20.readline().strip()
            d.espace = f20.readline().strip()
            f25 = open(d.ficsig, 'r')
            vals = f20.read_tokens(3)
            d.ficcor = vals[0]
            d.ficmat = vals[1]
            d.ficdc = vals[2]
            if d.ficcor != 'no':
                files_dict[17] = open(d.ficcor, 'w')
                files_dict[17].write(
                    '  tr     r      teta    dr    ds      dphi'
                    '   epsr    dwsw        wmpa     dw          dwdphi'
                    '       dz     epsz\n')
            if d.ficmat != 'no':
                files_dict[18] = open(d.ficmat, 'w')
            if d.ficdc != 'no':
                files_dict[19] = open(d.ficdc, 'w')
                files_dict[19].write(
                    '  tr     r      teta   drdc  dsdc  dphidc'
                    ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc\n')

            # Pre-lire les matrices sigma (calcul differe apres chargement champ magnetique)
            for k in range(d.nbrfam):
                cor_mat = np.zeros((6, 6), dtype=np.float64)
                for i in range(6):
                    vals_sig = f25.readline().split()
                    for j in range(6):
                        cor_mat[i, j] = fortran_float(vals_sig[j])
                sigma_matrices.append(cor_mat)
            f25.close()

        if d.np_ != d.n:
            print(f'n np : {d.n} {d.np_}')
            print('changer n ou np...')
            sys.exit(1)

        d.t0 = d.phi0[0] * d.pi / 180.0 / d.omhf
        d.t = d.t0

        vals = f20.readline().split()
        d.nbdt = int(vals[0])
        d.ifreq = int(vals[1])
        vals = f20.read_tokens(2)
        d.ficz = vals[0]
        d.ficxy = vals[1]
        vals = f20.read_tokens(3)
        d.ficb = vals[0]
        d.ficv1 = vals[1]
        d.ficv2 = vals[2]
        # Fin optionnelle : [b_factor] [ficloa ficsto]. Le Fortran ne lit que
        # ficloa ficsto ; b_factor n'est reconnu que s'il est numerique.
        tail = []
        while True:
            tok = f20.read_tokens(1)
            if not tok:
                break
            tail += tok
        d.b_factor = None
        if tail:
            try:
                d.b_factor = fortran_float(tail[0])
                tail = tail[1:]
            except ValueError:
                pass
        d.ficloa, d.ficsto = (tail + ['no', 'no'])[:2]

    # Lecture carte magnétique
    with open(d.ficb, 'r') as f21_raw:
      f21 = FortranFile(f21_raw)
      with f21:
        title = f21.readline().strip()
        vals = f21.readline().split()
        d.rmin = fortran_float(vals[0])
        d.rmax = fortran_float(vals[1])
        d.nr = int(vals[2])
        d.ntet = int(vals[3])

        if d.ntet != 91:
            print('lecture carte tosca.big ou carte mesuree ...')
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.tbz[i, j] = fortran_float(vals[j])
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.d2bz[i, j] = fortran_float(vals[j])
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.dbzx[i, j] = fortran_float(vals[j])
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.dbzy[i, j] = fortran_float(vals[j])
            d.htet = 2.0 * d.pi / (d.ntet - 1.0)
            ntetall = d.ntet
        else:
            print('carte tosca.small (1/8eme) ...')
            # Lecture avec inversion de l'index j
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.tbz[i, 91 - j] = fortran_float(vals[j])
            for i in range(d.nr):
                vals = f21.read_tokens(d.ntet)
                for j in range(d.ntet):
                    d.d2bz[i, 91 - j] = fortran_float(vals[j])
            try:
                for i in range(d.nr):
                    vals = f21.read_tokens(d.ntet)
                    for j in range(d.ntet):
                        d.dbzxlu[i, j] = fortran_float(vals[j])
                for i in range(d.nr):
                    vals = f21.read_tokens(d.ntet)
                    for j in range(d.ntet):
                        d.dbzylu[i, j] = fortran_float(vals[j])
            except Exception:
                print('attention fichier magnetique incomplet...')
                print('bx et by seront nuls...')
                print('vous resterez dans le plan median...')
                d.dbzxlu[:d.nr, :d.ntet] = 0
                d.dbzylu[:d.nr, :d.ntet] = 0

            # Symétries de la carte magnétique
            nr = d.nr
            d.tbz[:nr, 180:91:-1] = d.tbz[:nr, :89]
            d.tbz[:nr, 181:361] = d.tbz[:nr, 1:181]
            d.tbz[:nr, 361:721] = d.tbz[:nr, 1:361]
            d.d2bz[:nr, 180:91:-1] = d.d2bz[:nr, :89]
            d.d2bz[:nr, 181:361] = d.d2bz[:nr, 1:181]
            d.d2bz[:nr, 361:721] = d.d2bz[:nr, 1:361]

            r2 = np.sqrt(2.0) / 2.0
            d.dbzx[:nr, 90::-1] = r2 * (d.dbzxlu[:nr, :91] + d.dbzylu[:nr, :91])
            d.dbzy[:nr, 90::-1] = r2 * (d.dbzxlu[:nr, :91] - d.dbzylu[:nr, :91])
            d.dbzx[:nr, 180:91:-1] = d.dbzy[:nr, :89]
            d.dbzy[:nr, 180:91:-1] = d.dbzx[:nr, :89]
            d.dbzx[:nr, 181:361] = -d.dbzx[:nr, 179::-1]
            d.dbzy[:nr, 181:361] = d.dbzy[:nr, 179::-1]
            d.dbzx[:nr, 361:721] = d.dbzx[:nr, 359::-1]
            d.dbzy[:nr, 361:721] = -d.dbzy[:nr, 359::-1]
            d.htet = d.pi / 4.0 / (d.ntet - 1.0)
            ntetall = 721

    if d.ncoh[0] != 0 or d.ncoh[1] != 0:
        print('carte mag. + rajout harmonique 1 et/ou 2...')
        addh(d.rmin, d.rmax, d.nr, ntetall, d.tbz, d.d2bz,
             d.dbzx, d.dbzy, d.ncoh, d.coh, d.phihar)

    d.hr = (d.rmax - d.rmin) / (d.nr - 1.0)
    d.hrm1 = 1.0 / d.hr
    d.htetm1 = 1.0 / d.htet
    d.hrhtet = d.hr * d.htet
    d.hrhtetm1 = 1.0 / d.hrhtet

    # Sans b_factor, la carte est utilisee telle quelle (comme le Fortran).
    # Avec b_factor, elle est ramenee a max|Bz| = b_factor (T).
    if d.b_factor is not None:
        bz_max = np.max(np.abs(d.tbz[:d.nr, :ntetall]))
        if bz_max == 0.0:
            raise ValueError(f'carte magnetique nulle ({d.ficb}) : b_factor inapplicable')
        scale = d.b_factor / bz_max
        d.tbz *= scale
        d.d2bz *= scale
        d.dbzx *= scale
        d.dbzy *= scale
        print(f'carte magnetique : Bmax_lu = {bz_max:.6e} T, facteur B = {d.b_factor:.6f} T')

    # Correction relative dB/B
    d.tbz *= (1.0 + d.dbsb)
    d.d2bz *= (1.0 + d.dbsb)

    # Lecture potentiel centre
    if d.ficv1 != 'no':
        if d.repswap == 'y':
            fdtype = np.dtype(np.float64).newbyteorder('>')
            idtype = np.dtype(np.int32).newbyteorder('>')
        else:
            fdtype = np.float64
            idtype = np.int32
        with open(d.ficv1, 'rb') as f22:
            # Lecture entête Fortran (record marker)
            rec = f22.read(4)
            header_floats = np.frombuffer(f22.read(32), dtype=fdtype)
            header_ints = np.frombuffer(f22.read(8), dtype=idtype)
            d.xmin1 = float(header_floats[0]) / 1000.0
            d.xmax1 = float(header_floats[1]) / 1000.0
            d.ymin1 = float(header_floats[2]) / 1000.0
            d.ymax1 = float(header_floats[3]) / 1000.0
            nx1 = int(header_ints[0])
            ny1 = int(header_ints[1])
            f22.read(4)  # fin record
            for i in range(nx1):
                f22.read(4)
                row = np.frombuffer(f22.read(ny1 * 8), dtype=fdtype)
                d.pot1[i, :ny1] = row
                f22.read(4)
            try:
                for i in range(nx1):
                    f22.read(4)
                    row = np.frombuffer(f22.read(ny1 * 8), dtype=fdtype)
                    d.d2pot1[i, :ny1] = row
                    f22.read(4)
            except Exception:
                print('attention fichier potentiel centre incomplet...')
                print('ez sera nul...')
                d.d2pot1[:] = 0.0
    else:
        print('pas de fichier centre cyclo...')
        vals = parse_floats(read_input('xmin1 xmax1 ymin1 ymax1 nx1 ny1 ?'))
        d.xmin1, d.xmax1, d.ymin1, d.ymax1 = vals[:4]
        nx1 = int(vals[4])
        ny1 = int(vals[5])
        d.pot1[:] = 0.0
        d.d2pot1[:] = 0.0

    d.hx1 = (d.xmax1 - d.xmin1) / (nx1 - 1.0)
    d.hy1 = (d.ymax1 - d.ymin1) / (ny1 - 1.0)
    d.hx1hy1 = d.hx1 * d.hy1
    d.hx1m1 = 1.0 / d.hx1
    d.hy1m1 = 1.0 / d.hy1
    d.hx1hy1m1 = 1.0 / d.hx1hy1

    # Lecture potentiel périphérique
    with open(d.ficv2, 'r') as f23_raw:
      f23 = FortranFile(f23_raw)
      with f23:
        vals = f23.read_tokens(6)
        d.xmin2 = fortran_float(vals[0])
        d.xmax2 = fortran_float(vals[1])
        d.ymin2 = fortran_float(vals[2])
        d.ymax2 = fortran_float(vals[3])
        nx2 = int(fortran_float(vals[4]))
        ny2 = int(fortran_float(vals[5]))
        for i in range(nx2):
            vals = f23.read_tokens(ny2)
            for j in range(ny2):
                d.pot2[i, j] = fortran_float(vals[j])
        for i in range(nx2):
            vals = f23.read_tokens(ny2)
            for j in range(ny2):
                d.d2pot2[i, j] = fortran_float(vals[j])

    d.hx2 = (d.xmax2 - d.xmin2) / (nx2 - 1.0)
    d.hy2 = (d.ymax2 - d.ymin2) / (ny2 - 1.0)
    d.hx2hy2 = d.hx2 * d.hy2
    d.hx2m1 = 1.0 / d.hx2
    d.hy2m1 = 1.0 / d.hy2
    d.hx2hy2m1 = 1.0 / d.hx2hy2

    # Distribution des particules (iflcor) - APRES chargement champ magnetique
    if d.iflcor != 0:
        m = 6
        boite(d)
        calcb(d)
        for k in range(d.nbrfam):
            nbr = d.nbrion[k]
            if nbr <= 1:
                continue  # une seule particule, pas de distribution a appliquer
            ps = remp(m, nbr)
            d.dxs[:, :nbr] = ps
            rcor = choles(m, sigma_matrices[k])
            trcor = rcor.T
            d.dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    d.dxtmp[i, :nbr] += trcor[i, j] * d.dxs[j, :nbr]
            d.dxs[:6, :nbr] = d.dxtmp[:6, :nbr]

            i1 = d.index[k]
            i2 = i1 + d.nbrion[k]
            rc = -d.am0c[k] * d.p[i1] / d.q[i1] / d.bz[i1]

            if d.r0[k] > d.rmax1:
                print('prise en compte de sqrt(cos(phi))...')
                d.vec4[:nbr] = np.sqrt(np.cos(d.harm * d.dxs[4, :nbr] / rc))
                d.dxs[0, :nbr] = d.dxs[0, :nbr] * d.vec4[:nbr] - rc * (1.0 - d.vec4[:nbr])
            else:
                print('pas de sqrt(cos(phi))...')

            d.vec1[:nbr] = np.cos(d.dxs[4, :nbr] / rc)
            d.vec2[:nbr] = np.sin(d.dxs[4, :nbr] / rc)
            d.vec3[:nbr] = np.sin(d.dxs[4, :nbr] / rc / 2.0)**2
            d.dxd[0, :nbr] = d.dxs[0, :nbr] - 2.0 * (rc + d.dxs[0, :nbr]) * d.vec3[:nbr]
            d.dxs[5, :nbr] = np.sqrt(
                (d.p[i1] + d.dxs[5, :nbr])**2 -
                d.dxs[1, :nbr]**2 - d.dxs[3, :nbr]**2) - d.p[i1]
            d.dxd[4, :nbr] = (rc + d.dxs[0, :nbr]) * d.vec2[:nbr]
            d.dxd[1, :nbr] = d.dxs[1, :nbr] * d.vec1[:nbr] - d.dxs[5, :nbr] * d.vec2[:nbr]
            d.dxd[5, :nbr] = d.dxs[1, :nbr] * d.vec2[:nbr] + d.dxs[5, :nbr] * d.vec1[:nbr]
            d.dxd[2:4, :nbr] = d.dxs[2:4, :nbr]

            co = np.cos(d.teta0[k] - d.pr0[k])
            si = np.sin(d.teta0[k] - d.pr0[k])
            d.rot[:] = 0.0
            d.rot[0, 0] = co
            d.rot[0, 4] = -si
            d.rot[1, 1] = co
            d.rot[1, 5] = -si
            d.rot[2, 2] = 1.0
            d.rot[3, 3] = 1.0
            d.rot[4, 0] = si
            d.rot[4, 4] = co
            d.rot[5, 1] = si
            d.rot[5, 5] = co
            d.dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    d.dxtmp[i, :nbr] += d.rot[i, j] * d.dxd[j, :nbr]
            d.dxd[:6, :nbr] = d.dxtmp[:6, :nbr]

            d.x[i1:i2] += d.dxd[0, :nbr]
            d.y[i1:i2] += d.dxd[4, :nbr]
            d.vec[:nbr] = d.px[i1:i2] * d.vec1[:nbr] - d.py[i1:i2] * d.vec2[:nbr]
            d.py[i1:i2] = d.px[i1:i2] * d.vec2[:nbr] + d.py[i1:i2] * d.vec1[:nbr] + d.dxd[5, :nbr]
            d.px[i1:i2] = d.vec[:nbr] + d.dxd[1, :nbr]
            if d.espace[:2] in ('3d', '6d'):
                d.z[i1:i2] += d.dxd[2, :nbr]
                d.pz[i1:i2] += d.dxd[3, :nbr]
            if d.r0[k] > d.rmax1:
                d.px[i1:i2] *= d.vec4[:nbr]
                d.py[i1:i2] *= d.vec4[:nbr]
                d.pz[i1:i2] *= d.vec4[:nbr]

        if d.ficloa != 'no':
            print('attention : chargement paquet d ions...')
            with open(d.ficloa, 'rb') as f26:
                for i in range(d.n):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.x[i], d.y[i], d.z[i] = data
                for i in range(d.n):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.px[i], d.py[i], d.pz[i] = data

    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.tetnew = np.arccos(d.x[0] / rr)
    if d.y[0] < 0.0:
        d.tetnew = 2.0 * d.pi - d.tetnew
    d.tetnew = 180.0 / d.pi * d.tetnew
    d.numtou = 0
    if d.npaspt < 0:
        d.numtou = -d.nbdt // d.npaspt

    if d.ficz != 'no':
        files_dict[15] = open(d.ficz, 'w')
    if d.ficxy != 'no':
        files_dict[16] = open(d.ficxy, 'w')


# =============================================================================
# info - Affichage des informations
# =============================================================================
def info(d, iprinf_holder):
    """Affiche les informations de suivi."""

    if iprinf_holder[0] == 1:
        iprinf_holder[0] = 0
        d.ga0th = d.ga[0]
        d.w0th = d.w0[0]
        d.w0mevth = d.w0th / 1.60217733e-13 / d.a[0]
        v0th = d.clum * np.sqrt((d.ga0th - 1.0) * (d.ga0th + 1.0)) / d.ga0th
        d.r0th = v0th / d.omr
        d.b0th = -v0th * d.am0[0] * d.ga0th / d.q[0] / d.r0th
        d.ga1th = 1.0 / np.sqrt(1.0 - (d.ga0th**2 - 1.0) *
                                 d.r1th**2 / d.ga0th**2 / d.r0th**2)
        d.w1th = d.am0c2[0] * (d.ga1th - 1.0)
        d.w1mevth = d.w1th / 1.6021773e-13 / d.a[0]
        d.b1th = d.b0th / d.ga0th * d.ga1th
        potth = d.am0c2[0] * (d.ga1th - d.ga0th) / 4.0 / d.ntours / d.q[0]

        if d.ficz != 'no' and 15 in files_dict:
            f15 = files_dict[15]
            f15.write('Programme LIONS\n')
            f15.write('informations selon z\n')

        if d.ficxy != 'no' and 16 in files_dict:
            f16 = files_dict[16]
            f16.write('Programme LIONS\n')
            f16.write(f'{d.corps}\n')
            f16.write(f' masse atomique        :    {d.amat[0]:11.7f}\n')
            f16.write(f' A                     :          {d.a[0]:5.1f}\n')
            f16.write(f' q                     :          {d.qn[0]:5.1f}\n')
            r0mm = d.r0[0] * 1000.0
            teta0d = d.teta0[0] * 180.0 / d.pi
            z0mm = d.z0[0] * 1000.0
            pr0m = d.pr0[0] * 1000.0
            pz0m = d.pz0[0] * 1000.0
            f16.write(f' r0     : {r0mm:6.1f}\n')
            f16.write(f' teta0  : {teta0d:6.1f}\n')
            f16.write(f' z0     : {z0mm:6.1f}\n')
            f16.write(f' pr0    : {pr0m:8.2f}\n')
            f16.write(f' pz0    : {pz0m:8.2f}\n')
            f16.write(f' w0     : {d.w0mev1:14.10f}\n')
            f16.write(f' phi0   : {d.phi0[0]:11.2f}\n')
            f16.write(f' fhf    : {d.fhf:9.0f}\n')
            f16.write(f' harm   : {d.harm:9.0f}\n')
            f16.write(f' pot_d  : {-d.pothfd:9.1f}\n')
            f16.write(f' pot_g  : {-d.pothfg:9.1f}\n')
            f16.write(f' db/b   : {d.dbsb:9.6f}\n')
            if d.b_factor is not None:
                f16.write(f' B fact : {d.b_factor:9.6f} T\n')
            dtnano = d.dt * 1e9
            f16.write(f' nb part: {d.n:14d}\n')
            f16.write(f' nb iter: {d.nbdt:14d}\n')
            f16.write(f' nb i/t : {d.npaspt:14d}\n')
            f16.write(f' dt(ns) : {dtnano:9.6f}\n')
            f16.write('    ite   tr   teta      x       y      xc'
                      '      yc        bz '
                      '    es      er           r      tpr'
                      '       w           dphi\n')

    d.tetold = d.tetnew
    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.tetnew = np.arccos(np.clip(d.x[0] / rr, -1.0, 1.0))
    if d.y[0] < 0.0:
        d.tetnew = 2.0 * d.pi - d.tetnew
    d.tetnew = 180.0 / d.pi * d.tetnew

    if d.npaspt > 0:
        if d.tetnew < d.tetold:
            d.numtou += 1
    else:
        if d.tetnew > d.tetold:
            d.numtou -= 1

    d.tetall = d.tetnew + 360.0 * (d.numtou - 1)

    if (d.ite % d.ifreq != 0 and d.ite != 1 and
            d.ite != d.nbdt and d.ite != (d.nbdt - 1)):
        return

    phihf = 180.0 / d.pi * d.omhf * d.t
    dph = phihf - d.harm * d.tetall - 90.0
    wmevpa = d.am0c2[0] * (d.ga[0] - 1.0) / 1.60217733e-13 / d.a[0]
    xmm = d.x[0] * 1000.0
    ymm = d.y[0] * 1000.0
    zmm = d.z[0] * 1000.0
    rmm = rr * 1000.0
    pp = np.sqrt(d.px[0]**2 + d.py[0]**2)
    ppp = np.sqrt(d.px[0]**2 + d.py[0]**2 + d.pz[0]**2)
    pr_val = (xmm * d.px[0] + ymm * d.py[0]) / rmm
    tpr = np.arcsin(pr_val / ppp) * 1000.0
    xc1mm = (d.x[0] + d.py[0] * d.dt / d.qdtsm0c[0] / d.bz[0]) * 1000.0
    yc1mm = (d.y[0] - d.px[0] * d.dt / d.qdtsm0c[0] / d.bz[0]) * 1000.0
    eskvcm = (d.ex[0] * d.px[0] + d.ey[0] * d.py[0]) / pp / 100000.0
    erkvcm = (d.ex[0] * d.py[0] - d.ey[0] * d.px[0]) / pp / 100000.0
    ezkvcm = d.ez[0] / 100000.0
    tpz = np.arcsin(d.pz[0] / ppp) * 1000.0

    if d.ite % d.ifreq != 0 and d.ite != 1:
        pass  # skip normal output, go to end-of-run output
    else:
        if d.ficz != 'no' and 15 in files_dict:
            if abs(d.z[0]) > d.coupnuz:
                anuz_c = (-1.0 / d.am0[0] / d.ga[0] * (d.harm / d.omhf)**2 / d.z[0] *
                          d.q[0] * (d.ez[0] + d.clum / d.ga[0] *
                          (d.px[0] * d.by[0] - d.py[0] * d.bx[0])))
                d.anuz_cyclo = np.sqrt(abs(anuz_c))
            else:
                d.anuz_cyclo = 0.0
            files_dict[15].write(
                f' {d.ite:6d}  {d.numtou:3d}  {d.tetnew:5.1f}'
                f' {xmm:8.2f} {ymm:8.2f} {zmm:8.2f} {rmm:8.2f}'
                f'   {d.bx[0]:7.4f} {d.by[0]:7.4f} {d.bz[0]:7.4f}'
                f' {ezkvcm:7.3f}     {tpz:6.1f} {d.anuz_cyclo:12.5f}'
                f' {d.pot_ver_vrai:12.5f}\n')

        if d.ficxy != 'no' and 16 in files_dict:
            files_dict[16].write(
                f' {d.ite:7d}  {d.numtou:3d}  {d.tetnew:5.1f}'
                f' {xmm:7.1f} {ymm:7.1f} {xc1mm:7.1f} {yc1mm:7.1f}'
                f'   {d.bz[0]:7.4f} {eskvcm:7.3f} {erkvcm:7.3f}'
                f'     {rmm:7.1f}  {tpr:6.1f} {wmevpa:9.5f} {dph:15.2f}\n')

    # Corrélations (sautees avec les sorties z/xy, comme le goto 800 du Fortran)
    if ((d.ite % d.ifreq == 0 or d.ite == 1) and d.iflcor != 0
            and not (d.ficcor == 'no' and d.ficmat == 'no' and d.ficdc == 'no')):
        for k in range(d.nbrfam):
            nbr = d.nbrion[k]
            i1 = d.index[k]
            i2 = i1 + nbr
            d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
            rc = -d.am0c[k] * d.p[i1] / d.q[i1] / d.bz[i1]
            # Calcul phase space (simplifié, même logique que caldxs)
            d.dxd[0, :nbr] = d.x[i1:i2] - d.x[i1]
            d.dxd[4, :nbr] = d.y[i1:i2] - d.y[i1]
            ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
            tetai1 = np.arccos(np.clip(d.x[i1] / ri1, -1.0, 1.0))
            if d.y[i1] < 0.0:
                tetai1 = 2.0 * d.pi - tetai1
            pri1 = (d.x[i1] * d.px[i1] + d.y[i1] * d.py[i1]) / ri1
            tpri1 = np.arcsin(pri1 / d.p[i1])
            co = np.cos(tetai1 - tpri1)
            si = np.sin(tetai1 - tpri1)
            d.vec[:nbr] = d.dxd[0, :nbr] * co + d.dxd[4, :nbr] * si
            d.dxd[4, :nbr] = -d.dxd[0, :nbr] * si + d.dxd[4, :nbr] * co
            d.dxd[0, :nbr] = d.vec[:nbr]
            for ii in range(nbr):
                d.dxs[4, ii] = rc * np.arctan(d.dxd[4, ii] / (rc + d.dxd[0, ii]))
            d.vec[:nbr] = np.sin(d.dxs[4, :nbr] / 2.0 / rc)**2
            d.dxs[0, :nbr] = (d.dxd[0, :nbr] + 2.0 * rc * d.vec[:nbr]) / (1.0 - 2.0 * d.vec[:nbr])
            d.vec1[:nbr] = np.cos(d.dxs[4, :nbr] / rc)
            d.vec2[:nbr] = np.sin(d.dxs[4, :nbr] / rc)
            d.dxd[1, :nbr] = d.px[i1:i2] - d.px[i1] * d.vec1[:nbr] + d.py[i1] * d.vec2[:nbr]
            d.dxd[5, :nbr] = d.py[i1:i2] - d.px[i1] * d.vec2[:nbr] - d.py[i1] * d.vec1[:nbr]
            d.vec[:nbr] = d.dxd[1, :nbr] * co + d.dxd[5, :nbr] * si
            d.dxd[5, :nbr] = -d.dxd[1, :nbr] * si + d.dxd[5, :nbr] * co
            d.dxd[1, :nbr] = d.vec[:nbr]
            d.dxs[1, :nbr] = d.dxd[1, :nbr] * d.vec1[:nbr] + d.dxd[5, :nbr] * d.vec2[:nbr]
            d.dxs[5, :nbr] = -d.dxd[1, :nbr] * d.vec2[:nbr] + d.dxd[5, :nbr] * d.vec1[:nbr]
            d.dxs[2, :nbr] = d.z[i1:i2] - d.z[i1]
            d.dxs[3, :nbr] = d.pz[i1:i2] - d.pz[i1]
            d.dxs[5, :nbr] = np.sqrt(
                (d.p[i1] + d.dxs[5, :nbr])**2 +
                d.dxs[1, :nbr]**2 + d.dxs[3, :nbr]**2) - d.p[i1]

            d.cor[:] = 0.0
            for i in range(6):
                for j in range(6):
                    d.cor[i, j] = np.sum(d.dxs[i, :nbr] * d.dxs[j, :nbr] *
                                         d.ok[:nbr].astype(float)) / d.nbrok

            tetai1_deg = tetai1 * 180.0 / d.pi
            ri1mm = ri1 * 1000.0
            d.t1[:] = 0.0
            d.t1[0, 0] = -1.0
            d.t1[1, 1] = -1.0
            for i in range(2, 6):
                d.t1[i, i] = 1.0
            d.t1[0, 5] = ri1 / d.p[i1]
            d.t1[4, 1] = -d.t1[0, 5]
            d.cordc = d.t1 @ d.cor @ d.t1.T

            if d.ficcor != 'no' and 17 in files_dict:
                gai1 = np.sqrt(1.0 + d.p[i1]**2)
                dr = np.sqrt(d.cor[0, 0]) * d.rac8
                drmm = dr * 1000.0
                ds = np.sqrt(d.cor[4, 4]) * d.rac8
                dsmm = ds * 1000.0
                dz_val = np.sqrt(d.cor[2, 2]) * d.rac8
                dzmm = dz_val * 1000.0
                dphi = ds * d.omhf / d.clum / d.p[i1] * 180.0 / d.pi * gai1
                dpr = np.sqrt(d.cor[1, 1]) * d.rac8
                drdpr = 8.0 * np.sqrt(abs(d.cor[0, 0] * d.cor[1, 1] - d.cor[0, 1]**2))
                epsr = drdpr * 1e6 / d.p[i1]
                dzdpz = 8.0 * np.sqrt(abs(d.cor[2, 2] * d.cor[3, 3] - d.cor[2, 3]**2))
                epsz = dzdpz * 1e6 / d.p[i1]
                dp = np.sqrt(d.cor[5, 5]) * d.rac8
                dsdp = 8.0 * np.sqrt(abs(d.cor[4, 4] * d.cor[5, 5] - d.cor[4, 5]**2))
                dwsw = dp * d.p[i1] / gai1 / (gai1 - 1.0)
                wmpa = d.am0c2[k] * (gai1 - 1.0) / 1.60217733e-13 / d.a[k]
                dwmpa = dwsw * wmpa
                dwdphi = (d.omhf / d.clum / d.p[i1] * 180.0 / d.pi * gai1 *
                          d.p[i1] / gai1 / (gai1 - 1.0) * wmpa * dsdp * 1000.0)
                files_dict[17].write(
                    f' {d.numtou:3d}  {ri1mm:7.2f}  {tetai1_deg:5.1f}'
                    f' {drmm:6.2f} {dsmm:6.2f} {dphi:7.2f}'
                    f' {epsr:6.2f} {dwsw:12.7e} {wmpa:6.3f}'
                    f' {dwmpa:12.7e} {dwdphi:12.7e}'
                    f'  {dzmm:6.2f} {epsz:6.2f}\n')

    # Sauvegarde fin de simulation
    if d.ite == (d.nbdt - 1):
        d.xf = d.x[0]
        d.yf = d.y[0]
        d.zf = d.z[0]
        d.rf = rr
        d.tetf = d.tetnew
        d.phif = phihf
        if d.ficsto != 'no':
            print('attention : sauvegarde paquet d ions...')
            with open(d.ficsto, 'wb') as f26:
                for i in range(d.n):
                    f26.write(struct.pack('ddd', d.x[i], d.y[i], d.z[i]))

    if d.ite == d.nbdt:
        d.pxf = d.px[0]
        d.pyf = d.py[0]
        pzf = d.pz[0]
        d.pf = np.sqrt(d.pxf**2 + d.pyf**2 + pzf**2)
        copzf = np.cos(pzf)
        alphaf = np.arccos(d.pxf / d.pf / copzf)
        if d.pyf < 0.0:
            alphaf = 2.0 * d.pi - alphaf
        d.tprf = ((d.pi / 2.0 + d.tetf * d.pi / 180.0) - alphaf) * 1000.0
        d.wf = wmevpa
        d.tpzf = np.arcsin(pzf / d.pf) * 1000.0
        if d.ficsto != 'no':
            with open(d.ficsto, 'ab') as f26:
                for i in range(d.n):
                    f26.write(struct.pack('ddd', d.px[i], d.py[i], d.pz[i]))

        if d.ficxy != 'no' and 16 in files_dict:
            f16 = files_dict[16]
            f16.write('---------------------------------\n')
            f16.write('info pour le mouvement inverse...\n')
            f16.write(f'r0 teta0 z0 pr0 pz0 w0mev phi0 :\n')
            f16.write(f'{d.rf} {d.tetf} {d.zf} {d.tprf}\n')
            f16.write(f'{d.tpzf} {d.wf} {d.phif}\n')
            f16.write('info pour continuer le mouvement ...\n')
            f16.write(f'{rr} {d.tetnew} {d.z[0]} {d.tprf}\n')
            f16.write(f'{d.tpzf} {d.wf} {phihf}\n')
            f16.write('info utile pour DEDUCE_SIGMA ...\n')
            f16.write(f'amat a qn w0mev : \n')
            f16.write(f'{d.amat[0]} {d.a[0]} {d.qn[0]} {d.wf}\n')
            f16.write(f'champ magnetique (val. absolue) :\n{-d.bz[0]}\n')
            f16.write(f'fhf harmonique :\n{d.fhf} {d.harm}\n')


# =============================================================================
# initgra - Initialisation graphique (bidims)
# =============================================================================
def initgra(d, bidim, numbid=None):
    """Initialise les diagnostics graphiques (bidims). numbid=None : demande a l'utilisateur."""
    d.fenlog[:] = True
    d.ifreb = 1
    if numbid is None:
        numbid = int(read_input('numero premier bidim (0 si aucun bidim) :'))
    d.numbid = numbid
    if d.numbid == 0:
        return

    d.ifreb = int(read_input('frequence_de_trace par rapport aux pas de temps:'))
    d.baryc = read_input('trace des centres ou des barycentres (c ou b) :')
    d.ngravm = abs(d.npaspt // d.ifreb)
    d.ngrav = 0

    vals = parse_floats(read_input(
        'boite des bidims (xmin xmax ymin ymax en metres) :'))
    xbid = np.zeros((2, 2))
    xbid[0, 0] = vals[0]
    xbid[0, 1] = vals[1]
    xbid[1, 0] = vals[2]
    xbid[1, 1] = vals[3]

    nbx = np.full(20, 13, dtype=np.int32)
    nby = np.full(20, 13, dtype=np.int32)
    titre = [''] * 20
    labelx = [''] * 20
    labely = [''] * 20
    bornes = np.zeros((2, 2, 50), dtype=np.float64)

    titre[0] = 'trajectoire ion'
    labelx[0] = 'x'
    labely[0] = 'y'
    if d.baryc == 'c':
        titre[1] = 'centres de courbure'
    else:
        titre[1] = 'barycentres'
    labelx[1] = 'x'
    labely[1] = 'y'

    d.itebid = 1
    if d.n > 1:
        titre[2] = 'trajectoires ions'
        titre[3] = 'correlation dphi-dw/w'
        titre[4] = 'correlation r-dpr/p0'
        vals = parse_floats(read_input('bidim dphi-dw/w : dphi_max dw/w_max :'))
        dphimax = vals[0]
        dwswmax = vals[1]
        vals = parse_floats(read_input('bidim dr-dpr/p0 : dr_max dpr/p0_max :'))
        drmax = vals[0]
        dprsp0max = vals[1]
        vals = parse_floats(read_input('bidim dz-dpz/p0 : dz_max dpz/p0_max :'))
        dzmax = vals[0]
        dpzsp0max = vals[1]
        print(f'itemax : {d.nbdt}')
        d.itebid = int(read_input('debut_iteration_stockage_bidim (entre 1 et itemax) :'))

    vals = parse_floats(read_input('valeurs de la fenetre des barycentres en metres:'))
    xbamin, xbamax, ybamin, ybamax = vals[:4]

    bornes[0, 0, 0] = xbid[0, 0]
    bornes[1, 0, 0] = xbid[0, 1]
    bornes[0, 1, 0] = xbid[1, 0]
    bornes[1, 1, 0] = xbid[1, 1]
    bornes[0, 0, 1] = xbamin
    bornes[1, 0, 1] = xbamax
    bornes[0, 1, 1] = ybamin
    bornes[1, 1, 1] = ybamax

    if d.n > 1:
        nbr_bid = 17
    else:
        nbr_bid = 5

    ier = bidim.openbi(nbr_bid, d.numbid, nbx, nby, titre, labelx, labely, bornes)
    if ier != 0:
        print(f'bidim : openbi : stop erreur : {ier}')
        sys.exit(1)

    d.ax1[0] = 8191.0 / (xbid[0, 1] - xbid[0, 0])
    d.ax0[0] = -d.ax1[0] * xbid[0, 0]
    d.ay1[0] = 8191.0 / (xbid[1, 1] - xbid[1, 0])
    d.ay0[0] = -d.ay1[0] * xbid[1, 0]
    d.ax1[1] = 8191.0 / (xbamax - xbamin)
    d.ax0[1] = -d.ax1[1] * xbamin
    d.ay1[1] = 8191.0 / (ybamax - ybamin)
    d.ay0[1] = -d.ay1[1] * ybamin

    if d.n > 1:
        d.ax1[3] = 8191.0 / 2.0 / dphimax
        d.ax0[3] = d.ax1[3] * dphimax
        d.ay1[3] = 8191.0 / 2.0 / dwswmax
        d.ay0[3] = d.ay1[3] * dwswmax
        d.ax1[4] = 8191.0 / 2.0 / drmax
        d.ax0[4] = d.ax1[4] * drmax
        d.ay1[4] = 8191.0 / 2.0 / dprsp0max
        d.ay0[4] = d.ay1[4] * dprsp0max
        d.ax1[5] = 8191.0 / 2.0 / drmax
        d.ax0[5] = d.ax1[5] * drmax
        d.ay1[5] = 8191.0 / 2.0 / dphimax
        d.ay0[5] = d.ay1[5] * dphimax
        d.ax1[6] = 8191.0 * d.ifreb / d.nbdt
        d.ax0[6] = 0.0
        d.ay1[6] = 8191.0 / 2.0 / drmax
        d.ay0[6] = d.ay1[6] * drmax
        d.ax1[7] = 8191.0 / 2.0 / dzmax
        d.ax0[7] = d.ax1[7] * dzmax
        d.ay1[7] = 8191.0 / 2.0 / dpzsp0max
        d.ay0[7] = d.ay1[7] * dpzsp0max
        d.ax1[8] = 8191.0 * d.ifreb / d.nbdt
        d.ax0[8] = 0.0
        d.ay1[8] = 8191.0 / 2.0 / dzmax
        d.ay0[8] = d.ay1[8] * dzmax

    # Lecture repères
    ficrep = read_input('lecture du fichier_des_reperes (lignes, piliers etc...) :')
    d.nbrpie = 0
    if ficrep != 'no':
        try:
            with open(ficrep, 'r') as f24:
                vals = f24.readline().split()
                xr1, yr1, xr2, yr2 = fortran_float(vals[0]), fortran_float(vals[1]), fortran_float(vals[2]), fortran_float(vals[3])
                vals = f24.readline().split()
                vals = f24.readline().split()  # niter, omer
                nbpo = 360
                while True:
                    typcou = f24.readline().strip()
                    if typcou == 'fintyp':
                        break
                    elif typcou == 'droite':
                        nbrdro = int(f24.readline().split()[0])
                        for num in range(nbrdro):
                            vals = f24.readline().split()
                            xd1, yd1, xd2, yd2 = [float(v) for v in vals[:4]]
                            sxd = (xd2 - xd1) / (nbpo - 1.0)
                            syd = (yd2 - yd1) / (nbpo - 1.0)
                            for kk in range(nbpo):
                                xd = xd1 + kk * sxd
                                yd = yd1 + kk * syd
                                ix = int(d.ax1[0] * xd / 1000.0 + d.ax0[0])
                                iy = int(d.ay1[0] * yd / 1000.0 + d.ay0[0])
                                bidim.incbid(d.numbid, ix, iy)
                    elif typcou == 'cercle':
                        nbrcer = int(f24.readline().split()[0])
                        for num in range(nbrcer):
                            vals = f24.readline().split()
                            xce, yce, rce = fortran_float(vals[0]), fortran_float(vals[1]), fortran_float(vals[2])
                            tetce1 = fortran_float(vals[3]) * d.pi / 180.0
                            tetce2 = fortran_float(vals[4]) * d.pi / 180.0
                            npcer = int(vals[5])
                            flacer = fortran_float(vals[6])
                            dtetce = (tetce2 - tetce1) / (nbpo - 1.0)
                            if flacer < -10.0:
                                d.nbrpie += 1
                                d.xpie[d.nbrpie - 1] = xce / 1000.0
                                d.ypie[d.nbrpie - 1] = yce / 1000.0
                                d.rpie[d.nbrpie - 1] = rce / 1000.0
                            if flacer > -9000.0:
                                for kk in range(nbpo):
                                    tetce = tetce1 + kk * dtetce
                                    xd = xce + rce * np.cos(tetce)
                                    yd = yce + rce * np.sin(tetce)
                                    ix = int(d.ax1[0] * xd / 1000.0 + d.ax0[0])
                                    iy = int(d.ay1[0] * yd / 1000.0 + d.ay0[0])
                                    bidim.incbid(d.numbid, ix, iy)
                    else:
                        print(f'erreur lecture fichier reperes... typcou : {typcou}')
                        sys.exit(1)
        except Exception as e:
            print(f'erreur en lecture du fichier reperes: {e}')

    if d.nbrpie != 0:
        print(f'attention : il y a des pieges a ions...')
        print(f'nombre de pieges : {d.nbrpie}')
        d.itepiem = int(read_input('iteration max pour les pieges:'))
        print('flag de choix particules piegees: (0 ou 1)')
        print('0   : reste scotchee dans les pieges')
        print('1   : continue, confondue avec la particule centrale')
        d.iflpie = int(read_input(''))


# =============================================================================
# bid - Enregistrement des bidims
# =============================================================================
def bid(d, bidim):
    """Enregistre les données dans les bidims."""
    if d.ite < d.itebid:
        return

    for k in range(d.np_):
        if not d.fenlog[k]:
            continue
        ix = int(d.ax1[0] * d.x[k] + d.ax0[0])
        iy = int(d.ay1[0] * d.y[k] + d.ay0[0])
        if 0 <= ix <= 8191 and 0 <= iy <= 8191:
            if k == 0:
                bidim.incbid(d.numbid, ix, iy)
                if d.np_ > 1:
                    bidim.incbid(d.numbid + 2, ix, iy)
            else:
                bidim.incbid(d.numbid + 2, ix, iy)
        if k == 0:
            xck = d.x[k] + d.py[k] * d.dt / d.qdtsm0c[k] / d.bz[k]
            yck = d.y[k] - d.px[k] * d.dt / d.qdtsm0c[k] / d.bz[k]
            if d.baryc == 'c':
                ix2 = int(d.ax1[0] * xck + d.ax0[0])
                iy2 = int(d.ay1[0] * yck + d.ay0[0])
                if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                    bidim.incbid(d.numbid + 1, ix2, iy2)
            if d.baryc == 'b':
                d.ngrav += 1
                d.xgrav[d.ngrav - 1] = xck
                d.ygrav[d.ngrav - 1] = yck
                if d.ngrav == d.ngravm:
                    d.ngrav = 0
                    xbary = np.sum(d.xgrav[:d.ngravm]) / d.ngravm
                    ybary = np.sum(d.ygrav[:d.ngravm]) / d.ngravm
                    ix2 = int(d.ax1[1] * xbary + d.ax0[1])
                    iy2 = int(d.ay1[1] * ybary + d.ay0[1])
                    if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                        bidim.incbid(d.numbid + 1, ix2, iy2)

    if d.np_ > 1 and d.nbrfam == 1:
        caldxs(d)
        k = 0
        nbr = d.nbrion[k]
        i1 = d.index[k]
        gai1 = np.sqrt(1.0 + d.p[i1]**2)
        d.vec1[:nbr] = d.dxs[4, :nbr] * d.omhf / d.clum / d.p[i1] * 180.0 / d.pi * gai1
        d.vec2[:nbr] = d.dxs[5, :nbr] * d.p[i1] / gai1 / (gai1 - 1.0)
        d.vec3[:nbr] = d.dxs[0, :nbr] * 1000.0
        d.vec4[:nbr] = d.dxs[1, :nbr] / d.p[i1] * 1000.0
        d.vec5[:nbr] = d.dxs[2, :nbr] * 1000.0
        d.vec6[:nbr] = d.dxs[3, :nbr] / d.p[i1] * 1000.0

        for i in range(nbr):
            if not d.fenlog[i]:
                continue
            ix = int(d.ax1[3] * d.vec1[i] + d.ax0[3])
            iy = int(d.ay1[3] * d.vec2[i] + d.ay0[3])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 3, ix, iy)

            ix = int(d.ax1[4] * d.vec3[i] + d.ax0[4])
            iy = int(d.ay1[4] * d.vec4[i] + d.ay0[4])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 4, ix, iy)

            ix = int(d.ax1[7] * d.vec5[i] + d.ax0[7])
            iy = int(d.ay1[7] * d.vec6[i] + d.ay0[7])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 7, ix, iy)

        for i in range(nbr):
            if not d.fenlog[i]:
                continue
            ix = int(d.ax1[5] * d.vec3[i] + d.ax0[5])
            iy = int(d.ay1[5] * d.vec1[i] + d.ay0[5])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 5, ix, iy)

        ix = int(d.ax1[6] * d.ite / d.ifreb + d.ax0[6])
        ix = max(0, min(ix, 8191))
        for i in range(nbr):
            if not d.fenlog[i]:
                continue
            iy = int(d.ay1[6] * d.vec3[i] + d.ay0[6])
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 6, ix, iy)
            iy = int(d.ay1[8] * d.vec5[i] + d.ay0[8])
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.numbid + 8, ix, iy)


# =============================================================================
# Programme principal
# =============================================================================
def banner():
    print('----------------------------------------------------')
    print('LIONS       : Python (converti de Fortran 90)')
    print('auteur      : Patrick Bertrand - GANIL')
    print('version  1  : decembre 1998...decembre 2004')
    print('version  2  : mars 2017, modif seed, swap etc...')
    print('conversion Python : 2026')
    print('----------------------------------------------------')
    print()
    print('modif swap de mars 2017 :')
    print('on swap les octets ? (y/n) (y si anciens fichiers data unformatted...)')


def set_vertical_sorter(d, rmin=1000.0, rmax=2000.0, teta1=0.0, teta2=0.0,
                        pot=0.0, h=0.0, phi=0.0, anuz=0.0, coupnuz=0.0, dpotspot=0.0):
    """Parametres du trieur vertical (angles en degres). h=0 : pas de tri vertical."""
    d.rmin_ver = rmin
    d.rmax_ver = rmax
    d.teta1_ver = teta1 * np.pi / 180.0
    d.teta2_ver = teta2 * np.pi / 180.0
    d.pot_ver = pot
    d.h_ver = h
    d.phi_ver = phi * np.pi / 180.0
    d.anuz = anuz
    d.coupnuz = coupnuz
    d.dpotspot = dpotspot


def read_vertical_sorter(d):
    """Lit les parametres du trieur vertical sur l'entree standard."""
    print('--------------------')
    print('test trieur vertical :')
    print('--------------------')
    print('rmin_ver rmax_ver teta1_ver teta2_ver :')
    print('1000 2000 0 0 si pas de tri vertical ...')
    geom = parse_floats(input())
    print('pot_ver,h_ver,phi_ver anuz coupnuz dpotspot :')
    print('0 0 0 0 0 0 si pas de tri vertical...')
    set_vertical_sorter(d, *geom[:4], *parse_floats(input())[:6])
    print('-------------------')


def track(d, bidim=None, on_step=None):
    """Boucle de suivi : pas initial puis nbdt pas. Ferme les fichiers de sortie.

    on_step(d) est appele apres chaque pas ; s'il renvoie False, le suivi s'arrete.
    """
    d.ite = 0
    boite(d)
    calcbe(d)
    iprinf_holder = [1]
    info(d, iprinf_holder)
    try:
        if on_step is None or on_step(d) is not False:
            for ite in range(1, d.nbdt + 1):
                d.ite = ite
                piege(d)
                boite(d)
                calcbe(d)
                pousse(d)
                d.t = ite * d.dt + d.t0
                info(d, iprinf_holder)
                if d.numbid != 0 and (ite % d.ifreb == 0 or ite == 1):
                    bid(d, bidim)
                if on_step is not None and on_step(d) is False:
                    break
    finally:
        if d.numbid != 0:
            bidim.clobid()
        for f in files_dict.values():
            f.close()
        files_dict.clear()


def run(ficin, swap=False, vertical=None, on_step=None):
    """Simulation sans interaction ni bidim. Renvoie l'etat final (LionData).

    - ficin : fichier de donnees ; ses chemins sont relatifs au repertoire courant.
    - swap : inverser les octets du potentiel central (fichiers big-endian).
    - vertical : dict de parametres pour set_vertical_sorter, None = pas de tri.
    - on_step : voir track().
    """
    d = LionData()
    d.repswap = 'y' if swap else 'n'
    d.idebug = 0
    init(d, ficin)
    set_vertical_sorter(d, **(vertical or {}))
    initgra(d, None, numbid=0)
    track(d, on_step=on_step)
    return d


def write_trapped(d):
    """Propose d'ecrire les numeros des particules piegees."""
    if d.nbrpie != 0:
        nreste = np.count_nonzero(d.ok)
        print(f'nombre de particules non piegees : {nreste}')
        print('ecriture fichier des numero particules piegees ? (o/n)')
        reponse = input().strip()
        if reponse == 'o':
            print('choisissez nom du fichier : ')
            ficpie = input().strip()
            with open(ficpie, 'w') as f60:
                npi = 0
                for k in range(d.n):
                    if not d.ok[k]:
                        npi += 1
                        f60.write(f' {npi:5d} {k+1:5d}\n')


def lions():
    """Programme principal LIONS (interactif)."""
    banner()
    d = LionData()
    d.repswap = input().strip()
    d.idebug = 0
    init(d, read_input('fichier_de_donnees de depart :'))
    read_vertical_sorter(d)
    bidim = BidimManager()
    initgra(d, bidim)
    track(d, bidim)
    write_trapped(d)


if __name__ == '__main__':
    lions()
