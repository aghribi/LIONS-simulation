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
from lion_m import State

C_LIGHT = 299792458.0  # m/s
SQRT8 = np.sqrt(8.0)  # rms -> full width of a uniform 6D distribution


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
    d.in_box = ((d.r > d.map_rmin) & (d.r < (d.map_rmax - 0.02)) &
             (d.z > d.zmin) & (d.z < d.zmax))
    d.in_box = d.in_box & d.alive
    d.n_in_box = np.count_nonzero(d.in_box)


def box1(d):
    """Vérifie les particules dans la zone centrale."""
    d.in_central = d.in_box & (d.r < d.central_rmax) & (np.abs(d.y) < (d.central_ymax - d.central_dy))
    d.n_central = np.count_nonzero(d.in_central)


def box2(d, vecx, vecy):
    """Vérifie les particules dans la zone périphérique (vecx, vecy : repere du gap)."""
    d.in_outer = ((vecx > d.outer_xmin) & (vecx < d.outer_xmax) &
                (vecy > d.outer_ymin) & (vecy < d.outer_ymax) &
                d.in_box & ~d.in_central)
    d.n_outer = np.count_nonzero(d.in_outer)


def boxdf(d, i):
    """Vérifie les particules dans un déflecteur électrique."""
    # i est 0-based en Python
    cond_in = (((d.x - d.deflector_in_x[i]) * d.deflector_in_uy[i] - (d.y - d.deflector_in_y[i]) * d.deflector_in_ux[i]) <= 0.0)
    cond_out = (((d.x - d.deflector_out_x[i]) * d.deflector_out_uy[i] - (d.y - d.deflector_out_y[i]) * d.deflector_out_ux[i]) >= 0.0)

    if i > 0:
        d.in_element = d.extracted & cond_in & cond_out & d.alive
        d.n_in_element = np.count_nonzero(d.in_element)
        return

    if d.extraction_start_step == 0:
        d.in_element = (d.r > d.deflector_r[i]) & cond_in & cond_out
        d.in_element = d.in_element | (d.extracted & cond_in & cond_out)
        d.extracted = d.extracted | d.in_element
    else:
        d.in_element = (d.step >= d.extraction_start_step) & cond_in & cond_out
        d.in_element = d.in_element | (d.extracted & cond_in & cond_out)
        d.extracted = d.extracted | d.in_element

    d.extracted = d.extracted & d.alive
    d.in_element = d.in_element & d.alive
    d.n_in_element = np.count_nonzero(d.in_element)


def boxdm(d, i):
    """Vérifie les particules dans un dipôle magnétique."""
    cond_in = (((d.x - d.dipole_x[i]) * d.dipole_in_uy[i] - (d.y - d.dipole_y[i]) * d.dipole_in_ux[i]) <= 0.0)
    cond_out = (((d.x - d.dipole_x[i]) * d.dipole_out_uy[i] - (d.y - d.dipole_y[i]) * d.dipole_out_ux[i]) >= 0.0)
    d.in_element = d.extracted & cond_in & cond_out & d.alive
    d.n_in_element = np.count_nonzero(d.in_element)


def boxgex(d, i):
    """Vérifie les particules aux grilles d'extraction."""
    cond_cross = (((d.x - d.grid_x[i]) * d.grid_uy[i] -
                   (d.y - d.grid_y[i]) * d.grid_ux[i]) <= 0.0)

    if d.n_deflectors == 0 and d.n_dipoles == 0 and d.n_grids == 1:
        if d.step < d.extraction_start_step:
            d.crossed_grid[:] = False
            return
        d.crossed_grid = (d.grid_reached == 0) & d.alive & cond_cross
        reset_mask = (d.grid_reached == 1) & (
            ((d.x - d.grid_x[i]) * d.grid_uy[i] -
             (d.y - d.grid_y[i]) * d.grid_ux[i]) > 0.0)
        d.grid_reached[reset_mask] = 0
        return

    d.crossed_grid = d.extracted & d.alive & (d.grid_reached == i) & cond_cross


# =============================================================================
# piege - Piégeage des ions
# =============================================================================
def piege(d):
    """Gère le piégeage des ions dans les obstacles."""
    if d.n_traps == 0:
        return
    mask_z = (d.z < d.zmin) | (d.z > d.zmax)
    d.alive[mask_z] = False

    if d.step <= d.trap_last_step:
        for k in range(d.n_traps):
            rpie2 = d.trap_r[k]**2
            mask_pie = ((d.x - d.trap_x[k])**2 + (d.y - d.trap_y[k])**2) <= rpie2
            d.alive[mask_pie] = False

    d.n_alive = np.count_nonzero(d.alive)
    if d.trapped_follow == 1:
        not_ok = ~d.alive
        d.x[not_ok] = d.x[0]
        d.y[not_ok] = d.y[0]
        d.z[not_ok] = d.z[0]
        d.px[not_ok] = d.px[0]
        d.py[not_ok] = d.py[0]
        d.pz[not_ok] = d.pz[0]
        d.gamma[not_ok] = d.gamma[0]


# =============================================================================
# pousse - Boris pusher (avance des particules)
# =============================================================================
def pousse(d):
    """Avance les particules avec le Boris pusher."""
    if d.trapped_follow == 0:
        m = d.alive  # masque
        vec_m = d.boris_b_coef[m] / d.gamma[m]
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

        fx_m = d.boris_e_coef[m] * d.ex[m]
        fy_m = d.boris_e_coef[m] * d.ey[m]
        fz_m = d.boris_e_coef[m] * d.ez[m]

        v4 = vec_m * (fx_m + 2.0 * d.px[m])
        v5 = vec_m * (fy_m + 2.0 * d.py[m])
        v6 = vec_m * (fz_m + 2.0 * d.pz[m])
        d.px[m] = d.px[m] - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx_m
        d.py[m] = d.py[m] + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy_m
        d.pz[m] = d.pz[m] + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz_m

        d.gamma[m] = np.sqrt(1.0 + d.px[m]**2 + d.py[m]**2 + d.pz[m]**2)
        vec_m = d.c_dt / d.gamma[m]
        d.x[m] = d.x[m] + vec_m * d.px[m]
        d.y[m] = d.y[m] + vec_m * d.py[m]
        d.z[m] = d.z[m] + vec_m * d.pz[m]
    else:
        # Pas de masque, toutes les particules
        vec = d.boris_b_coef / d.gamma
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

        fx = d.boris_e_coef * d.ex
        fy = d.boris_e_coef * d.ey
        fz = d.boris_e_coef * d.ez

        v4 = vec * (fx + 2.0 * d.px)
        v5 = vec * (fy + 2.0 * d.py)
        v6 = vec * (fz + 2.0 * d.pz)
        d.px = d.px - (v2 + v3) * v4 + s12 * v5 + s13 * v6 + fx
        d.py = d.py + s21 * v4 - (v1 + v3) * v5 + s23 * v6 + fy
        d.pz = d.pz + s31 * v4 + s32 * v5 - (v1 + v2) * v6 + fz

        d.gamma = np.sqrt(1.0 + d.px**2 + d.py**2 + d.pz**2)
        vec = d.c_dt / d.gamma
        d.x = d.x + vec * d.px
        d.y = d.y + vec * d.py
        d.z = d.z + vec * d.pz


# =============================================================================
# calcbe - Calcul complet du champ B et E
# =============================================================================
def calcbe(d):
    """Calcul des champs magnétique et électrique pour toutes les particules."""
    hrhtetm1 = 1.0 / (d.map_dr * d.map_dtheta)
    hrm1 = 1.0 / d.map_dr
    htetm1 = 1.0 / d.map_dtheta
    hx1hy1m1 = 1.0 / (d.central_dx * d.central_dy)
    hx1m1 = 1.0 / d.central_dx
    hx2hy2m1 = 1.0 / (d.outer_dx * d.outer_dy)
    hx2m1 = 1.0 / d.outer_dx
    hy1m1 = 1.0 / d.central_dy
    hy2m1 = 1.0 / d.outer_dy
    ivec = np.zeros(d.n_particles, dtype=np.int32)
    jvec = np.zeros(d.n_particles, dtype=np.int32)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    veca = np.zeros(d.n_particles)
    vecb1 = np.zeros(d.n_particles)
    vecb2 = np.zeros(d.n_particles)
    vecb3 = np.zeros(d.n_particles)
    vecb4 = np.zeros(d.n_particles)
    vecb = np.zeros(d.n_particles)
    vecd1 = np.zeros(d.n_particles)
    vecd2 = np.zeros(d.n_particles)
    vecp1 = np.zeros(d.n_particles)
    vecp2 = np.zeros(d.n_particles)
    vecp3 = np.zeros(d.n_particles)
    vecp4 = np.zeros(d.n_particles)
    vecpc = np.zeros(d.n_particles)
    vecx = np.zeros(d.n_particles)
    vecy = np.zeros(d.n_particles)
    vecz = np.zeros(d.n_particles)
    vex1 = np.zeros(d.n_particles)
    vex2 = np.zeros(d.n_particles)
    vex3 = np.zeros(d.n_particles)
    vex4 = np.zeros(d.n_particles)
    vey1 = np.zeros(d.n_particles)
    vey2 = np.zeros(d.n_particles)
    vey3 = np.zeros(d.n_particles)
    vey4 = np.zeros(d.n_particles)
    np_ = d.n_particles

    # Calcul de B
    nbiz = np.count_nonzero((d.x / d.r > 1.0) | (d.x / d.r < -1.0))
    if nbiz != 0:
        print(f'nbiz : {nbiz}')
        for i in range(np_):
            rap = d.x[i] / d.r[i]
            if rap > 1.0 or rap < -1.0:
                print(f'{i} {d.x[i]} {d.y[i]} {d.r[i]} {rap}')

    veca = d.x / d.r
    mask_normal = (veca > -1.0) & (veca < 1.0)
    mask_ge1 = veca >= 1.0
    mask_lem1 = veca <= -1.0
    d.theta[mask_normal] = np.arccos(veca[mask_normal])
    d.theta[mask_ge1] = np.arccos(1.0)
    d.theta[mask_lem1] = np.arccos(-1.0)
    mask_yneg = d.y < 0.0
    d.theta[mask_yneg] = 2.0 * np.pi - d.theta[mask_yneg]

    ivec[:] = 0
    jvec[:] = 0
    mask = d.in_box
    ivec[mask] = ((d.r[mask] - d.map_rmin) * hrm1).astype(np.int32)
    jvec[mask] = (d.theta[mask] * htetm1).astype(np.int32)
    veca[:] = d.r - d.map_rmin - ivec * d.map_dr
    vecb[:] = d.theta - jvec * d.map_dtheta
    vecd1[:] = veca - d.map_dr
    vecd2[:] = vecb - d.map_dtheta
    vecz[:] = d.z * d.z * 0.5

    for k in range(np_):
        ik = ivec[k]
        jk = jvec[k]
        d.bz[k] = (vecd1[k] * vecd2[k] * (d.bz_map[ik, jk] -
                   vecz[k] * d.d2bz_map[ik, jk]) -
                   veca[k] * vecd2[k] * (d.bz_map[ik+1, jk] -
                   vecz[k] * d.d2bz_map[ik+1, jk]) +
                   veca[k] * vecb[k] * (d.bz_map[ik+1, jk+1] -
                   vecz[k] * d.d2bz_map[ik+1, jk+1]) -
                   vecd1[k] * vecb[k] * (d.bz_map[ik, jk+1] -
                   vecz[k] * d.d2bz_map[ik, jk+1])) * hrhtetm1
        d.bx[k] = (vecd1[k] * vecd2[k] * d.dbzdx_map[ik, jk] -
                   veca[k] * vecd2[k] * d.dbzdx_map[ik+1, jk] +
                   veca[k] * vecb[k] * d.dbzdx_map[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.dbzdx_map[ik, jk+1]
                   ) * d.z[k] * hrhtetm1
        d.by[k] = (vecd1[k] * vecd2[k] * d.dbzdy_map[ik, jk] -
                   veca[k] * vecd2[k] * d.dbzdy_map[ik+1, jk] +
                   veca[k] * vecb[k] * d.dbzdy_map[ik+1, jk+1] -
                   vecd1[k] * vecb[k] * d.dbzdy_map[ik, jk+1]
                   ) * d.z[k] * hrhtetm1

    mask_bz0 = (d.in_box & (d.bz == 0.0)) | (~d.in_box)
    d.bx[mask_bz0] = 0.0
    d.by[mask_bz0] = 0.0
    d.bz[mask_bz0] = -0.00001

    mask_rle = d.r <= d.map_rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.bz_map[0, 0]

    # Calcul de E au centre
    d.ex[:] = 0.0
    d.ey[:] = 0.0
    d.ez[:] = 0.0
    box1(d)

    if d.n_central > 0:
        m_c = d.in_central
        ivec[m_c] = ((d.x[m_c] - d.central_xmin) * hx1m1).astype(np.int32)
        jvec[m_c] = ((d.y[m_c] - d.central_ymin) * hy1m1).astype(np.int32)
        veca[m_c] = d.x[m_c] - d.central_xmin - ivec[m_c] * d.central_dx
        vecb[m_c] = d.y[m_c] - d.central_ymin - jvec[m_c] * d.central_dy
        vecd1[m_c] = veca[m_c] - d.central_dx
        vecd2[m_c] = vecb[m_c] - d.central_dy
        vecp1[m_c] = vecd1[m_c] * vecd2[m_c]
        vecp2[m_c] = -veca[m_c] * vecd2[m_c]
        vecp3[m_c] = veca[m_c] * vecb[m_c]
        vecp4[m_c] = -vecd1[m_c] * vecb[m_c]
        vecz[m_c] = d.z[m_c]**2 * 0.5

        for k in range(np_):
            if d.in_central[k]:
                ik = ivec[k]
                jk = jvec[k]
                vecb1[k] = d.central_d2pot[ik, jk]
                vecb2[k] = d.central_d2pot[ik+1, jk]
                vecb3[k] = d.central_d2pot[ik+1, jk+1]
                vecb4[k] = d.central_d2pot[ik, jk+1]
                vz = vecz[k]
                vex1[k] = (-(d.central_pot[ik+1, jk] - d.central_pot[ik-1, jk]) +
                             vz * (d.central_d2pot[ik+1, jk] - d.central_d2pot[ik-1, jk])
                             ) * 0.5 / d.central_dx
                vex2[k] = (-(d.central_pot[ik+2, jk] - d.central_pot[ik, jk]) +
                             vz * (d.central_d2pot[ik+2, jk] - d.central_d2pot[ik, jk])
                             ) * 0.5 / d.central_dx
                vex3[k] = (-(d.central_pot[ik+2, jk+1] - d.central_pot[ik, jk+1]) +
                             vz * (d.central_d2pot[ik+2, jk+1] - d.central_d2pot[ik, jk+1])
                             ) * 0.5 / d.central_dx
                vex4[k] = (-(d.central_pot[ik+1, jk+1] - d.central_pot[ik-1, jk+1]) +
                             vz * (d.central_d2pot[ik+1, jk+1] - d.central_d2pot[ik-1, jk+1])
                             ) * 0.5 / d.central_dx
                vey1[k] = (-(d.central_pot[ik, jk+1] - d.central_pot[ik, jk-1]) +
                             vz * (d.central_d2pot[ik, jk+1] - d.central_d2pot[ik, jk-1])
                             ) * 0.5 / d.central_dy
                vey2[k] = (-(d.central_pot[ik+1, jk+1] - d.central_pot[ik+1, jk-1]) +
                             vz * (d.central_d2pot[ik+1, jk+1] - d.central_d2pot[ik+1, jk-1])
                             ) * 0.5 / d.central_dy
                vey3[k] = (-(d.central_pot[ik+1, jk+2] - d.central_pot[ik+1, jk]) +
                             vz * (d.central_d2pot[ik+1, jk+2] - d.central_d2pot[ik+1, jk])
                             ) * 0.5 / d.central_dy
                vey4[k] = (-(d.central_pot[ik, jk+2] - d.central_pot[ik, jk]) +
                             vz * (d.central_d2pot[ik, jk+2] - d.central_d2pot[ik, jk])
                             ) * 0.5 / d.central_dy

        nto = int(d.omega_rf * d.t / (2.0 * np.pi))
        omhft = d.omega_rf * d.t - nto * 2.0 * np.pi
        pcright = d.v_dee_right * np.cos(omhft + d.phase_right)
        pcleft = d.v_dee_left * np.cos(omhft + d.phase_left)
        vecpc[d.x >= 0.0] = pcright
        vecpc[d.x < 0.0] = pcleft

        m_c = d.in_central
        d.ex[m_c] = vecpc[m_c] * (
            vex1[m_c] * vecp1[m_c] + vex2[m_c] * vecp2[m_c] +
            vex3[m_c] * vecp3[m_c] + vex4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.ey[m_c] = vecpc[m_c] * (
            vey1[m_c] * vecp1[m_c] + vey2[m_c] * vecp2[m_c] +
            vey3[m_c] * vecp3[m_c] + vey4[m_c] * vecp4[m_c]
        ) * hx1hy1m1
        d.ez[m_c] = d.z[m_c] * vecpc[m_c] * (
            vecp1[m_c] * vecb1[m_c] + vecp2[m_c] * vecb2[m_c] +
            vecp3[m_c] * vecb3[m_c] + vecp4[m_c] * vecb4[m_c]
        ) * hx1hy1m1

    # Calcul de E périphérique
    if (d.n_in_box == d.n_central) and (d.n_dipoles == 0) and (d.n_deflectors == 0) and (d.n_grids == 1):
        pass  # goto 300
    elif d.n_in_box == d.n_central:
        _calcbe_vertical(d)
        _calcbe_extraction(d)
        return
    else:
        dc1 = np.cos(d.dee_half_angle)
        ds1 = np.sin(d.dee_half_angle)
        mask_nc = ~d.in_central
        vecx[mask_nc] = np.abs(d.x[mask_nc]) * dc1 + (np.abs(d.y[mask_nc]) - d.gap_offset) * ds1
        vecy[mask_nc] = -np.abs(d.x[mask_nc]) * ds1 + (np.abs(d.y[mask_nc]) - d.gap_offset) * dc1

        box2(d, vecx, vecy)
        nto = int(d.omega_rf * d.t / (2.0 * np.pi))
        omhft = d.omega_rf * d.t - nto * 2.0 * np.pi

        if d.n_outer > 0:
            if d.n_central == 0:
                pcright = d.v_dee_right * np.cos(omhft + d.phase_right)
                pcleft = d.v_dee_left * np.cos(omhft + d.phase_left)
                vecpc[d.x >= 0.0] = pcright
                vecpc[d.x < 0.0] = pcleft

            m_p = d.in_outer
            ivec[m_p] = ((vecx[m_p] - d.outer_xmin) * hx2m1).astype(np.int32)
            jvec[m_p] = ((vecy[m_p] - d.outer_ymin) * hy2m1).astype(np.int32)
            veca[m_p] = vecx[m_p] - d.outer_xmin - ivec[m_p] * d.outer_dx
            vecb[m_p] = vecy[m_p] - d.outer_ymin - jvec[m_p] * d.outer_dy
            vecd1[m_p] = veca[m_p] - d.outer_dx
            vecd2[m_p] = vecb[m_p] - d.outer_dy
            vecp1[m_p] = vecd1[m_p] * vecd2[m_p]
            vecp2[m_p] = -veca[m_p] * vecd2[m_p]
            vecp3[m_p] = veca[m_p] * vecb[m_p]
            vecp4[m_p] = -vecd1[m_p] * vecb[m_p]
            vecz[m_p] = d.z[m_p]**2 * 0.5

            for k in range(np_):
                if d.in_outer[k]:
                    ik = ivec[k]
                    jk = jvec[k]
                    vec1[k] = d.outer_pot[ik, jk]
                    vec2[k] = d.outer_pot[ik+1, jk]
                    vec3[k] = d.outer_pot[ik+1, jk+1]
                    vec4[k] = d.outer_pot[ik, jk+1]
                    vecb1[k] = d.outer_d2pot[ik, jk]
                    vecb2[k] = d.outer_d2pot[ik+1, jk]
                    vecb3[k] = d.outer_d2pot[ik+1, jk+1]
                    vecb4[k] = d.outer_d2pot[ik, jk+1]

            d.ex[m_p] = -vecpc[m_p] * (
                vecd2[m_p] * (vec1[m_p] - vec2[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb2[m_p])) +
                vecb[m_p] * (vec3[m_p] - vec4[m_p] +
                                   vecz[m_p] * (vecb3[m_p] - vecb4[m_p]))
            ) * hx2hy2m1
            d.ey[m_p] = -vecpc[m_p] * (
                vecd1[m_p] * (vec1[m_p] - vec4[m_p] -
                                vecz[m_p] * (vecb1[m_p] - vecb4[m_p])) -
                veca[m_p] * (vec2[m_p] - vec3[m_p] +
                               vecz[m_p] * (vecb2[m_p] - vecb3[m_p]))
            ) * hx2hy2m1

            # Rotation inverse
            vec_tmp = d.ex[m_p] * dc1 - d.ey[m_p] * ds1
            d.ey[m_p] = d.ex[m_p] * ds1 + d.ey[m_p] * dc1
            d.ex[m_p] = vec_tmp

            d.ez[m_p] = d.z[m_p] * vecpc[m_p] * (
                vecp1[m_p] * vecb1[m_p] + vecp2[m_p] * vecb2[m_p] +
                vecp3[m_p] * vecb3[m_p] + vecp4[m_p] * vecb4[m_p]
            ) * hx2hy2m1

            # Symétries
            m1 = m_p & (d.x < 0.0) & (d.y > 0.0)
            d.ex[m1] = d.harmonic_sign * d.ex[m1]
            d.ey[m1] = -d.harmonic_sign * d.ey[m1]

            m2 = m_p & (d.x > 0.0) & (d.y < 0.0)
            d.ey[m2] = -d.ey[m2]

            m3 = m_p & (d.x < 0.0) & (d.y < 0.0)
            d.ex[m3] = d.harmonic_sign * d.ex[m3]
            d.ey[m3] = d.harmonic_sign * d.ey[m3]

            m4 = m_p & (d.x < 0.0)
            d.ez[m4] = -d.harmonic_sign * d.ez[m4]

    # Tri vertical
    _calcbe_vertical(d)

    # Extraction
    _calcbe_extraction(d)


def _calcbe_vertical(d):
    """Ajout du test tri vertical."""
    d.sorter_voltage_seen = 0.0
    if d.sorter_gap == 0.0:
        return
    mask_ver = ((d.r >= d.sorter_rmin) & (d.r <= d.sorter_rmax) &
                (d.theta >= d.sorter_theta1) & (d.theta <= d.sorter_theta2))
    d.ez[mask_ver] = (d.ez[mask_ver] -
        d.sorter_voltage * (1.0 - d.sorter_dv) / d.sorter_gap / 2.0 *
        np.cos(2.0 * d.omega_rf * d.t * (1.0 + d.sorter_nuz / 2.0 / d.harmonic) +
               2.0 * d.sorter_phase - np.pi / 2.0) -
        d.sorter_voltage * (1.0 + d.sorter_dv) / d.sorter_gap / 2.0 *
        np.cos(2.0 * d.omega_rf * d.t * (1.0 - d.sorter_nuz / 2.0 / d.harmonic) +
               2.0 * d.sorter_phase - np.pi / 2.0))
    if (d.r[0] >= d.sorter_rmin and d.r[0] <= d.sorter_rmax and
            d.theta[0] >= d.sorter_theta1 and d.theta[0] <= d.sorter_theta2):
        d.sorter_voltage_seen = -d.ez[0] * d.sorter_gap


def _calcbe_extraction(d):
    """Calcul de l'extraction."""
    for i in range(d.n_deflectors):
        boxdf(d, i)
        if d.n_in_element != 0:
            mask = d.in_element
            vec_p = np.sqrt(d.px[mask]**2 + d.py[mask]**2)
            d.ex[mask] = d.deflector_e[i] * d.py[mask] / vec_p
            d.ey[mask] = -d.deflector_e[i] * d.px[mask] / vec_p

    for i in range(d.n_dipoles):
        boxdm(d, i)
        if d.n_in_element != 0:
            mask = d.in_element
            vec_dm = np.sqrt((d.x[mask] - d.dipole_x[i])**2 + (d.y[mask] - d.dipole_y[i])**2)
            vec1_dm = vec_dm - d.dipole_r[i]
            d.bx[mask] += d.z[mask] * d.dipole_gradient[i] * (d.x[mask] - d.dipole_x[i]) / vec_dm
            d.by[mask] += d.z[mask] * d.dipole_gradient[i] * (d.y[mask] - d.dipole_y[i]) / vec_dm
            d.bz[mask] += vec1_dm * d.dipole_gradient[i]

    # Grilles d'extraction
    if d.extraction < 0.0 and d.n_grids > 0:
        if d.step >= d.grid_start_step:
            for kg in range(d.n_grids):
                boxgex(d, kg)
                d.grid_reached[d.crossed_grid] = kg + 1  # 1-based in data
                for i in range(d.n_particles):
                    if d.crossed_grid[i]:
                        if (70 + kg + 1) in files_dict:
                            files_dict[70 + kg + 1].write(
                                f' {i+1:5d} {d.step:6d}'
                                f' {d.x[i]:15.8e} {d.y[i]:15.8e} {d.z[i]:15.8e}'
                                f' {d.px[i]:15.8e} {d.py[i]:15.8e} {d.pz[i]:15.8e}'
                                f' {d.bz[i]:8.5f}\n')


# =============================================================================
# calcb - Calcul de B simplifié (utilisé pendant init)
# =============================================================================
def calcb(d):
    """Calcul simplifié du champ magnétique (pendant l'initialisation)."""
    hrhtet = d.map_dr * d.map_dtheta
    ivec = np.zeros(d.n_particles, dtype=np.int32)
    jvec = np.zeros(d.n_particles, dtype=np.int32)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    veca = np.zeros(d.n_particles)
    vecb1 = np.zeros(d.n_particles)
    vecb2 = np.zeros(d.n_particles)
    vecb3 = np.zeros(d.n_particles)
    vecb4 = np.zeros(d.n_particles)
    vecb = np.zeros(d.n_particles)
    vecd1 = np.zeros(d.n_particles)
    vecd2 = np.zeros(d.n_particles)
    vecp1 = np.zeros(d.n_particles)
    vecp2 = np.zeros(d.n_particles)
    vecp3 = np.zeros(d.n_particles)
    vecp4 = np.zeros(d.n_particles)
    vecx1 = np.zeros(d.n_particles)
    vecx2 = np.zeros(d.n_particles)
    vecx3 = np.zeros(d.n_particles)
    vecx4 = np.zeros(d.n_particles)
    vecy1 = np.zeros(d.n_particles)
    vecy2 = np.zeros(d.n_particles)
    vecy3 = np.zeros(d.n_particles)
    vecy4 = np.zeros(d.n_particles)
    vecz = np.zeros(d.n_particles)
    d.theta = np.arccos(np.clip(d.x / d.r, -1.0, 1.0))
    d.theta[d.y < 0.0] = 2.0 * np.pi - d.theta[d.y < 0.0]

    mask = d.in_box
    ivec[mask] = ((d.r[mask] - d.map_rmin) / d.map_dr).astype(np.int32)
    jvec[mask] = (d.theta[mask] / d.map_dtheta).astype(np.int32)
    veca[mask] = d.r[mask] - d.map_rmin - ivec[mask] * d.map_dr
    vecb[mask] = d.theta[mask] - jvec[mask] * d.map_dtheta
    vecd1[mask] = veca[mask] - d.map_dr
    vecd2[mask] = vecb[mask] - d.map_dtheta
    vecp1[mask] = vecd1[mask] * vecd2[mask]
    vecp2[mask] = -veca[mask] * vecd2[mask]
    vecp3[mask] = veca[mask] * vecb[mask]
    vecp4[mask] = -vecd1[mask] * vecb[mask]
    vecz[mask] = d.z[mask]**2 / 2.0

    for k in range(d.n_particles):
        if d.in_box[k]:
            ik = ivec[k]
            jk = jvec[k]
            vk = vecz[k]
            vecb1[k] = d.bz_map[ik, jk]
            vecb2[k] = d.bz_map[ik+1, jk]
            vecb3[k] = d.bz_map[ik+1, jk+1]
            vecb4[k] = d.bz_map[ik, jk+1]
            vec1[k] = -vk * d.d2bz_map[ik, jk]
            vec2[k] = -vk * d.d2bz_map[ik+1, jk]
            vec3[k] = -vk * d.d2bz_map[ik+1, jk+1]
            vec4[k] = -vk * d.d2bz_map[ik, jk+1]
            vecx1[k] = d.dbzdx_map[ik, jk]
            vecx2[k] = d.dbzdx_map[ik+1, jk]
            vecx3[k] = d.dbzdx_map[ik+1, jk+1]
            vecx4[k] = d.dbzdx_map[ik, jk+1]
            vecy1[k] = d.dbzdy_map[ik, jk]
            vecy2[k] = d.dbzdy_map[ik+1, jk]
            vecy3[k] = d.dbzdy_map[ik+1, jk+1]
            vecy4[k] = d.dbzdy_map[ik, jk+1]

    vec1[mask] += vecb1[mask]
    vec2[mask] += vecb2[mask]
    vec3[mask] += vecb3[mask]
    vec4[mask] += vecb4[mask]
    d.bz[mask] = (vecp1[mask] * vec1[mask] + vecp2[mask] * vec2[mask] +
                  vecp3[mask] * vec3[mask] + vecp4[mask] * vec4[mask]) / hrhtet
    d.bx[mask] = d.z[mask] * (
        vecp1[mask] * vecx1[mask] + vecp2[mask] * vecx2[mask] +
        vecp3[mask] * vecx3[mask] + vecp4[mask] * vecx4[mask]) / hrhtet
    d.by[mask] = d.z[mask] * (
        vecp1[mask] * vecy1[mask] + vecp2[mask] * vecy2[mask] +
        vecp3[mask] * vecy3[mask] + vecp4[mask] * vecy4[mask]) / hrhtet

    mask_rle = d.r <= d.map_rmin
    d.bx[mask_rle] = 0.0
    d.by[mask_rle] = 0.0
    d.bz[mask_rle] = d.bz_map[0, 0]


# =============================================================================
# caldxs - Calcul de l'espace de phase
# =============================================================================
def caldxs(d):
    """Calcul des coordonnées dans l'espace de phase."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    for k in range(d.n_families):
        nbr = d.family_size[k]
        i1 = d.family_start[k]
        i2 = i1 + nbr
        d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
        rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]
        d.ps_cartesian[0, :nbr] = d.x[i1:i2] - d.x[i1]
        d.ps_cartesian[4, :nbr] = d.y[i1:i2] - d.y[i1]
        ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
        tetai1 = np.arccos(d.x[i1] / ri1)
        if d.y[i1] < 0.0:
            tetai1 = 2.0 * np.pi - tetai1
        pri1 = (d.x[i1] * d.px[i1] + d.y[i1] * d.py[i1]) / ri1
        tpri1 = np.arcsin(pri1 / d.p[i1])
        co = np.cos(tetai1 - tpri1)
        si = np.sin(tetai1 - tpri1)
        vec[:nbr] = d.ps_cartesian[0, :nbr] * co + d.ps_cartesian[4, :nbr] * si
        d.ps_cartesian[4, :nbr] = -d.ps_cartesian[0, :nbr] * si + d.ps_cartesian[4, :nbr] * co
        d.ps_cartesian[0, :nbr] = vec[:nbr]
        for ii in range(nbr):
            d.ps[4, ii] = rc * np.arctan(d.ps_cartesian[4, ii] / (rc + d.ps_cartesian[0, ii]))
        vec[:nbr] = np.sin(d.ps[4, :nbr] / 2.0 / rc)**2
        d.ps[0, :nbr] = (d.ps_cartesian[0, :nbr] + 2.0 * rc * vec[:nbr]) / (1.0 - 2.0 * vec[:nbr])
        vec1[:nbr] = np.cos(d.ps[4, :nbr] / rc)
        vec2[:nbr] = np.sin(d.ps[4, :nbr] / rc)
        d.ps_cartesian[1, :nbr] = d.px[i1:i2] - d.px[i1] * vec1[:nbr] + d.py[i1] * vec2[:nbr]
        d.ps_cartesian[5, :nbr] = d.py[i1:i2] - d.px[i1] * vec2[:nbr] - d.py[i1] * vec1[:nbr]
        vec[:nbr] = d.ps_cartesian[1, :nbr] * co + d.ps_cartesian[5, :nbr] * si
        d.ps_cartesian[5, :nbr] = -d.ps_cartesian[1, :nbr] * si + d.ps_cartesian[5, :nbr] * co
        d.ps_cartesian[1, :nbr] = vec[:nbr]
        d.ps[1, :nbr] = d.ps_cartesian[1, :nbr] * vec1[:nbr] + d.ps_cartesian[5, :nbr] * vec2[:nbr]
        d.ps[5, :nbr] = -d.ps_cartesian[1, :nbr] * vec2[:nbr] + d.ps_cartesian[5, :nbr] * vec1[:nbr]
        d.ps[2, :nbr] = d.z[i1:i2] - d.z[i1]
        d.ps[3, :nbr] = d.pz[i1:i2] - d.pz[i1]
        d.ps[5, :nbr] = np.sqrt(
            (d.p[i1] + d.ps[5, :nbr])**2 +
            d.ps[1, :nbr]**2 + d.ps[3, :nbr]**2) - d.p[i1]

        for i in range(6):
            for j in range(6):
                d.sigma[i, j] = np.sum(d.ps[i, :nbr] * d.ps[j, :nbr]) / nbr

        tetai1 = tetai1 * 180.0 / np.pi
        ri1mm = ri1 * 1000.0
        d.dc_transform[:] = 0.0
        d.dc_transform[0, 0] = -1.0
        d.dc_transform[1, 1] = -1.0
        for i in range(2, 6):
            d.dc_transform[i, i] = 1.0
        d.dc_transform[0, 5] = ri1 / d.p[i1]
        d.dc_transform[4, 1] = -d.dc_transform[0, 5]
        d.sigma_dc = d.dc_transform @ d.sigma @ d.dc_transform.T


# =============================================================================
# calgri / calgric - Calcul sur grille
# =============================================================================
def calgri(d):
    """Calcul des coordonnées sur grille."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    vec5 = np.zeros(d.n_particles)
    vec6 = np.zeros(d.n_particles)
    for k in range(d.n_families):
        nbr = d.family_size[k]
        i1 = d.family_start[k]
        i2 = i1 + nbr
        d.p[i1:i2] = np.sqrt(d.px[i1:i2]**2 + d.py[i1:i2]**2 + d.pz[i1:i2]**2)
        vec[:nbr] = -d.mass_c[k] * d.p[i1:i2] / d.charge[i1:i2] / d.bz[i1]
        vec1[i1:i2] = d.x[i1:i2] + d.py[i1:i2] * d.dt / d.boris_e_coef[i1] / d.bz[i1]
        vec2[i1:i2] = d.y[i1:i2] - d.px[i1:i2] * d.dt / d.boris_e_coef[i1] / d.bz[i1]
        ux = (d.x[i1] - vec1[i1]) / vec[0]
        uy = (d.y[i1] - vec2[i1]) / vec[0]
        vec3[i1:i2] = ((d.x[i1] - vec1[i1:i2]) * ux +
                         (d.y[i1] - vec2[i1:i2]) * uy)
        vec4[i1:i2] = ((d.x[i1] - vec1[i1:i2])**2 +
                         (d.y[i1] - vec2[i1:i2])**2 - vec[:nbr]**2)
        d.proj_offset[i1:i2] = -vec3[i1:i2] + np.sqrt(np.abs(vec3[i1:i2]**2 - vec4[i1:i2]))
        d.proj_x[i1:i2] = d.x[i1] + d.proj_offset[i1:i2] * ux
        d.proj_y[i1:i2] = d.y[i1] + d.proj_offset[i1:i2] * uy
        vec3[i1:i2] = np.sqrt((d.proj_x[i1:i2] - vec1[i1:i2])**2 +
                                (d.proj_y[i1:i2] - vec2[i1:i2])**2)
        vec4[i1:i2] = np.sqrt((d.proj_x[i1:i2] - vec1[i1])**2 +
                                (d.proj_y[i1:i2] - vec2[i1])**2)
        vec5[i1:i2] = (
            (d.proj_x[i1:i2] - vec1[i1:i2]) * (d.proj_y[i1:i2] - vec2[i1]) -
            (d.proj_y[i1:i2] - vec2[i1:i2]) * (d.proj_x[i1:i2] - vec1[i1])
        ) / vec3[i1:i2] / vec4[i1:i2]
        d.proj_angle[i1:i2] = np.arcsin(vec5[i1:i2])
        vec4[i1:i2] = np.sqrt((d.x[i1:i2] - vec1[i1:i2])**2 +
                                (d.y[i1:i2] - vec2[i1:i2])**2)
        vec6[i1:i2] = (
            (d.proj_x[i1:i2] - vec1[i1:i2]) * (d.y[i1:i2] - vec2[i1:i2]) -
            (d.proj_y[i1:i2] - vec2[i1:i2]) * (d.x[i1:i2] - vec1[i1:i2])
        ) / vec3[i1:i2] / vec4[i1:i2]
        d.proj_phase[i1:i2] = d.harmonic * np.arcsin(vec6[i1:i2])
    return vec


def calgric(d):
    """Calcul sur grille pour corrélations."""
    nbr = d.family_size[0]
    vec = calgri(d)
    d.ps_projected[0, :nbr] = d.proj_offset[:nbr]
    d.ps_projected[1, :nbr] = d.proj_angle[:nbr] * d.p[0]
    d.ps_projected[2, :nbr] = d.ps[2, :nbr]
    d.ps_projected[3, :nbr] = d.ps[3, :nbr]
    d.ps_projected[4, :nbr] = d.proj_phase[:nbr] * vec[:nbr] / d.harmonic
    d.ps_projected[5, :nbr] = d.ps[5, :nbr]


# =============================================================================
# Variable globale pour les fichiers ouverts
# =============================================================================
files_dict = {}


def distribute_bunch(d, sigma_matrices):
    """Bunch from the sigma matrices (after the field map is read), then optional bunch load."""
    dxtmp = np.zeros((6, d.n_particles))
    rot = np.zeros((6, 6))
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    if d.correlations != 0:
        m = 6
        boite(d)
        calcb(d)
        for k in range(d.n_families):
            nbr = d.family_size[k]
            if nbr <= 1:
                continue  # une seule particule, pas de distribution a appliquer
            ps = remp(m, nbr)
            d.ps[:, :nbr] = ps
            rcor = choles(m, sigma_matrices[k])
            trcor = rcor.T
            dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    dxtmp[i, :nbr] += trcor[i, j] * d.ps[j, :nbr]
            d.ps[:6, :nbr] = dxtmp[:6, :nbr]

            i1 = d.family_start[k]
            i2 = i1 + d.family_size[k]
            rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]

            if d.r0[k] > d.central_rmax:
                print('prise en compte de sqrt(cos(phi))...')
                vec4[:nbr] = np.sqrt(np.cos(d.harmonic * d.ps[4, :nbr] / rc))
                d.ps[0, :nbr] = d.ps[0, :nbr] * vec4[:nbr] - rc * (1.0 - vec4[:nbr])
            else:
                print('pas de sqrt(cos(phi))...')

            vec1[:nbr] = np.cos(d.ps[4, :nbr] / rc)
            vec2[:nbr] = np.sin(d.ps[4, :nbr] / rc)
            vec3[:nbr] = np.sin(d.ps[4, :nbr] / rc / 2.0)**2
            d.ps_cartesian[0, :nbr] = d.ps[0, :nbr] - 2.0 * (rc + d.ps[0, :nbr]) * vec3[:nbr]
            d.ps[5, :nbr] = np.sqrt(
                (d.p[i1] + d.ps[5, :nbr])**2 -
                d.ps[1, :nbr]**2 - d.ps[3, :nbr]**2) - d.p[i1]
            d.ps_cartesian[4, :nbr] = (rc + d.ps[0, :nbr]) * vec2[:nbr]
            d.ps_cartesian[1, :nbr] = d.ps[1, :nbr] * vec1[:nbr] - d.ps[5, :nbr] * vec2[:nbr]
            d.ps_cartesian[5, :nbr] = d.ps[1, :nbr] * vec2[:nbr] + d.ps[5, :nbr] * vec1[:nbr]
            d.ps_cartesian[2:4, :nbr] = d.ps[2:4, :nbr]

            co = np.cos(d.theta0[k] - d.pr0[k])
            si = np.sin(d.theta0[k] - d.pr0[k])
            rot[:] = 0.0
            rot[0, 0] = co
            rot[0, 4] = -si
            rot[1, 1] = co
            rot[1, 5] = -si
            rot[2, 2] = 1.0
            rot[3, 3] = 1.0
            rot[4, 0] = si
            rot[4, 4] = co
            rot[5, 1] = si
            rot[5, 5] = co
            dxtmp[:6, :nbr] = 0.0
            for i in range(6):
                for j in range(6):
                    dxtmp[i, :nbr] += rot[i, j] * d.ps_cartesian[j, :nbr]
            d.ps_cartesian[:6, :nbr] = dxtmp[:6, :nbr]

            d.x[i1:i2] += d.ps_cartesian[0, :nbr]
            d.y[i1:i2] += d.ps_cartesian[4, :nbr]
            vec[:nbr] = d.px[i1:i2] * vec1[:nbr] - d.py[i1:i2] * vec2[:nbr]
            d.py[i1:i2] = d.px[i1:i2] * vec2[:nbr] + d.py[i1:i2] * vec1[:nbr] + d.ps_cartesian[5, :nbr]
            d.px[i1:i2] = vec[:nbr] + d.ps_cartesian[1, :nbr]
            if d.phase_space[:2] in ('3d', '6d'):
                d.z[i1:i2] += d.ps_cartesian[2, :nbr]
                d.pz[i1:i2] += d.ps_cartesian[3, :nbr]
            if d.r0[k] > d.central_rmax:
                d.px[i1:i2] *= vec4[:nbr]
                d.py[i1:i2] *= vec4[:nbr]
                d.pz[i1:i2] *= vec4[:nbr]

        if d.bunch_load_file != 'no':
            print('attention : chargement paquet d ions...')
            with open(d.bunch_load_file, 'rb') as f26:
                for i in range(d.n_particles):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.x[i], d.y[i], d.z[i] = data
                for i in range(d.n_particles):
                    data = np.fromfile(f26, dtype=np.float64, count=3)
                    if len(data) == 3:
                        d.px[i], d.py[i], d.pz[i] = data


# =============================================================================
# init - Initialisation du programme
# =============================================================================
def init(d, ficin):
    """Initialisation : lecture du fichier de données ficin, cartes de champ, etc."""
    dbzxlu = np.zeros((250, 91))
    dbzylu = np.zeros((250, 91))
    files_dict.clear()


    d.deck_file = ficin
    with open(d.deck_file, 'r') as f20_raw:
      f20 = FortranFile(f20_raw)
      with f20:  # dummy context for indentation
        d.title = f20.readline().strip().strip('"').strip("'")
        vals = f20.readline().split()
        d.steps_per_turn = int(vals[0])
        vals = f20.readline().split()
        d.design_r_final = fortran_float(vals[0])
        d.design_turns = int(vals[1])
        vals = f20.readline().split()
        d.rf_frequency = fortran_float(vals[0])
        d.dee_half_angle = fortran_float(vals[2])
        d.harmonic = fortran_float(vals[3])
        d.db_over_b = fortran_float(vals[4])
        vals = f20.readline().split()
        d.phase_right = fortran_float(vals[0]) * np.pi / 180.0
        d.phase_left = fortran_float(vals[1]) * np.pi / 180.0
        vals = f20.readline().split()
        d.v_dee_right = -fortran_float(vals[0])
        d.v_dee_left = -fortran_float(vals[1])
        vals = f20.readline().split()
        d.harmonic_n_coef[0] = int(vals[0])
        if d.harmonic_n_coef[0] != 0:
            vals = f20.readline().split()
            for i in range(d.harmonic_n_coef[0]):
                d.harmonic_coef[i, 0] = fortran_float(vals[i])
            d.harmonic_phase[0] = fortran_float(f20.readline().split()[0])
        vals = f20.readline().split()
        d.harmonic_n_coef[1] = int(vals[0])
        if d.harmonic_n_coef[1] != 0:
            vals = f20.readline().split()
            for i in range(d.harmonic_n_coef[1]):
                d.harmonic_coef[i, 1] = fortran_float(vals[i])
            d.harmonic_phase[1] = fortran_float(f20.readline().split()[0])

        d.gap_offset = 0.0
        if d.dee_half_angle < 0.0:
            print('attention : gap en biais...')
            print('on remet teta1 positif')
            print('et on va lire le decalage...')
            d.dee_half_angle = -d.dee_half_angle
            d.gap_offset = fortran_float(f20.readline().split()[0])
            print(f'decalage yomega : {d.gap_offset}')

        vals = f20.readline().split()
        d.zmin = fortran_float(vals[0])
        d.zmax = fortran_float(vals[1])
        d.central_rmax = fortran_float(vals[2])
        vals = f20.readline().split()
        d.extraction = fortran_float(vals[0])
        d.cut_bunch = fortran_float(vals[1])
        d.grids_on = fortran_float(vals[2])
        d.grid_start = fortran_float(vals[3])
        d.n_deflectors = 0
        d.n_dipoles = 0

        if d.extraction < 0.0:
            if d.cut_bunch <= 0.0:
                d.extraction_start_step = 0
                print('extraction : le paquet sera coupe')
            else:
                d.extraction_start_step = int(d.cut_bunch + 0.5)
                print('extraction : le paquet sera entier')
            d.n_deflectors = int(f20.readline().split()[0])
            if d.n_deflectors > 0:
                for i in range(d.n_deflectors):
                    vals = f20.readline().split()
                    d.deflector_r[i] = fortran_float(vals[0])
                    d.deflector_e[i] = fortran_float(vals[1])
                    d.deflector_in_x[i] = fortran_float(vals[2])
                    d.deflector_in_y[i] = fortran_float(vals[3])
                    d.deflector_out_x[i] = fortran_float(vals[4])
                    d.deflector_out_y[i] = fortran_float(vals[5])
                    tetdfid = fortran_float(vals[6])
                    tetdfod = fortran_float(vals[7])
                    d.deflector_in_angle[i] = tetdfid * np.pi / 180.0
                    d.deflector_out_angle[i] = tetdfod * np.pi / 180.0
                    d.deflector_in_ux[i] = np.cos(d.deflector_in_angle[i])
                    d.deflector_in_uy[i] = np.sin(d.deflector_in_angle[i])
                    d.deflector_out_ux[i] = np.cos(d.deflector_out_angle[i])
                    d.deflector_out_uy[i] = np.sin(d.deflector_out_angle[i])
            d.n_dipoles = int(f20.readline().split()[0])
            if d.n_dipoles > 0:
                for i in range(d.n_dipoles):
                    vals = f20.readline().split()
                    d.dipole_r[i] = fortran_float(vals[0])
                    d.dipole_gradient[i] = fortran_float(vals[1])
                    d.dipole_x[i] = fortran_float(vals[2])
                    d.dipole_y[i] = fortran_float(vals[3])
                    tetdmid = fortran_float(vals[4])
                    tetdmod = fortran_float(vals[5])
                    d.dipole_in_angle[i] = tetdmid * np.pi / 180.0
                    d.dipole_out_angle[i] = tetdmod * np.pi / 180.0
                    d.dipole_in_ux[i] = np.cos(d.dipole_in_angle[i])
                    d.dipole_in_uy[i] = np.sin(d.dipole_in_angle[i])
                    d.dipole_out_ux[i] = np.cos(d.dipole_out_angle[i])
                    d.dipole_out_uy[i] = np.sin(d.dipole_out_angle[i])
            if d.grids_on > 0.0:
                print('grilles extraction...')
                d.grid_start_step = int(d.grid_start + 0.5)
                d.n_grids = int(f20.readline().split()[0])
                for kg in range(d.n_grids):
                    vals = f20.readline().split()
                    d.grid_angle_deg[kg] = fortran_float(vals[0])
                    d.grid_x[kg] = fortran_float(vals[1]) / 1000.0
                    d.grid_y[kg] = fortran_float(vals[2]) / 1000.0
                    d.grid_xc[kg] = fortran_float(vals[3]) / 1000.0
                    d.grid_yc[kg] = fortran_float(vals[4]) / 1000.0
                    d.grid_length[kg] = np.sqrt((d.grid_x[kg] - d.grid_xc[kg])**2 +
                                          (d.grid_y[kg] - d.grid_yc[kg])**2)
                    d.grid_ux[kg] = (d.grid_x[kg] - d.grid_xc[kg]) / d.grid_length[kg]
                    d.grid_uy[kg] = (d.grid_y[kg] - d.grid_yc[kg]) / d.grid_length[kg]
                d.grid_file_prefix = f20.readline().strip()
                for kg in range(d.n_grids):
                    ficfingex = str(kg + 1)
                    d.grid_files[kg] = d.grid_file_prefix + '.' + ficfingex
                    fgex = open(d.grid_files[kg], 'w')
                    files_dict[70 + kg + 1] = fgex
                    fgex.write('teta x y xc yc : \n')
                    fgex.write(f'{d.grid_angle_deg[kg]} {d.grid_x[kg]} {d.grid_y[kg]}'
                               f' {d.grid_xc[kg]} {d.grid_yc[kg]}\n')

        d.dt = d.harmonic / d.rf_frequency / d.steps_per_turn
        d.omega_rf = 2.0 * np.pi * d.rf_frequency
        d.omega_rev = d.omega_rf / d.harmonic
        d.c_dt = C_LIGHT * d.dt
        d.dee_half_angle = d.dee_half_angle * np.pi / 180.0
        d.cut_bunch = d.cut_bunch * np.pi / 180.0
        d.grids_on = d.grids_on * np.pi / 180.0
        iharm = int(d.harmonic + 0.5)
        d.harmonic_sign = 1.0
        if (iharm // 2) * 2 == iharm:
            d.harmonic_sign = -1.0

        d.n_families = int(f20.readline().split()[0])
        if d.n_families > 1:
            d.n_particles = int(f20.readline().split()[0])
        d.family_start[0] = 0  # 0-based
        n_loaded = 0

        for k in range(d.n_families):
            vals = f20.readline().split()
            d.family_size[k] = int(vals[0])
            d.atomic_mass[k] = fortran_float(vals[1])
            d.mass_number[k] = fortran_float(vals[2])
            d.charge_state[k] = fortran_float(vals[3])
            if k == 0:
                if d.n_families == 1:
                    d.n_particles = d.family_size[0]
                d.allocate()
                d.extracted[:] = False
                d.grid_reached[:] = 0
                d.alive[:] = True
                d.n_alive = d.n_particles
                d.n_traps = 0
                d.trapped_follow = 1

            vals = f20.readline().split()
            d.r0[k] = fortran_float(vals[0])
            d.theta0[k] = fortran_float(vals[1])
            d.z0[k] = fortran_float(vals[2])
            d.pr0[k] = fortran_float(vals[3])
            d.pz0[k] = fortran_float(vals[4])
            w0mev = fortran_float(vals[5])
            d.phi0_deg[k] = fortran_float(vals[6])
            if k == 0:
                d.w0_mev_per_u = w0mev

            d.theta0[k] *= np.pi / 180.0
            d.pr0[k] *= 0.001
            d.pz0[k] *= 0.001
            n_loaded += d.family_size[k]
            if n_loaded > d.n_particles:
                print(f'n np : {d.n_particles} {n_loaded}')
                print('changer n ou np...')
                sys.exit(1)
            d.family_start[k + 1] = d.family_start[k] + d.family_size[k]
            i1 = d.family_start[k]
            i2 = d.family_start[k] + d.family_size[k]
            d.x[i1:i2] = d.r0[k] * np.cos(d.theta0[k])
            d.y[i1:i2] = d.r0[k] * np.sin(d.theta0[k])
            d.z[i1:i2] = d.z0[k]
            d.w0[k] = w0mev * 1.60217733e-13 * d.mass_number[k]
            am0k = 1.6605402e-27 * (d.atomic_mass[k] - 0.0005486 * d.charge_state[k])
            d.mass[i1:i2] = am0k
            d.mass_c[i1:i2] = d.mass[i1:i2] * C_LIGHT
            d.rest_energy[i1:i2] = d.mass[i1:i2] * C_LIGHT**2
            d.gamma[i1:i2] = 1.0 + d.w0[k] / d.rest_energy[i1:i2]
            alpha = (np.pi / 2.0 + d.theta0[k]) - d.pr0[k]
            d.p[i1:i2] = np.sqrt((d.gamma[i1:i2] - 1.0) * (d.gamma[i1:i2] + 1.0))
            d.px[i1:i2] = d.p[i1:i2] * np.cos(alpha) * np.cos(d.pz0[k])
            d.py[i1:i2] = d.p[i1:i2] * np.sin(alpha) * np.cos(d.pz0[k])
            d.pz[i1:i2] = d.p[i1:i2] * np.sin(d.pz0[k])
            d.charge[i1:i2] = d.charge_state[k] * 1.60217733e-19
            d.boris_b_coef[i1:i2] = d.charge[i1:i2] * d.dt / 2.0 / d.mass[i1:i2]
            d.boris_e_coef[i1:i2] = d.charge[i1:i2] * d.dt / d.mass_c[i1:i2]

        if d.grids_on > 0.0 and d.n_grids > 0:
            for kg in range(d.n_grids):
                f = files_dict[70 + kg + 1]
                f.write('caracteristiques particules :\n')
                f.write(f'{d.n_particles} {d.mass_number[0]} {d.atomic_mass[0]} {d.charge_state[0]}\n')
                f.write('caracteristiques HF :\n')
                f.write(f'{d.steps_per_turn} {d.rf_frequency} {d.harmonic}\n')

        d.correlations = int(f20.readline().split()[0])
        d.sigma_file = "pas de fichier correl..."
        sigma_matrices = []
        if d.correlations != 0:
            m = 6
            d.sigma_file = f20.readline().strip()
            d.phase_space = f20.readline().strip()
            f25 = open(d.sigma_file, 'r')
            vals = f20.read_tokens(3)
            d.correlation_file = vals[0]
            d.matrix_file = vals[1]
            d.dc_matrix_file = vals[2]
            if d.correlation_file != 'no':
                files_dict[17] = open(d.correlation_file, 'w')
                files_dict[17].write(
                    '  tr     r      teta    dr    ds      dphi'
                    '   epsr    dwsw        wmpa     dw          dwdphi'
                    '       dz     epsz\n')
            if d.matrix_file != 'no':
                files_dict[18] = open(d.matrix_file, 'w')
            if d.dc_matrix_file != 'no':
                files_dict[19] = open(d.dc_matrix_file, 'w')
                files_dict[19].write(
                    '  tr     r      teta   drdc  dsdc  dphidc'
                    ' epsrdc  dwswdc        wmpa   dwdc        dwdphidc\n')

            # Pre-lire les matrices sigma (calcul differe apres chargement champ magnetique)
            for k in range(d.n_families):
                cor_mat = np.zeros((6, 6), dtype=np.float64)
                for i in range(6):
                    vals_sig = f25.readline().split()
                    for j in range(6):
                        cor_mat[i, j] = fortran_float(vals_sig[j])
                sigma_matrices.append(cor_mat)
            f25.close()

        if n_loaded != d.n_particles:
            print(f'n np : {d.n_particles} {n_loaded}')
            print('changer n ou np...')
            sys.exit(1)

        d.t0 = d.phi0_deg[0] * np.pi / 180.0 / d.omega_rf
        d.t = d.t0

        vals = f20.readline().split()
        d.n_steps = int(vals[0])
        d.output_every = int(vals[1])
        vals = f20.read_tokens(2)
        d.z_output_file = vals[0]
        d.xy_output_file = vals[1]
        vals = f20.read_tokens(3)
        d.bfield_file = vals[0]
        d.central_potential_file = vals[1]
        d.outer_potential_file = vals[2]
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
        d.bunch_load_file, d.bunch_save_file = (tail + ['no', 'no'])[:2]

    # Lecture carte magnétique
    with open(d.bfield_file, 'r') as f21_raw:
      f21 = FortranFile(f21_raw)
      with f21:
        title = f21.readline().strip()
        vals = f21.readline().split()
        d.map_rmin = fortran_float(vals[0])
        d.map_rmax = fortran_float(vals[1])
        d.map_nr = int(vals[2])
        d.map_ntheta = int(vals[3])

        if d.map_ntheta != 91:
            print('lecture carte tosca.big ou carte mesuree ...')
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.bz_map[i, j] = fortran_float(vals[j])
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.d2bz_map[i, j] = fortran_float(vals[j])
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.dbzdx_map[i, j] = fortran_float(vals[j])
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.dbzdy_map[i, j] = fortran_float(vals[j])
            d.map_dtheta = 2.0 * np.pi / (d.map_ntheta - 1.0)
            ntetall = d.map_ntheta
        else:
            print('carte tosca.small (1/8eme) ...')
            # Lecture avec inversion de l'index j
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.bz_map[i, 91 - j] = fortran_float(vals[j])
            for i in range(d.map_nr):
                vals = f21.read_tokens(d.map_ntheta)
                for j in range(d.map_ntheta):
                    d.d2bz_map[i, 91 - j] = fortran_float(vals[j])
            try:
                for i in range(d.map_nr):
                    vals = f21.read_tokens(d.map_ntheta)
                    for j in range(d.map_ntheta):
                        dbzxlu[i, j] = fortran_float(vals[j])
                for i in range(d.map_nr):
                    vals = f21.read_tokens(d.map_ntheta)
                    for j in range(d.map_ntheta):
                        dbzylu[i, j] = fortran_float(vals[j])
            except Exception:
                print('attention fichier magnetique incomplet...')
                print('bx et by seront nuls...')
                print('vous resterez dans le plan median...')
                dbzxlu[:d.map_nr, :d.map_ntheta] = 0
                dbzylu[:d.map_nr, :d.map_ntheta] = 0

            # Symétries de la carte magnétique
            nr = d.map_nr
            d.bz_map[:nr, 180:91:-1] = d.bz_map[:nr, :89]
            d.bz_map[:nr, 181:361] = d.bz_map[:nr, 1:181]
            d.bz_map[:nr, 361:721] = d.bz_map[:nr, 1:361]
            d.d2bz_map[:nr, 180:91:-1] = d.d2bz_map[:nr, :89]
            d.d2bz_map[:nr, 181:361] = d.d2bz_map[:nr, 1:181]
            d.d2bz_map[:nr, 361:721] = d.d2bz_map[:nr, 1:361]

            r2 = np.sqrt(2.0) / 2.0
            d.dbzdx_map[:nr, 90::-1] = r2 * (dbzxlu[:nr, :91] + dbzylu[:nr, :91])
            d.dbzdy_map[:nr, 90::-1] = r2 * (dbzxlu[:nr, :91] - dbzylu[:nr, :91])
            d.dbzdx_map[:nr, 180:91:-1] = d.dbzdy_map[:nr, :89]
            d.dbzdy_map[:nr, 180:91:-1] = d.dbzdx_map[:nr, :89]
            d.dbzdx_map[:nr, 181:361] = -d.dbzdx_map[:nr, 179::-1]
            d.dbzdy_map[:nr, 181:361] = d.dbzdy_map[:nr, 179::-1]
            d.dbzdx_map[:nr, 361:721] = d.dbzdx_map[:nr, 359::-1]
            d.dbzdy_map[:nr, 361:721] = -d.dbzdy_map[:nr, 359::-1]
            d.map_dtheta = np.pi / 4.0 / (d.map_ntheta - 1.0)
            ntetall = 721

    if d.harmonic_n_coef[0] != 0 or d.harmonic_n_coef[1] != 0:
        print('carte mag. + rajout harmonique 1 et/ou 2...')
        addh(d.map_rmin, d.map_rmax, d.map_nr, ntetall, d.bz_map, d.d2bz_map,
             d.dbzdx_map, d.dbzdy_map, d.harmonic_n_coef, d.harmonic_coef, d.harmonic_phase)

    d.map_dr = (d.map_rmax - d.map_rmin) / (d.map_nr - 1.0)

    # Sans b_factor, la carte est utilisee telle quelle (comme le Fortran).
    # Avec b_factor, elle est ramenee a max|Bz| = b_factor (T).
    if d.b_factor is not None:
        bz_max = np.max(np.abs(d.bz_map[:d.map_nr, :ntetall]))
        if bz_max == 0.0:
            raise ValueError(f'carte magnetique nulle ({d.bfield_file}) : b_factor inapplicable')
        scale = d.b_factor / bz_max
        d.bz_map *= scale
        d.d2bz_map *= scale
        d.dbzdx_map *= scale
        d.dbzdy_map *= scale
        print(f'carte magnetique : Bmax_lu = {bz_max:.6e} T, facteur B = {d.b_factor:.6f} T')

    # Correction relative dB/B
    d.bz_map *= (1.0 + d.db_over_b)
    d.d2bz_map *= (1.0 + d.db_over_b)

    # Lecture potentiel centre
    if d.central_potential_file != 'no':
        if d.swap_bytes == 'y':
            fdtype = np.dtype(np.float64).newbyteorder('>')
            idtype = np.dtype(np.int32).newbyteorder('>')
        else:
            fdtype = np.float64
            idtype = np.int32
        with open(d.central_potential_file, 'rb') as f22:
            # Lecture entête Fortran (record marker)
            rec = f22.read(4)
            header_floats = np.frombuffer(f22.read(32), dtype=fdtype)
            header_ints = np.frombuffer(f22.read(8), dtype=idtype)
            d.central_xmin = float(header_floats[0]) / 1000.0
            d.central_xmax = float(header_floats[1]) / 1000.0
            d.central_ymin = float(header_floats[2]) / 1000.0
            d.central_ymax = float(header_floats[3]) / 1000.0
            nx1 = int(header_ints[0])
            ny1 = int(header_ints[1])
            f22.read(4)  # fin record
            for i in range(nx1):
                f22.read(4)
                row = np.frombuffer(f22.read(ny1 * 8), dtype=fdtype)
                d.central_pot[i, :ny1] = row
                f22.read(4)
            try:
                for i in range(nx1):
                    f22.read(4)
                    row = np.frombuffer(f22.read(ny1 * 8), dtype=fdtype)
                    d.central_d2pot[i, :ny1] = row
                    f22.read(4)
            except Exception:
                print('attention fichier potentiel centre incomplet...')
                print('ez sera nul...')
                d.central_d2pot[:] = 0.0
    else:
        print('pas de fichier centre cyclo...')
        vals = parse_floats(read_input('xmin1 xmax1 ymin1 ymax1 nx1 ny1 ?'))
        d.central_xmin, d.central_xmax, d.central_ymin, d.central_ymax = vals[:4]
        nx1 = int(vals[4])
        ny1 = int(vals[5])
        d.central_pot[:] = 0.0
        d.central_d2pot[:] = 0.0

    d.central_dx = (d.central_xmax - d.central_xmin) / (nx1 - 1.0)
    d.central_dy = (d.central_ymax - d.central_ymin) / (ny1 - 1.0)

    # Lecture potentiel périphérique
    with open(d.outer_potential_file, 'r') as f23_raw:
      f23 = FortranFile(f23_raw)
      with f23:
        vals = f23.read_tokens(6)
        d.outer_xmin = fortran_float(vals[0])
        d.outer_xmax = fortran_float(vals[1])
        d.outer_ymin = fortran_float(vals[2])
        d.outer_ymax = fortran_float(vals[3])
        nx2 = int(fortran_float(vals[4]))
        ny2 = int(fortran_float(vals[5]))
        for i in range(nx2):
            vals = f23.read_tokens(ny2)
            for j in range(ny2):
                d.outer_pot[i, j] = fortran_float(vals[j])
        for i in range(nx2):
            vals = f23.read_tokens(ny2)
            for j in range(ny2):
                d.outer_d2pot[i, j] = fortran_float(vals[j])

    d.outer_dx = (d.outer_xmax - d.outer_xmin) / (nx2 - 1.0)
    d.outer_dy = (d.outer_ymax - d.outer_ymin) / (ny2 - 1.0)

    distribute_bunch(d, sigma_matrices)

    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.theta_deg = np.arccos(d.x[0] / rr)
    if d.y[0] < 0.0:
        d.theta_deg = 2.0 * np.pi - d.theta_deg
    d.theta_deg = 180.0 / np.pi * d.theta_deg
    d.turn = 0
    if d.steps_per_turn < 0:
        d.turn = -d.n_steps // d.steps_per_turn

    if d.z_output_file != 'no':
        files_dict[15] = open(d.z_output_file, 'w')
    if d.xy_output_file != 'no':
        files_dict[16] = open(d.xy_output_file, 'w')


# =============================================================================
# info - Affichage des informations
# =============================================================================
def info(d, iprinf_holder):
    """Affiche les informations de suivi."""
    vec = np.zeros(d.n_particles)
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)

    if iprinf_holder[0] == 1:
        iprinf_holder[0] = 0

        if d.z_output_file != 'no' and 15 in files_dict:
            f15 = files_dict[15]
            f15.write('Programme LIONS\n')
            f15.write('informations selon z\n')

        if d.xy_output_file != 'no' and 16 in files_dict:
            f16 = files_dict[16]
            f16.write('Programme LIONS\n')
            f16.write(f'{d.title}\n')
            f16.write(f' masse atomique        :    {d.atomic_mass[0]:11.7f}\n')
            f16.write(f' A                     :          {d.mass_number[0]:5.1f}\n')
            f16.write(f' q                     :          {d.charge_state[0]:5.1f}\n')
            r0mm = d.r0[0] * 1000.0
            teta0d = d.theta0[0] * 180.0 / np.pi
            z0mm = d.z0[0] * 1000.0
            pr0m = d.pr0[0] * 1000.0
            pz0m = d.pz0[0] * 1000.0
            f16.write(f' r0     : {r0mm:6.1f}\n')
            f16.write(f' teta0  : {teta0d:6.1f}\n')
            f16.write(f' z0     : {z0mm:6.1f}\n')
            f16.write(f' pr0    : {pr0m:8.2f}\n')
            f16.write(f' pz0    : {pz0m:8.2f}\n')
            f16.write(f' w0     : {d.w0_mev_per_u:14.10f}\n')
            f16.write(f' phi0   : {d.phi0_deg[0]:11.2f}\n')
            f16.write(f' fhf    : {d.rf_frequency:9.0f}\n')
            f16.write(f' harm   : {d.harmonic:9.0f}\n')
            f16.write(f' pot_d  : {-d.v_dee_right:9.1f}\n')
            f16.write(f' pot_g  : {-d.v_dee_left:9.1f}\n')
            f16.write(f' db/b   : {d.db_over_b:9.6f}\n')
            if d.b_factor is not None:
                f16.write(f' B fact : {d.b_factor:9.6f} T\n')
            dtnano = d.dt * 1e9
            f16.write(f' nb part: {d.n_particles:14d}\n')
            f16.write(f' nb iter: {d.n_steps:14d}\n')
            f16.write(f' nb i/t : {d.steps_per_turn:14d}\n')
            f16.write(f' dt(ns) : {dtnano:9.6f}\n')
            f16.write('    ite   tr   teta      x       y      xc'
                      '      yc        bz '
                      '    es      er           r      tpr'
                      '       w           dphi\n')

    d.theta_prev_deg = d.theta_deg
    rr = np.sqrt(d.x[0]**2 + d.y[0]**2)
    d.theta_deg = np.arccos(np.clip(d.x[0] / rr, -1.0, 1.0))
    if d.y[0] < 0.0:
        d.theta_deg = 2.0 * np.pi - d.theta_deg
    d.theta_deg = 180.0 / np.pi * d.theta_deg

    if d.steps_per_turn > 0:
        if d.theta_deg < d.theta_prev_deg:
            d.turn += 1
    else:
        if d.theta_deg > d.theta_prev_deg:
            d.turn -= 1

    d.theta_total_deg = d.theta_deg + 360.0 * (d.turn - 1)

    if (d.step % d.output_every != 0 and d.step != 1 and
            d.step != d.n_steps and d.step != (d.n_steps - 1)):
        return

    phihf = 180.0 / np.pi * d.omega_rf * d.t
    dph = phihf - d.harmonic * d.theta_total_deg - 90.0
    wmevpa = d.rest_energy[0] * (d.gamma[0] - 1.0) / 1.60217733e-13 / d.mass_number[0]
    xmm = d.x[0] * 1000.0
    ymm = d.y[0] * 1000.0
    zmm = d.z[0] * 1000.0
    rmm = rr * 1000.0
    pp = np.sqrt(d.px[0]**2 + d.py[0]**2)
    ppp = np.sqrt(d.px[0]**2 + d.py[0]**2 + d.pz[0]**2)
    pr_val = (xmm * d.px[0] + ymm * d.py[0]) / rmm
    tpr = np.arcsin(pr_val / ppp) * 1000.0
    xc1mm = (d.x[0] + d.py[0] * d.dt / d.boris_e_coef[0] / d.bz[0]) * 1000.0
    yc1mm = (d.y[0] - d.px[0] * d.dt / d.boris_e_coef[0] / d.bz[0]) * 1000.0
    eskvcm = (d.ex[0] * d.px[0] + d.ey[0] * d.py[0]) / pp / 100000.0
    erkvcm = (d.ex[0] * d.py[0] - d.ey[0] * d.px[0]) / pp / 100000.0
    ezkvcm = d.ez[0] / 100000.0
    tpz = np.arcsin(d.pz[0] / ppp) * 1000.0

    if d.step % d.output_every != 0 and d.step != 1:
        pass  # skip normal output, go to end-of-run output
    else:
        if d.z_output_file != 'no' and 15 in files_dict:
            if abs(d.z[0]) > d.nuz_z_cut:
                anuz_c = (-1.0 / d.mass[0] / d.gamma[0] * (d.harmonic / d.omega_rf)**2 / d.z[0] *
                          d.charge[0] * (d.ez[0] + C_LIGHT / d.gamma[0] *
                          (d.px[0] * d.by[0] - d.py[0] * d.bx[0])))
                d.nuz = np.sqrt(abs(anuz_c))
            else:
                d.nuz = 0.0
            files_dict[15].write(
                f' {d.step:6d}  {d.turn:3d}  {d.theta_deg:5.1f}'
                f' {xmm:8.2f} {ymm:8.2f} {zmm:8.2f} {rmm:8.2f}'
                f'   {d.bx[0]:7.4f} {d.by[0]:7.4f} {d.bz[0]:7.4f}'
                f' {ezkvcm:7.3f}     {tpz:6.1f} {d.nuz:12.5f}'
                f' {d.sorter_voltage_seen:12.5f}\n')

        if d.xy_output_file != 'no' and 16 in files_dict:
            files_dict[16].write(
                f' {d.step:7d}  {d.turn:3d}  {d.theta_deg:5.1f}'
                f' {xmm:7.1f} {ymm:7.1f} {xc1mm:7.1f} {yc1mm:7.1f}'
                f'   {d.bz[0]:7.4f} {eskvcm:7.3f} {erkvcm:7.3f}'
                f'     {rmm:7.1f}  {tpr:6.1f} {wmevpa:9.5f} {dph:15.2f}\n')

    # Corrélations (sautees avec les sorties z/xy, comme le goto 800 du Fortran)
    if ((d.step % d.output_every == 0 or d.step == 1) and d.correlations != 0
            and not (d.correlation_file == 'no' and d.matrix_file == 'no' and d.dc_matrix_file == 'no')):
        for k in range(d.n_families):
            nbr = d.family_size[k]
            i1 = d.family_start[k]
            i2 = i1 + nbr
            d.p[i1] = np.sqrt(d.px[i1]**2 + d.py[i1]**2 + d.pz[i1]**2)
            rc = -d.mass_c[k] * d.p[i1] / d.charge[i1] / d.bz[i1]
            # Calcul phase space (simplifié, même logique que caldxs)
            d.ps_cartesian[0, :nbr] = d.x[i1:i2] - d.x[i1]
            d.ps_cartesian[4, :nbr] = d.y[i1:i2] - d.y[i1]
            ri1 = np.sqrt(d.x[i1]**2 + d.y[i1]**2)
            tetai1 = np.arccos(np.clip(d.x[i1] / ri1, -1.0, 1.0))
            if d.y[i1] < 0.0:
                tetai1 = 2.0 * np.pi - tetai1
            pri1 = (d.x[i1] * d.px[i1] + d.y[i1] * d.py[i1]) / ri1
            tpri1 = np.arcsin(pri1 / d.p[i1])
            co = np.cos(tetai1 - tpri1)
            si = np.sin(tetai1 - tpri1)
            vec[:nbr] = d.ps_cartesian[0, :nbr] * co + d.ps_cartesian[4, :nbr] * si
            d.ps_cartesian[4, :nbr] = -d.ps_cartesian[0, :nbr] * si + d.ps_cartesian[4, :nbr] * co
            d.ps_cartesian[0, :nbr] = vec[:nbr]
            for ii in range(nbr):
                d.ps[4, ii] = rc * np.arctan(d.ps_cartesian[4, ii] / (rc + d.ps_cartesian[0, ii]))
            vec[:nbr] = np.sin(d.ps[4, :nbr] / 2.0 / rc)**2
            d.ps[0, :nbr] = (d.ps_cartesian[0, :nbr] + 2.0 * rc * vec[:nbr]) / (1.0 - 2.0 * vec[:nbr])
            vec1[:nbr] = np.cos(d.ps[4, :nbr] / rc)
            vec2[:nbr] = np.sin(d.ps[4, :nbr] / rc)
            d.ps_cartesian[1, :nbr] = d.px[i1:i2] - d.px[i1] * vec1[:nbr] + d.py[i1] * vec2[:nbr]
            d.ps_cartesian[5, :nbr] = d.py[i1:i2] - d.px[i1] * vec2[:nbr] - d.py[i1] * vec1[:nbr]
            vec[:nbr] = d.ps_cartesian[1, :nbr] * co + d.ps_cartesian[5, :nbr] * si
            d.ps_cartesian[5, :nbr] = -d.ps_cartesian[1, :nbr] * si + d.ps_cartesian[5, :nbr] * co
            d.ps_cartesian[1, :nbr] = vec[:nbr]
            d.ps[1, :nbr] = d.ps_cartesian[1, :nbr] * vec1[:nbr] + d.ps_cartesian[5, :nbr] * vec2[:nbr]
            d.ps[5, :nbr] = -d.ps_cartesian[1, :nbr] * vec2[:nbr] + d.ps_cartesian[5, :nbr] * vec1[:nbr]
            d.ps[2, :nbr] = d.z[i1:i2] - d.z[i1]
            d.ps[3, :nbr] = d.pz[i1:i2] - d.pz[i1]
            d.ps[5, :nbr] = np.sqrt(
                (d.p[i1] + d.ps[5, :nbr])**2 +
                d.ps[1, :nbr]**2 + d.ps[3, :nbr]**2) - d.p[i1]

            d.sigma[:] = 0.0
            for i in range(6):
                for j in range(6):
                    d.sigma[i, j] = np.sum(d.ps[i, :nbr] * d.ps[j, :nbr] *
                                         d.alive[:nbr].astype(float)) / d.n_alive

            tetai1_deg = tetai1 * 180.0 / np.pi
            ri1mm = ri1 * 1000.0
            d.dc_transform[:] = 0.0
            d.dc_transform[0, 0] = -1.0
            d.dc_transform[1, 1] = -1.0
            for i in range(2, 6):
                d.dc_transform[i, i] = 1.0
            d.dc_transform[0, 5] = ri1 / d.p[i1]
            d.dc_transform[4, 1] = -d.dc_transform[0, 5]
            d.sigma_dc = d.dc_transform @ d.sigma @ d.dc_transform.T

            if d.correlation_file != 'no' and 17 in files_dict:
                gai1 = np.sqrt(1.0 + d.p[i1]**2)
                dr = np.sqrt(d.sigma[0, 0]) * SQRT8
                drmm = dr * 1000.0
                ds = np.sqrt(d.sigma[4, 4]) * SQRT8
                dsmm = ds * 1000.0
                dz_val = np.sqrt(d.sigma[2, 2]) * SQRT8
                dzmm = dz_val * 1000.0
                dphi = ds * d.omega_rf / C_LIGHT / d.p[i1] * 180.0 / np.pi * gai1
                dpr = np.sqrt(d.sigma[1, 1]) * SQRT8
                drdpr = 8.0 * np.sqrt(abs(d.sigma[0, 0] * d.sigma[1, 1] - d.sigma[0, 1]**2))
                epsr = drdpr * 1e6 / d.p[i1]
                dzdpz = 8.0 * np.sqrt(abs(d.sigma[2, 2] * d.sigma[3, 3] - d.sigma[2, 3]**2))
                epsz = dzdpz * 1e6 / d.p[i1]
                dp = np.sqrt(d.sigma[5, 5]) * SQRT8
                dsdp = 8.0 * np.sqrt(abs(d.sigma[4, 4] * d.sigma[5, 5] - d.sigma[4, 5]**2))
                dwsw = dp * d.p[i1] / gai1 / (gai1 - 1.0)
                wmpa = d.rest_energy[k] * (gai1 - 1.0) / 1.60217733e-13 / d.mass_number[k]
                dwmpa = dwsw * wmpa
                dwdphi = (d.omega_rf / C_LIGHT / d.p[i1] * 180.0 / np.pi * gai1 *
                          d.p[i1] / gai1 / (gai1 - 1.0) * wmpa * dsdp * 1000.0)
                files_dict[17].write(
                    f' {d.turn:3d}  {ri1mm:7.2f}  {tetai1_deg:5.1f}'
                    f' {drmm:6.2f} {dsmm:6.2f} {dphi:7.2f}'
                    f' {epsr:6.2f} {dwsw:12.7e} {wmpa:6.3f}'
                    f' {dwmpa:12.7e} {dwdphi:12.7e}'
                    f'  {dzmm:6.2f} {epsz:6.2f}\n')

    # Sauvegarde fin de simulation
    if d.step == (d.n_steps - 1):
        d.end_x = d.x[0]
        d.end_y = d.y[0]
        d.end_z = d.z[0]
        d.end_r = rr
        d.end_theta_deg = d.theta_deg
        d.end_phase_deg = phihf
        if d.bunch_save_file != 'no':
            print('attention : sauvegarde paquet d ions...')
            with open(d.bunch_save_file, 'wb') as f26:
                for i in range(d.n_particles):
                    f26.write(struct.pack('ddd', d.x[i], d.y[i], d.z[i]))

    if d.step == d.n_steps:
        pxf = d.px[0]
        pyf = d.py[0]
        pzf = d.pz[0]
        pf = np.sqrt(pxf**2 + pyf**2 + pzf**2)
        copzf = np.cos(pzf)
        alphaf = np.arccos(pxf / pf / copzf)
        if pyf < 0.0:
            alphaf = 2.0 * np.pi - alphaf
        tprf = ((np.pi / 2.0 + d.end_theta_deg * np.pi / 180.0) - alphaf) * 1000.0
        wf = wmevpa
        tpzf = np.arcsin(pzf / pf) * 1000.0
        if d.bunch_save_file != 'no':
            with open(d.bunch_save_file, 'ab') as f26:
                for i in range(d.n_particles):
                    f26.write(struct.pack('ddd', d.px[i], d.py[i], d.pz[i]))

        if d.xy_output_file != 'no' and 16 in files_dict:
            f16 = files_dict[16]
            f16.write('---------------------------------\n')
            f16.write('info pour le mouvement inverse...\n')
            f16.write(f'r0 teta0 z0 pr0 pz0 w0mev phi0 :\n')
            f16.write(f'{d.end_r} {d.end_theta_deg} {d.end_z} {tprf}\n')
            f16.write(f'{tpzf} {wf} {d.end_phase_deg}\n')
            f16.write('info pour continuer le mouvement ...\n')
            f16.write(f'{rr} {d.theta_deg} {d.z[0]} {tprf}\n')
            f16.write(f'{tpzf} {wf} {phihf}\n')
            f16.write('info utile pour DEDUCE_SIGMA ...\n')
            f16.write(f'amat a qn w0mev : \n')
            f16.write(f'{d.atomic_mass[0]} {d.mass_number[0]} {d.charge_state[0]} {wf}\n')
            f16.write(f'champ magnetique (val. absolue) :\n{-d.bz[0]}\n')
            f16.write(f'fhf harmonique :\n{d.rf_frequency} {d.harmonic}\n')


# =============================================================================
# initgra - Initialisation graphique (bidims)
# =============================================================================
def initgra(d, bidim, numbid=None):
    """Initialise les diagnostics graphiques (bidims). numbid=None : demande a l'utilisateur."""
    d.in_window[:] = True
    d.bidim_every = 1
    if numbid is None:
        numbid = int(read_input('numero premier bidim (0 si aucun bidim) :'))
    d.bidim_first = numbid
    if d.bidim_first == 0:
        return

    d.bidim_every = int(read_input('frequence_de_trace par rapport aux pas de temps:'))
    d.centre_mode = read_input('trace des centres ou des barycentres (c ou b) :')
    d.bary_size = abs(d.steps_per_turn // d.bidim_every)
    d.bary_count = 0

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
    if d.centre_mode == 'c':
        titre[1] = 'centres de courbure'
    else:
        titre[1] = 'barycentres'
    labelx[1] = 'x'
    labely[1] = 'y'

    d.bidim_start_step = 1
    if d.n_particles > 1:
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
        print(f'itemax : {d.n_steps}')
        d.bidim_start_step = int(read_input('debut_iteration_stockage_bidim (entre 1 et itemax) :'))

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

    if d.n_particles > 1:
        nbr_bid = 17
    else:
        nbr_bid = 5

    ier = bidim.openbi(nbr_bid, d.bidim_first, nbx, nby, titre, labelx, labely, bornes)
    if ier != 0:
        print(f'bidim : openbi : stop erreur : {ier}')
        sys.exit(1)

    d.bidim_ax1[0] = 8191.0 / (xbid[0, 1] - xbid[0, 0])
    d.bidim_ax0[0] = -d.bidim_ax1[0] * xbid[0, 0]
    d.bidim_ay1[0] = 8191.0 / (xbid[1, 1] - xbid[1, 0])
    d.bidim_ay0[0] = -d.bidim_ay1[0] * xbid[1, 0]
    d.bidim_ax1[1] = 8191.0 / (xbamax - xbamin)
    d.bidim_ax0[1] = -d.bidim_ax1[1] * xbamin
    d.bidim_ay1[1] = 8191.0 / (ybamax - ybamin)
    d.bidim_ay0[1] = -d.bidim_ay1[1] * ybamin

    if d.n_particles > 1:
        d.bidim_ax1[3] = 8191.0 / 2.0 / dphimax
        d.bidim_ax0[3] = d.bidim_ax1[3] * dphimax
        d.bidim_ay1[3] = 8191.0 / 2.0 / dwswmax
        d.bidim_ay0[3] = d.bidim_ay1[3] * dwswmax
        d.bidim_ax1[4] = 8191.0 / 2.0 / drmax
        d.bidim_ax0[4] = d.bidim_ax1[4] * drmax
        d.bidim_ay1[4] = 8191.0 / 2.0 / dprsp0max
        d.bidim_ay0[4] = d.bidim_ay1[4] * dprsp0max
        d.bidim_ax1[5] = 8191.0 / 2.0 / drmax
        d.bidim_ax0[5] = d.bidim_ax1[5] * drmax
        d.bidim_ay1[5] = 8191.0 / 2.0 / dphimax
        d.bidim_ay0[5] = d.bidim_ay1[5] * dphimax
        d.bidim_ax1[6] = 8191.0 * d.bidim_every / d.n_steps
        d.bidim_ax0[6] = 0.0
        d.bidim_ay1[6] = 8191.0 / 2.0 / drmax
        d.bidim_ay0[6] = d.bidim_ay1[6] * drmax
        d.bidim_ax1[7] = 8191.0 / 2.0 / dzmax
        d.bidim_ax0[7] = d.bidim_ax1[7] * dzmax
        d.bidim_ay1[7] = 8191.0 / 2.0 / dpzsp0max
        d.bidim_ay0[7] = d.bidim_ay1[7] * dpzsp0max
        d.bidim_ax1[8] = 8191.0 * d.bidim_every / d.n_steps
        d.bidim_ax0[8] = 0.0
        d.bidim_ay1[8] = 8191.0 / 2.0 / dzmax
        d.bidim_ay0[8] = d.bidim_ay1[8] * dzmax

    # Lecture repères
    ficrep = read_input('lecture du fichier_des_reperes (lignes, piliers etc...) :')
    d.n_traps = 0
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
                                ix = int(d.bidim_ax1[0] * xd / 1000.0 + d.bidim_ax0[0])
                                iy = int(d.bidim_ay1[0] * yd / 1000.0 + d.bidim_ay0[0])
                                bidim.incbid(d.bidim_first, ix, iy)
                    elif typcou == 'cercle':
                        nbrcer = int(f24.readline().split()[0])
                        for num in range(nbrcer):
                            vals = f24.readline().split()
                            xce, yce, rce = fortran_float(vals[0]), fortran_float(vals[1]), fortran_float(vals[2])
                            tetce1 = fortran_float(vals[3]) * np.pi / 180.0
                            tetce2 = fortran_float(vals[4]) * np.pi / 180.0
                            npcer = int(vals[5])
                            flacer = fortran_float(vals[6])
                            dtetce = (tetce2 - tetce1) / (nbpo - 1.0)
                            if flacer < -10.0:
                                d.n_traps += 1
                                d.trap_x[d.n_traps - 1] = xce / 1000.0
                                d.trap_y[d.n_traps - 1] = yce / 1000.0
                                d.trap_r[d.n_traps - 1] = rce / 1000.0
                            if flacer > -9000.0:
                                for kk in range(nbpo):
                                    tetce = tetce1 + kk * dtetce
                                    xd = xce + rce * np.cos(tetce)
                                    yd = yce + rce * np.sin(tetce)
                                    ix = int(d.bidim_ax1[0] * xd / 1000.0 + d.bidim_ax0[0])
                                    iy = int(d.bidim_ay1[0] * yd / 1000.0 + d.bidim_ay0[0])
                                    bidim.incbid(d.bidim_first, ix, iy)
                    else:
                        print(f'erreur lecture fichier reperes... typcou : {typcou}')
                        sys.exit(1)
        except Exception as e:
            print(f'erreur en lecture du fichier reperes: {e}')

    if d.n_traps != 0:
        print(f'attention : il y a des pieges a ions...')
        print(f'nombre de pieges : {d.n_traps}')
        d.trap_last_step = int(read_input('iteration max pour les pieges:'))
        print('flag de choix particules piegees: (0 ou 1)')
        print('0   : reste scotchee dans les pieges')
        print('1   : continue, confondue avec la particule centrale')
        d.trapped_follow = int(read_input(''))


# =============================================================================
# bid - Enregistrement des bidims
# =============================================================================
def bid(d, bidim):
    """Enregistre les données dans les bidims."""
    vec1 = np.zeros(d.n_particles)
    vec2 = np.zeros(d.n_particles)
    vec3 = np.zeros(d.n_particles)
    vec4 = np.zeros(d.n_particles)
    vec5 = np.zeros(d.n_particles)
    vec6 = np.zeros(d.n_particles)
    if d.step < d.bidim_start_step:
        return

    for k in range(d.n_particles):
        if not d.in_window[k]:
            continue
        ix = int(d.bidim_ax1[0] * d.x[k] + d.bidim_ax0[0])
        iy = int(d.bidim_ay1[0] * d.y[k] + d.bidim_ay0[0])
        if 0 <= ix <= 8191 and 0 <= iy <= 8191:
            if k == 0:
                bidim.incbid(d.bidim_first, ix, iy)
                if d.n_particles > 1:
                    bidim.incbid(d.bidim_first + 2, ix, iy)
            else:
                bidim.incbid(d.bidim_first + 2, ix, iy)
        if k == 0:
            xck = d.x[k] + d.py[k] * d.dt / d.boris_e_coef[k] / d.bz[k]
            yck = d.y[k] - d.px[k] * d.dt / d.boris_e_coef[k] / d.bz[k]
            if d.centre_mode == 'c':
                ix2 = int(d.bidim_ax1[0] * xck + d.bidim_ax0[0])
                iy2 = int(d.bidim_ay1[0] * yck + d.bidim_ay0[0])
                if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                    bidim.incbid(d.bidim_first + 1, ix2, iy2)
            if d.centre_mode == 'b':
                d.bary_count += 1
                d.bary_x[d.bary_count - 1] = xck
                d.bary_y[d.bary_count - 1] = yck
                if d.bary_count == d.bary_size:
                    d.bary_count = 0
                    xbary = np.sum(d.bary_x[:d.bary_size]) / d.bary_size
                    ybary = np.sum(d.bary_y[:d.bary_size]) / d.bary_size
                    ix2 = int(d.bidim_ax1[1] * xbary + d.bidim_ax0[1])
                    iy2 = int(d.bidim_ay1[1] * ybary + d.bidim_ay0[1])
                    if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                        bidim.incbid(d.bidim_first + 1, ix2, iy2)

    if d.n_particles > 1 and d.n_families == 1:
        caldxs(d)
        k = 0
        nbr = d.family_size[k]
        i1 = d.family_start[k]
        gai1 = np.sqrt(1.0 + d.p[i1]**2)
        vec1[:nbr] = d.ps[4, :nbr] * d.omega_rf / C_LIGHT / d.p[i1] * 180.0 / np.pi * gai1
        vec2[:nbr] = d.ps[5, :nbr] * d.p[i1] / gai1 / (gai1 - 1.0)
        vec3[:nbr] = d.ps[0, :nbr] * 1000.0
        vec4[:nbr] = d.ps[1, :nbr] / d.p[i1] * 1000.0
        vec5[:nbr] = d.ps[2, :nbr] * 1000.0
        vec6[:nbr] = d.ps[3, :nbr] / d.p[i1] * 1000.0

        for i in range(nbr):
            if not d.in_window[i]:
                continue
            ix = int(d.bidim_ax1[3] * vec1[i] + d.bidim_ax0[3])
            iy = int(d.bidim_ay1[3] * vec2[i] + d.bidim_ay0[3])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 3, ix, iy)

            ix = int(d.bidim_ax1[4] * vec3[i] + d.bidim_ax0[4])
            iy = int(d.bidim_ay1[4] * vec4[i] + d.bidim_ay0[4])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 4, ix, iy)

            ix = int(d.bidim_ax1[7] * vec5[i] + d.bidim_ax0[7])
            iy = int(d.bidim_ay1[7] * vec6[i] + d.bidim_ay0[7])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 7, ix, iy)

        for i in range(nbr):
            if not d.in_window[i]:
                continue
            ix = int(d.bidim_ax1[5] * vec3[i] + d.bidim_ax0[5])
            iy = int(d.bidim_ay1[5] * vec1[i] + d.bidim_ay0[5])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 5, ix, iy)

        ix = int(d.bidim_ax1[6] * d.step / d.bidim_every + d.bidim_ax0[6])
        ix = max(0, min(ix, 8191))
        for i in range(nbr):
            if not d.in_window[i]:
                continue
            iy = int(d.bidim_ay1[6] * vec3[i] + d.bidim_ay0[6])
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 6, ix, iy)
            iy = int(d.bidim_ay1[8] * vec5[i] + d.bidim_ay0[8])
            iy = max(0, min(iy, 8191))
            bidim.incbid(d.bidim_first + 8, ix, iy)


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
    d.sorter_rmin = rmin
    d.sorter_rmax = rmax
    d.sorter_theta1 = teta1 * np.pi / 180.0
    d.sorter_theta2 = teta2 * np.pi / 180.0
    d.sorter_voltage = pot
    d.sorter_gap = h
    d.sorter_phase = phi * np.pi / 180.0
    d.sorter_nuz = anuz
    d.nuz_z_cut = coupnuz
    d.sorter_dv = dpotspot


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
    d.step = 0
    boite(d)
    calcbe(d)
    iprinf_holder = [1]
    info(d, iprinf_holder)
    try:
        if on_step is None or on_step(d) is not False:
            for ite in range(1, d.n_steps + 1):
                d.step = ite
                piege(d)
                boite(d)
                calcbe(d)
                pousse(d)
                d.t = ite * d.dt + d.t0
                info(d, iprinf_holder)
                if d.bidim_first != 0 and (ite % d.bidim_every == 0 or ite == 1):
                    bid(d, bidim)
                if on_step is not None and on_step(d) is False:
                    break
    finally:
        if d.bidim_first != 0:
            bidim.clobid()
        for f in files_dict.values():
            f.close()
        files_dict.clear()


def run(ficin, swap=False, vertical=None, on_step=None):
    """Simulation sans interaction ni bidim. Renvoie l'etat final (State).

    - ficin : fichier de donnees ; ses chemins sont relatifs au repertoire courant.
    - swap : inverser les octets du potentiel central (fichiers big-endian).
    - vertical : dict de parametres pour set_vertical_sorter, None = pas de tri.
    - on_step : voir track().
    """
    d = State()
    d.swap_bytes = 'y' if swap else 'n'
    init(d, ficin)
    set_vertical_sorter(d, **(vertical or {}))
    initgra(d, None, numbid=0)
    track(d, on_step=on_step)
    return d


def write_trapped(d):
    """Propose d'ecrire les numeros des particules piegees."""
    if d.n_traps != 0:
        nreste = np.count_nonzero(d.alive)
        print(f'nombre de particules non piegees : {nreste}')
        print('ecriture fichier des numero particules piegees ? (o/n)')
        reponse = input().strip()
        if reponse == 'o':
            print('choisissez nom du fichier : ')
            ficpie = input().strip()
            with open(ficpie, 'w') as f60:
                npi = 0
                for k in range(d.n_particles):
                    if not d.alive[k]:
                        npi += 1
                        f60.write(f' {npi:5d} {k+1:5d}\n')


def lions():
    """Programme principal LIONS (interactif)."""
    banner()
    d = State()
    d.swap_bytes = input().strip()
    init(d, read_input('fichier_de_donnees de depart :'))
    read_vertical_sorter(d)
    bidim = BidimManager()
    initgra(d, bidim)
    track(d, bidim)
    write_trapped(d)


if __name__ == '__main__':
    lions()
