"""2D histograms (bidims) in the GANIL binary format."""

import struct
import sys

import numpy as np

from .constants import C_LIGHT
from .beam import phase_space_coords
from .io import fortran_float, parse_floats, prompt


class BidimWriter:
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


def setup_bidims(d, bidim, numbid=None):
    """Initialise les diagnostics graphiques (bidims). numbid=None : demande a l'utilisateur."""
    d.in_window[:] = True
    d.bidim_every = 1
    if numbid is None:
        numbid = int(prompt('numero premier bidim (0 si aucun bidim) :'))
    d.bidim_first = numbid
    if d.bidim_first == 0:
        return

    d.bidim_every = int(prompt('frequence_de_trace par rapport aux pas de temps:'))
    d.centre_mode = prompt('trace des centres ou des barycentres (c ou b) :')
    d.bary_size = abs(d.steps_per_turn // d.bidim_every)
    d.bary_count = 0

    vals = parse_floats(prompt(
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
        vals = parse_floats(prompt('bidim dphi-dw/w : dphi_max dw/w_max :'))
        dphimax = vals[0]
        dwswmax = vals[1]
        vals = parse_floats(prompt('bidim dr-dpr/p0 : dr_max dpr/p0_max :'))
        drmax = vals[0]
        dprsp0max = vals[1]
        vals = parse_floats(prompt('bidim dz-dpz/p0 : dz_max dpz/p0_max :'))
        dzmax = vals[0]
        dpzsp0max = vals[1]
        print(f'itemax : {d.n_steps}')
        d.bidim_start_step = int(prompt('debut_iteration_stockage_bidim (entre 1 et itemax) :'))

    vals = parse_floats(prompt('valeurs de la fenetre des barycentres en metres:'))
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
    ficrep = prompt('lecture du fichier_des_reperes (lignes, piliers etc...) :')
    d.n_traps = 0
    if ficrep != 'no':
        try:
            with open(ficrep, 'r') as f24:
                for _ in range(3):  # display window, grid size, niter omer: unused
                    f24.readline()
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
        print('attention : il y a des pieges a ions...')
        print(f'nombre de pieges : {d.n_traps}')
        d.trap_last_step = int(prompt('iteration max pour les pieges:'))
        print('flag de choix particules piegees: (0 ou 1)')
        print('0   : reste scotchee dans les pieges')
        print('1   : continue, confondue avec la particule centrale')
        d.trapped_follow = int(prompt(''))


def fill_bidims(d, bidim):
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
        phase_space_coords(d)
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
