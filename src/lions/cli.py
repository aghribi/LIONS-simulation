"""Tracking loop, scripted run and the interactive program."""


import numpy as np

from .state import State
from .apertures import apply_traps, check_box
from .bidim import BidimWriter, fill_bidims, setup_bidims
from .diagnostics import write_outputs
from .fields import compute_fields
from .io import load, parse_floats, prompt
from .pusher import boris_push


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
    check_box(d)
    compute_fields(d)
    write_outputs(d)
    try:
        if on_step is None or on_step(d) is not False:
            for ite in range(1, d.n_steps + 1):
                d.step = ite
                apply_traps(d)
                check_box(d)
                compute_fields(d)
                boris_push(d)
                d.t = ite * d.dt + d.t0
                write_outputs(d)
                if d.bidim_first != 0 and (ite % d.bidim_every == 0 or ite == 1):
                    fill_bidims(d, bidim)
                if on_step is not None and on_step(d) is False:
                    break
    finally:
        if d.bidim_first != 0:
            bidim.clobid()
        for f in d.outputs.values():
            f.close()
        d.outputs.clear()


def run(ficin, swap=False, vertical=None, on_step=None):
    """Simulation sans interaction ni bidim. Renvoie l'etat final (State).

    - ficin : fichier de donnees ; ses chemins sont relatifs au repertoire courant.
    - swap : inverser les octets du potentiel central (fichiers big-endian).
    - vertical : dict de parametres pour set_vertical_sorter, None = pas de tri.
    - on_step : voir track().
    """
    d = State()
    d.swap_bytes = 'y' if swap else 'n'
    load(d, ficin)
    set_vertical_sorter(d, **(vertical or {}))
    setup_bidims(d, None, numbid=0)
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


def main():
    """Programme principal LIONS (interactif)."""
    banner()
    d = State()
    d.swap_bytes = input().strip()
    load(d, prompt('fichier_de_donnees de depart :'))
    read_vertical_sorter(d)
    bidim = BidimWriter()
    setup_bidims(d, bidim)
    track(d, bidim)
    write_trapped(d)
