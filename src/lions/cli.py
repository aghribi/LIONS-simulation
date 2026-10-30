"""Tracking loop, scripted run and the interactive program."""


import numpy as np

from .state import State
from .apertures import apply_traps, check_box
from .bidim import BidimWriter, fill_bidims, setup_bidims
from .diagnostics import write_outputs
from .fields import compute_fields
from .fortran_text import parse_floats
from .io import load, prompt
from .pusher import boris_push


def banner():
    print('----------------------------------------------------')
    print('LIONS       : Python (converted from Fortran 90)')
    print('author      : Patrick Bertrand - GANIL')
    print('version  1  : December 1998 to December 2004')
    print('version  2  : March 2017, seed and byte swap changes')
    print('Python port : 2026')
    print('----------------------------------------------------')
    print()
    print('byte swap (March 2017):')
    print('swap bytes? (y/n) (y for old unformatted data files)')


def set_vertical_sorter(d, rmin=1000.0, rmax=2000.0, theta1=0.0, theta2=0.0,
                        voltage=0.0, gap=0.0, phase=0.0, nuz=0.0, nuz_z_cut=0.0, dv=0.0):
    """Vertical sorter: region (m, deg), drive voltage (V), gap (m), phase (deg).

    gap = 0: no vertical sorting. Prompt order of the Fortran: rmin rmax theta1 theta2,
    then voltage gap phase nuz nuz_z_cut dv (pot_ver h_ver phi_ver anuz coupnuz dpotspot).
    """
    d.sorter_rmin = rmin
    d.sorter_rmax = rmax
    d.sorter_theta1 = theta1 * np.pi / 180.0
    d.sorter_theta2 = theta2 * np.pi / 180.0
    d.sorter_voltage = voltage
    d.sorter_gap = gap
    d.sorter_phase = phase * np.pi / 180.0
    d.sorter_nuz = nuz
    d.nuz_z_cut = nuz_z_cut
    d.sorter_dv = dv


def read_vertical_sorter(d):
    """Read the vertical sorter parameters from stdin."""
    print('--------------------')
    print('vertical sorter test:')
    print('--------------------')
    print('rmin_ver rmax_ver theta1_ver theta2_ver:')
    print('1000 2000 0 0 for no vertical sorting')
    geom = parse_floats(input())
    print('pot_ver h_ver phi_ver anuz coupnuz dpotspot:')
    print('0 0 0 0 0 0 for no vertical sorting')
    set_vertical_sorter(d, *geom[:4], *parse_floats(input())[:6])
    print('-------------------')


def track(d, bidim=None, on_step=None):
    """Tracking loop: initial step, then n_steps steps. Closes the output files.

    on_step(d) is called after each step; tracking stops if it returns False.
    """
    d.step = 0
    check_box(d)
    compute_fields(d)
    write_outputs(d)
    try:
        if on_step is None or on_step(d) is not False:
            for step in range(1, d.n_steps + 1):
                d.step = step
                apply_traps(d)
                check_box(d)
                compute_fields(d)
                boris_push(d)
                d.t = step * d.dt + d.t0
                write_outputs(d)
                if d.bidim_first != 0 and (step % d.bidim_every == 0 or step == 1):
                    fill_bidims(d, bidim)
                if on_step is not None and on_step(d) is False:
                    break
    finally:
        if d.bidim_first != 0:
            bidim.close()
        for f in d.outputs.values():
            f.close()
        d.outputs.clear()


def run(deck, swap=False, vertical=None, on_step=None):
    """Run without prompts or bidims. Returns the final State.

    - deck: Deck, or path of an input deck; its file paths are relative to the current directory.
    - swap: swap the bytes of the central potential (big-endian files).
    - vertical: keyword arguments of set_vertical_sorter, None = no vertical sorting.
    - on_step: see track().
    """
    d = State()
    d.swap_bytes = 'y' if swap else 'n'
    load(d, deck)
    set_vertical_sorter(d, **(vertical or {}))
    setup_bidims(d, None, numbid=0)
    track(d, on_step=on_step)
    return d


def write_trapped(d):
    """Offer to write the indices of the trapped particles."""
    if d.n_traps != 0:
        n_free = np.count_nonzero(d.alive)
        print(f'particles not trapped: {n_free}')
        print("write the trapped particle numbers to a file? ('o' for yes, 'n' for no)")
        answer = input().strip()
        if answer == 'o':
            print('file name:')
            path = input().strip()
            with open(path, 'w') as f:
                n_trapped = 0
                for k in range(d.n_particles):
                    if not d.alive[k]:
                        n_trapped += 1
                        f.write(f' {n_trapped:5d} {k+1:5d}\n')


def main():
    """Interactive LIONS program."""
    banner()
    d = State()
    d.swap_bytes = input().strip()
    load(d, prompt('input deck file:'))
    read_vertical_sorter(d)
    bidim = BidimWriter()
    setup_bidims(d, bidim)
    track(d, bidim)
    write_trapped(d)
