"""Deck reader and writer: round trips, every optional section, run from a Deck object."""
from pathlib import Path

import numpy as np
import pytest

import synthetic as syn
from lions.deck import Deck, Deflector, Dipole, Family, FieldHarmonic, Grid, read_deck, write_deck

EXAMPLES = sorted((Path(__file__).resolve().parent.parent / 'examples').glob('don.*'))


def full_deck():
    """Every optional section of the format set."""
    fam = dict(atomic_mass=1.00782503, mass_number=1.0, charge_state=1.0, r0=0.1445,
               theta0_deg=0.0, z0=0.002, pr0_mrad=10.0, pz0_mrad=5.0, w0_mev_per_u=1.0,
               phi0_deg=-250.0)
    return Deck(
        title='all sections', steps_per_turn=360, design_r_final=1.5, design_turns=250,
        rf_frequency=14.4e6, rf_voltage=-51890.0, dee_half_angle_deg=20.0, harmonic=3.0,
        db_over_b=1e-4, phase_right_deg=1.5, phase_left_deg=-2.0, v_dee_right=51890.0,
        v_dee_left=50000.0, zmin=-0.1, zmax=0.1, central_rmax=0.164,
        families=[Family(size=3, **fam), Family(size=2, **fam)], n_particles=5,
        n_steps=3600, output_every=30, z_output_file='z.out', xy_output_file='out/xy.out',
        bfield_file='maps/b field.dat', central_potential_file='pot1.bin',
        outer_potential_file='pot2.dat',
        harmonic1=FieldHarmonic([1e-3, 2e-3], 45.0), harmonic2=FieldHarmonic([3e-4], 10.0),
        gap_offset=0.01, extraction=-1.0, cut_bunch=12.0, grids_on=1.0, grid_start=100.0,
        deflectors=[Deflector(1.513, 6.99e6, 0.0, 0.0, 0.0, 0.0, 250.0, 267.0)],
        dipoles=[Dipole(0.15, 0.5, 0.0, 0.0, 120.0, 180.0)],
        grids=[Grid(200.0, -939.69, -342.02, 0.0, 0.0)], grid_file_prefix='grid',
        correlations=1, sigma_file='sigma.dat', phase_space='6d', correlation_file='cor.out',
        matrix_file='no', dc_matrix_file='no', b_factor=1.56,
        bunch_load_file='bunch.bin', bunch_save_file='final.bin')


@pytest.mark.parametrize('path', EXAMPLES, ids=lambda p: p.name)
def test_example_round_trip(tmp_path, path):
    deck = read_deck(path)
    write_deck(deck, tmp_path / 'a.dat')
    assert read_deck(tmp_path / 'a.dat') == deck
    write_deck(read_deck(tmp_path / 'a.dat'), tmp_path / 'b.dat')
    assert (tmp_path / 'a.dat').read_text() == (tmp_path / 'b.dat').read_text()


def test_all_sections_round_trip(tmp_path):
    deck = full_deck()
    write_deck(deck, tmp_path / 'deck.dat')
    assert read_deck(tmp_path / 'deck.dat') == deck


def test_example_values():
    deck = read_deck(EXAMPLES[0])
    assert (deck.steps_per_turn, deck.rf_frequency, deck.harmonic) == (1080, 14.4e6, 3.0)
    assert deck.families[0].size == 100 and deck.families[0].charge_state == 22.0
    assert [e.r for e in deck.deflectors] == [1.513]
    assert deck.bfield_file == 'cime700rec.1804.pbob1.lio.dat'
    assert deck.b_factor is None and deck.bunch_save_file == 'no'


def test_run_from_deck_object(workdir):
    syn.write_bmap('bmap.dat', -1.0)
    syn.write_deck('deck.dat', fhf=syn.revolution_frequency(1.0, -1.0),
                   r0=syn.larmor_radius(1.0, -1.0), w0=1.0, turns=1)
    from_file = syn.quiet_run('deck.dat')
    from_object = syn.quiet_run(read_deck('deck.dat'))
    for name in ('x', 'y', 'z', 'px', 'py', 'pz'):
        np.testing.assert_array_equal(getattr(from_object.bunch, name), getattr(from_file.bunch, name))
