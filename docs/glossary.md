# Fortran glossary

Python names of the LIONS state (`lions.state.State`, ex `LionData`) and their Fortran
counterparts in `fortran/lions_2.0_ep.f90`. Units are those held in memory (SI, angles in
rad) unless stated. Per-particle arrays have length `n_particles`; per-family arrays are
indexed by family.

## Run and files

| Fortran | Python | Meaning |
|---|---|---|
| `corps` | `title` | Deck title |
| `ficin` | `deck_file` | Input deck |
| `ficb` | `bfield_file` | Magnetic field map |
| `ficv1` | `central_potential_file` | Central-region potential (binary) |
| `ficv2` | `outer_potential_file` | Outer-region potential (text) |
| `ficz`, `ficxy` | `z_output_file`, `xy_output_file` | Text outputs, `no` = off |
| `ficsig` | `sigma_file` | Initial 6D sigma matrices |
| `espace` | `phase_space` | `6d` or `3d` |
| `ficcor`, `ficmat`, `ficdc` | `correlation_file`, `matrix_file`, `dc_matrix_file` | Correlation outputs |
| `ficloa`, `ficsto` | `bunch_load_file`, `bunch_save_file` | Bunch files |
| `ficdebgex`, `ficgex` | `grid_file_prefix`, `grid_files` | Extraction grid crossing files |
| `repswap` | `swap_bytes` | Byte swap of the central potential (bool) |
| `npaspt` | `steps_per_turn` | Time steps per turn (dt = harm / fhf / npaspt); < 0 tracks backwards |
| `nbdt` | `n_steps` | Number of time steps |
| `ifreq` | `output_every` | Output period (steps) |
| `ite` | `step` | Current step |
| `t0`, `t`, `dt` | `t0`, `t`, `dt` | Time (s) |
| `cdt` | `c_dt` | c dt (m) |
| `numtou` | `turn` | Turn counter of particle 0 |
| `tetnew`, `tetold`, `tetall` | `theta_deg`, `theta_prev_deg`, `theta_total_deg` | Azimuth of particle 0 (deg), cumulated |
| `r1th`, `ntours` | `design_r_final`, `design_turns` | Design final radius (m) and turn count |
| `iflcor` | `correlations` | Correlated initial bunch (bool) |

## RF and magnet

| Fortran | Python | Meaning |
|---|---|---|
| `fhf`, `omhf`, `omr` | `rf_frequency`, `omega_rf`, `omega_rev` | RF frequency (Hz), angular RF and revolution frequencies |
| `harm`, `heps` | `harmonic`, `harmonic_sign` | Harmonic number; -1 if even, +1 if odd |
| `pothfd`, `pothfg` | `v_dee_right`, `v_dee_left` | Dee voltages (V), sign as used in the field |
| `dphid`, `dphig` | `phase_right`, `phase_left` | Dee phase offsets (rad) |
| `teta1`, `yomega` | `dee_half_angle`, `gap_offset` | Dee half-angle (rad), skewed-gap offset (m) |
| `dbsb`, `b_factor` | `db_over_b`, `b_factor` | Relative field correction; peak field override (T) or None |
| `tbz`, `d2bz`, `dbzx`, `dbzy` | `bz_map`, `d2bz_map`, `dbzdx_map`, `dbzdy_map` | Median-plane maps on (r, theta) |
| `nr`, `ntet`, `rmin`, `rmax`, `hr`, `htet` | `map_nr`, `map_ntheta`, `map_rmin`, `map_rmax`, `map_dr`, `map_dtheta` | Map grid |
| `ncoh`, `coh`, `phihar` | `harmonic_n_coef`, `harmonic_coef`, `harmonic_phase` | First and second field harmonics added to the map |
| `pot1`, `d2pot1` | `central_pot`, `central_d2pot` | Central-region potential and its z curvature |
| `xmin1`, `xmax1`, `ymin1`, `ymax1`, `hx1`, `hy1` | `central_xmin`, `central_xmax`, `central_ymin`, `central_ymax`, `central_dx`, `central_dy` | Central grid (m) |
| `rmax1` | `central_rmax` | Central / outer region boundary (m) |
| `pot2`, `d2pot2` | `outer_pot`, `outer_d2pot` | Outer-region normalised gap potential |
| `xmin2`, `xmax2`, `ymin2`, `ymax2`, `hx2`, `hy2` | `outer_xmin`, `outer_xmax`, `outer_ymin`, `outer_ymax`, `outer_dx`, `outer_dy` | Outer grid, gap frame (m) |
| `zmin`, `zmax` | `zmin`, `zmax` | Vertical box (m) |

## Particles

| Fortran | Python | Meaning |
|---|---|---|
| `n` | `n_particles` | Number of particles (`np` merged into it) |
| `nbrfam`, `nbrion`, `index` | `n_families`, `family_size`, `family_start` | Families |
| `amat`, `a`, `qn` | `atomic_mass`, `mass_number`, `charge_state` | Per family (u, -, e) |
| `r0`, `teta0`, `z0`, `pr0`, `pz0`, `w0`, `phi0` | `r0`, `theta0`, `z0`, `pr0`, `pz0`, `w0`, `phi0_deg` | Initial conditions per family (m, rad, m, rad, rad, J, deg) |
| `w0mev1` | `w0_mev_per_u` | Initial energy of family 0 (MeV/u) |
| `x`, `y`, `z`, `px`, `py`, `pz` | same | Positions (m), momenta (beta gamma) |
| `ga`, `p`, `r`, `teta` | `gamma`, `p`, `r`, `theta` | Per particle |
| `am0`, `q`, `am0c`, `am0c2` | `mass`, `charge`, `mass_c`, `rest_energy` | Per particle (kg, C, kg m/s, J) |
| `qdts2m0`, `qdtsm0c` | `boris_b_coef`, `boris_e_coef` | q dt / 2m0 and q dt / m0 c |
| `bx`, `by`, `bz`, `ex`, `ey`, `ez` | same | Fields at the particles (T, V/m) |
| `ok` | `alive` | Not trapped |
| `box`, `centre`, `periph` | `in_box`, `in_central`, `in_outer` | Region masks |
| `extra`, `sortant` | `in_element`, `extracted` | Inside the current extraction element; extracted |
| `franchi`, `indgex` | `crossed_grid`, `grid_reached` | Grid crossing at this step; last grid crossed |
| `fenlog` | `in_window` | Diagnostic window |
| `nbrok`, `nbox`, `nbox1`, `nbox2`, `nextra` | `n_alive`, `n_in_box`, `n_central`, `n_outer`, `n_in_element` | Counts |

## Phase space and correlations

| Fortran | Python | Meaning |
|---|---|---|
| `dxs` | `ps` | Phase-space deviations from particle 0: r, pr, z, pz, s, dp |
| `dxd` | `ps_cartesian` | Same, Cartesian frame |
| `dxsgri`, `rgri`, `tprgri`, `phigri`, `xgri`, `ygri` | `ps_projected`, `proj_offset`, `proj_angle`, `proj_phase`, `proj_x`, `proj_y` | Projection on the radial line through particle 0 |
| `cor`, `cordc`, `t1` | `sigma`, `sigma_dc`, `dc_transform` | 6x6 sigma matrix, with dispersion correction |

## Extraction

| Fortran | Python | Meaning |
|---|---|---|
| `rex` | `extraction` | Deck flag: < 0 enables deflectors, dipoles, grids |
| `tet1ex`, `iteext` | `cut_bunch`, `extraction_start_step` | Bunch cut at the deflector radius, or whole bunch from a step |
| `tet2ex`, `eex`, `itegex` | `grids_on`, -, `grid_start_step` | Grids enabled, from a step |
| `ndf`, `rdf`, `edf` | `n_deflectors`, `deflector_r`, `deflector_e` | Electrostatic deflectors (m, V/m) |
| `xdfi`, `ydfi`, `tetdfi`, `uxdfi`, `uydfi` | `deflector_in_x`, `deflector_in_y`, `deflector_in_angle`, `deflector_in_ux`, `deflector_in_uy` | Entrance line |
| `xdfo`, `ydfo`, `tetdfo`, `uxdfo`, `uydfo` | `deflector_out_x`, ... | Exit line |
| `ndm`, `rdm`, `bdm`, `xdm`, `ydm` | `n_dipoles`, `dipole_r`, `dipole_gradient`, `dipole_x`, `dipole_y` | Magnetic dipoles (m, T/m) |
| `tetdmi`, `uxdmi`, `uydmi`, `tetdmo`, ... | `dipole_in_angle`, `dipole_in_ux`, `dipole_in_uy`, `dipole_out_angle`, ... | Entrance and exit lines |
| `ngex`, `tetagex` | `n_grids`, `grid_angle_deg` | Extraction grids |
| `xgex`, `ygex`, `xgexc`, `ygexc`, `drgex`, `uxgex`, `uygex` | `grid_x`, `grid_y`, `grid_xc`, `grid_yc`, `grid_length`, `grid_ux`, `grid_uy` | Grid lines (m) |

## Traps, vertical sorter, bidims

| Fortran | Python | Meaning |
|---|---|---|
| `nbrpie`, `xpie`, `ypie`, `rpie` | `n_traps`, `trap_x`, `trap_y`, `trap_r` | Trap circles (m) |
| `itepiem`, `iflpie` | `trap_last_step`, `trapped_follow` | Traps active up to a step; 1 = trapped particles follow particle 0 |
| `rmin_ver`, `rmax_ver`, `teta1_ver`, `teta2_ver` | `sorter_rmin`, `sorter_rmax`, `sorter_theta1`, `sorter_theta2` | Vertical sorter region |
| `pot_ver`, `h_ver`, `phi_ver`, `dpotspot` | `sorter_voltage`, `sorter_gap`, `sorter_phase`, `sorter_dv` | Vertical sorter drive |
| `anuz`, `coupnuz`, `anuz_cyclo`, `pot_ver_vrai` | `sorter_nuz`, `nuz_z_cut`, `nuz`, `sorter_voltage_seen` | |
| `numbid`, `ifreb`, `itebid`, `baryc` | `bidim_first`, `bidim_every`, `bidim_start_step`, `centre_mode` | 2D histograms |
| `ax0`, `ax1`, `ay0`, `ay1` | `bidim_ax0`, `bidim_ax1`, `bidim_ay0`, `bidim_ay1` | Histogram scalings |
| `xgrav`, `ygrav`, `ngrav`, `ngravm` | `bary_x`, `bary_y`, `bary_count`, `bary_size` | Orbit-centre barycentre buffer |

## Removed from the state

- Constants: `pi`, `clum`, `rac8` (module constants).
- Scratch arrays, now local variables: `vec`, `vec1`-`vec6`, `veca`, `vecb_arr`, `vecb1`-`vecb4`, `vecd1`, `vecd2`, `vecp1`-`vecp4`, `vecx`, `vecy`, `vecz`, `vecx1`-`vecx4`, `vecy1`-`vecy4`, `vecpc`, `vex1`-`vex4`, `vey1`-`vey4`, `ivec`, `jvec`, `dxtmp`, `rot`, `dbzxlu`, `dbzylu`, inverse grid steps (`hrm1`, `htetm1`, `hrhtet`, `hrhtetm1`, `hx1m1`, ...).
- End-of-run values computed at the last step, now local: `pxf`, `pyf`, `pf`, `tprf`, `wf`, `tpzf`. Those saved one step earlier stay in the state as `end_x`, `end_y`, `end_z`, `end_r`, `end_theta_deg`, `end_phase_deg` (`xf`, `yf`, `zf`, `rf`, `tetf`, `phif`).
- Design values computed and never written: `w0th`, `w0mevth`, `ga0th`, `r0th`, `b0th`, `w1th`, `w1mevth`, `ga1th`, `b1th`.
- Read and never used: `pothf`, `am02c2`, `idebug`, `np` (equal to `n`).
- Never used (Fortran leftovers or unported features): `a1`-`a3`, `s12`-`s32`, `fx`-`fz`, `vecd3`, `vecd4`, `iflag`, `ifirst`, `nbrbid`, `psi`, `dpsi`, `rcor`, `trcor`, `corgri`, `corfen`, `corfeng`, `nfen`, `fendwsw`, `ficinf`, `ficpie`, `ficfingex`, `reponse`, `rdeb1`, `dbdeb1`, `prf`, `rpcmax`, `rbib1`, `rbib2`, `tetabib1`, `tetabib2`, `rbibtest`, `dbbib`, `xrec`, `yrec`, `courec`, `rinflu`, `covg`, `covd`, `coef_ver`, `coef_cos`.
