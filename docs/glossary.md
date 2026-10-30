# Fortran glossary

Python names of the LIONS state and their Fortran counterparts in `fortran/lions_2.0_ep.f90`.
`lions.state.State` groups the data in components: `run`, `files`, `rf`, `field`, `central`, `outer`,
`box`, `bunch`, `extraction`, `traps`, `sorter`, `bidims`, `outputs`. Names below are paths from a
`State` `d` (for example `d.rf.v_right`). Units are those held in memory (SI, angles in rad) unless
stated. Per-particle arrays have length `bunch.n_particles`. Families, deflectors, dipoles, grids
and trap circles are lists of objects (`IonFamily`, `DeflectorElement`, `DipoleElement`,
`GridElement`, `TrapCircle`); the Fortran parallel arrays of fixed size are gone.

## Run and files

| Fortran | Python | Meaning |
|---|---|---|
| `corps` | `run.title` | Deck title |
| `ficin` | `run.deck_file` | Input deck |
| `ficb` | `files.bfield` | Magnetic field map |
| `ficv1` | `files.central_potential` | Central-region potential (binary) |
| `ficv2` | `files.outer_potential` | Outer-region potential (text) |
| `ficz`, `ficxy` | `files.z_output`, `files.xy_output` | Text outputs, `no` = off |
| `ficsig` | `files.sigma` | Initial 6D sigma matrices |
| `espace` | `bunch.phase_space` | `6d` or `3d` |
| `ficcor`, `ficmat`, `ficdc` | `files.correlation`, `files.matrix`, `files.dc_matrix` | Correlation outputs |
| `ficloa`, `ficsto` | `files.bunch_load`, `files.bunch_save` | Bunch files |
| `ficdebgex`, `ficgex` | `extraction.grid_file_prefix`, `extraction.grids[k].file` | Extraction grid crossing files |
| `repswap` | `run.swap_bytes` | Byte swap of the central potential (bool) |
| `npaspt` | `run.steps_per_turn` | Time steps per turn (dt = harm / fhf / npaspt); < 0 tracks backwards |
| `nbdt` | `run.n_steps` | Number of time steps |
| `ifreq` | `run.output_every` | Output period (steps) |
| `ite` | `run.step` | Current step |
| `t0`, `t`, `dt` | `run.t0`, `run.t`, `run.dt` | Time (s) |
| `cdt` | `run.c_dt` | c dt (m) |
| `numtou` | `run.turn` | Turn counter of particle 0 |
| `tetnew`, `tetold`, `tetall` | `run.theta_deg`, `run.theta_prev_deg`, `run.theta_total_deg` | Azimuth of particle 0 (deg), cumulated |
| `r1th`, `ntours` | `run.design_r_final`, `run.design_turns` | Design final radius (m) and turn count |
| `iflcor` | `bunch.correlations` | Correlated initial bunch (bool) |

## RF and magnet

| Fortran | Python | Meaning |
|---|---|---|
| `fhf`, `omhf`, `omr` | `rf.frequency`, `rf.omega`, `rf.omega_rev` | RF frequency (Hz), angular RF and revolution frequencies |
| `harm`, `heps` | `rf.harmonic`, `rf.harmonic_sign` | Harmonic number; -1 if even, +1 if odd |
| `pothfd`, `pothfg` | `rf.v_right`, `rf.v_left` | Dee voltages (V), sign as used in the field |
| `dphid`, `dphig` | `rf.phase_right`, `rf.phase_left` | Dee phase offsets (rad) |
| `teta1`, `yomega` | `rf.dee_half_angle`, `rf.gap_offset` | Dee half-angle (rad), skewed-gap offset (m) |
| `dbsb`, `b_factor` | `field.db_over_b`, `field.b_factor` | Relative field correction; peak field override (T) or None |
| `tbz`, `d2bz`, `dbzx`, `dbzy` | `field.bz`, `field.d2bz`, `field.dbzdx`, `field.dbzdy` | Median-plane maps on (r, theta) |
| `nr`, `ntet`, `rmin`, `rmax`, `hr`, `htet` | `field.nr`, `field.ntheta`, `field.rmin`, `field.rmax`, `field.dr`, `field.dtheta` | Map grid |
| `ncoh`, `coh`, `phihar` | `field.harmonic_n_coef`, `field.harmonic_coef`, `field.harmonic_phase` | First and second field harmonics added to the map |
| `pot1`, `d2pot1` | `central.pot`, `central.d2pot` | Central-region potential and its z curvature |
| `xmin1`, `xmax1`, `ymin1`, `ymax1`, `hx1`, `hy1` | `central.xmin`, `central.xmax`, `central.ymin`, `central.ymax`, `central.dx`, `central.dy` | Central grid (m) |
| `rmax1` | `central.rmax` | Central / outer region boundary (m) |
| `pot2`, `d2pot2` | `outer.pot`, `outer.d2pot` | Outer-region normalised gap potential |
| `xmin2`, `xmax2`, `ymin2`, `ymax2`, `hx2`, `hy2` | `outer.xmin`, `outer.xmax`, `outer.ymin`, `outer.ymax`, `outer.dx`, `outer.dy` | Outer grid, gap frame (m) |
| `zmin`, `zmax` | `box.zmin`, `box.zmax` | Vertical box (m) |

## Particles

| Fortran | Python | Meaning |
|---|---|---|
| `n` | `bunch.n_particles` | Number of particles (`np` merged into it) |
| `nbrfam`, `nbrion`, `index` | `len(bunch.families)`, `bunch.families[k].size`, `bunch.families[k].start` | Families |
| `amat`, `a`, `qn` | `bunch.families[k].atomic_mass`, `bunch.families[k].mass_number`, `bunch.families[k].charge_state` | Per family (u, -, e) |
| `r0`, `teta0`, `z0`, `pr0`, `pz0`, `w0`, `phi0` | `bunch.families[k].r0`, `bunch.families[k].theta0`, `bunch.families[k].z0`, `bunch.families[k].pr0`, `bunch.families[k].pz0`, `bunch.families[k].w0`, `bunch.families[k].phi0_deg` | Initial conditions per family (m, rad, m, rad, rad, J, deg) |
| `w0mev1` | `bunch.w0_mev_per_u` | Initial energy of family 0 (MeV/u) |
| `x`, `y`, `z`, `px`, `py`, `pz` | `bunch.x`, `bunch.y`, `bunch.z`, `bunch.px`, `bunch.py`, `bunch.pz` | Positions (m), momenta (beta gamma) |
| `ga`, `p`, `r`, `teta` | `bunch.gamma`, `bunch.p`, `bunch.r`, `bunch.theta` | Per particle |
| `am0`, `q`, `am0c`, `am0c2` | `bunch.mass`, `bunch.charge`, `bunch.mass_c`, `bunch.rest_energy` | Per particle (kg, C, kg m/s, J) |
| `qdts2m0`, `qdtsm0c` | `bunch.boris_b_coef`, `bunch.boris_e_coef` | q dt / 2m0 and q dt / m0 c |
| `bx`, `by`, `bz`, `ex`, `ey`, `ez` | `bunch.bx`, `bunch.by`, `bunch.bz`, `bunch.ex`, `bunch.ey`, `bunch.ez` | Fields at the particles (T, V/m) |
| `ok` | `bunch.alive` | Not trapped |
| `box`, `centre`, `periph` | `bunch.in_box`, `bunch.in_central`, `bunch.in_outer` | Region masks |
| `extra`, `sortant` | `bunch.in_element`, `bunch.extracted` | Inside the current extraction element; extracted |
| `franchi`, `indgex` | `bunch.crossed_grid`, `bunch.grid_reached` | Grid crossing at this step; last grid crossed |
| `fenlog` | `bunch.in_window` | Diagnostic window |
| `nbrok`, `nbox`, `nbox1`, `nbox2`, `nextra` | `bunch.n_alive`, `bunch.n_in_box`, `bunch.n_central`, `bunch.n_outer`, `bunch.n_in_element` | Counts |

## Phase space and correlations

| Fortran | Python | Meaning |
|---|---|---|
| `dxs` | `bunch.ps` | Phase-space deviations from particle 0: r, pr, z, pz, s, dp |
| `dxd` | `bunch.ps_cartesian` | Same, Cartesian frame |
| `dxsgri`, `rgri`, `tprgri`, `phigri`, `xgri`, `ygri` | `bunch.ps_projected`, `bunch.proj_offset`, `bunch.proj_angle`, `bunch.proj_phase`, `bunch.proj_x`, `bunch.proj_y` | Projection on the radial line through particle 0 |
| `cor`, `cordc`, `t1` | `bunch.sigma`, `bunch.sigma_dc`, `bunch.dc_transform` | 6x6 sigma matrix, with dispersion correction |

## Extraction

| Fortran | Python | Meaning |
|---|---|---|
| `rex` | `extraction.mode` | Deck flag: < 0 enables deflectors, dipoles, grids |
| `tet1ex`, `iteext` | `extraction.cut_bunch`, `extraction.start_step` | Bunch cut at the deflector radius, or whole bunch from a step |
| `tet2ex`, `eex`, `itegex` | `extraction.grids_on`, -, `extraction.grid_start_step` | Grids enabled, from a step |
| `ndf`, `rdf`, `edf` | `len(extraction.deflectors)`, `extraction.deflectors[i].r`, `extraction.deflectors[i].e` | Electrostatic deflectors (m, V/m) |
| `xdfi`, `ydfi`, `tetdfi`, `uxdfi`, `uydfi` | `extraction.deflectors[i].in_x`, `extraction.deflectors[i].in_y`, `extraction.deflectors[i].in_angle`, `extraction.deflectors[i].in_ux`, `extraction.deflectors[i].in_uy` | Entrance line |
| `xdfo`, `ydfo`, `tetdfo`, `uxdfo`, `uydfo` | `extraction.deflectors[i].out_x`, `.out_y`, `.out_angle`, `.out_ux`, `.out_uy` | Exit line |
| `ndm`, `rdm`, `bdm`, `xdm`, `ydm` | `len(extraction.dipoles)`, `extraction.dipoles[i].r`, `extraction.dipoles[i].gradient`, `extraction.dipoles[i].x`, `extraction.dipoles[i].y` | Magnetic dipoles (m, T/m) |
| `tetdmi`, `uxdmi`, `uydmi`, `tetdmo`, ... | `extraction.dipoles[i].in_angle`, `extraction.dipoles[i].in_ux`, `extraction.dipoles[i].in_uy`, `extraction.dipoles[i].out_angle`, `.out_ux`, `.out_uy` | Entrance and exit lines |
| `ngex`, `tetagex` | `len(extraction.grids)`, `extraction.grids[k].angle_deg` | Extraction grids |
| `xgex`, `ygex`, `xgexc`, `ygexc`, `drgex`, `uxgex`, `uygex` | `extraction.grids[k].x`, `extraction.grids[k].y`, `extraction.grids[k].xc`, `extraction.grids[k].yc`, `extraction.grids[k].length`, `extraction.grids[k].ux`, `extraction.grids[k].uy` | Grid lines (m) |

## Traps, vertical sorter, bidims

| Fortran | Python | Meaning |
|---|---|---|
| `nbrpie`, `xpie`, `ypie`, `rpie` | `len(traps.circles)`, `traps.circles[k].x`, `traps.circles[k].y`, `traps.circles[k].r` | Trap circles (m) |
| `itepiem`, `iflpie` | `traps.last_step`, `traps.follow` | Traps active up to a step; 1 = trapped particles follow particle 0 |
| `rmin_ver`, `rmax_ver`, `teta1_ver`, `teta2_ver` | `sorter.rmin`, `sorter.rmax`, `sorter.theta1`, `sorter.theta2` | Vertical sorter region |
| `pot_ver`, `h_ver`, `phi_ver`, `dpotspot` | `sorter.voltage`, `sorter.gap`, `sorter.phase`, `sorter.dv` | Vertical sorter drive |
| `anuz`, `coupnuz`, `anuz_cyclo`, `pot_ver_vrai` | `sorter.nuz`, `sorter.nuz_z_cut`, `outputs.nuz`, `sorter.voltage_seen` | |
| `numbid`, `ifreb`, `itebid`, `baryc` | `bidims.first`, `bidims.every`, `bidims.start_step`, `bidims.centre_mode` | 2D histograms |
| `ax0`, `ax1`, `ay0`, `ay1` | `bidims.ax0`, `bidims.ax1`, `bidims.ay0`, `bidims.ay1` | Histogram scalings |
| `xgrav`, `ygrav`, `ngrav`, `ngravm` | `bidims.bary_x`, `bidims.bary_y`, `bidims.bary_count`, `bidims.bary_size` | Orbit-centre barycentre buffer |

## Removed from the state

- Constants: `pi`, `clum`, `rac8` (module constants).
- Scratch arrays, now local variables: `vec`, `vec1`-`vec6`, `veca`, `vecb_arr`, `vecb1`-`vecb4`, `vecd1`, `vecd2`, `vecp1`-`vecp4`, `vecx`, `vecy`, `vecz`, `vecx1`-`vecx4`, `vecy1`-`vecy4`, `vecpc`, `vex1`-`vex4`, `vey1`-`vey4`, `ivec`, `jvec`, `dxtmp`, `rot`, `dbzxlu`, `dbzylu`, inverse grid steps (`hrm1`, `htetm1`, `hrhtet`, `hrhtetm1`, `hx1m1`, ...).
- End-of-run values computed at the last step, now local: `pxf`, `pyf`, `pf`, `tprf`, `wf`, `tpzf`. Those saved one step earlier stay in the state as `end_x`, `end_y`, `end_z`, `end_r`, `end_theta_deg`, `end_phase_deg` (`xf`, `yf`, `zf`, `rf`, `tetf`, `phif`).
- Design values computed and never written: `w0th`, `w0mevth`, `ga0th`, `r0th`, `b0th`, `w1th`, `w1mevth`, `ga1th`, `b1th`.
- Read and never used: `pothf`, `am02c2`, `idebug`, `np` (equal to `n`).
- Never used (Fortran leftovers or unported features): `a1`-`a3`, `s12`-`s32`, `fx`-`fz`, `vecd3`, `vecd4`, `iflag`, `ifirst`, `nbrbid`, `psi`, `dpsi`, `rcor`, `trcor`, `corgri`, `corfen`, `corfeng`, `nfen`, `fendwsw`, `ficinf`, `ficpie`, `ficfingex`, `reponse`, `rdeb1`, `dbdeb1`, `prf`, `rpcmax`, `rbib1`, `rbib2`, `tetabib1`, `tetabib2`, `rbibtest`, `dbbib`, `xrec`, `yrec`, `courec`, `rinflu`, `covg`, `covd`, `coef_ver`, `coef_cos`.
