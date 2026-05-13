##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:20:51 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_incendi_rhoGTPase.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 12146 atoms
Valist_getStatistics:  Max atom coordinate:  (62.164, 42.94, 61.694)
Valist_getStatistics:  Min atom coordinate:  (-34.318, -32.355, -47.939)
Valist_getStatistics:  Molecule center:  (13.923, 5.2925, 6.8775)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 13.923 5.2925 6.8775
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 13.923 5.2925 6.8775
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 13.923 5.2925 6.8775
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 63.923 55.2925 56.8775
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -36.077 -44.7075 -43.1225
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 53.923 45.2925 46.8775
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -26.077 -34.7075 -33.1225
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 53.923 45.2925 46.8775
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -26.077 -34.7075 -33.1225
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 69.9651
Vpbe_ctor2:  solute dimensions = 99.194 x 77.584 x 112.357
Vpbe_ctor2:  solute charge = -4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (107.558, 86.371, 120.709)
Vclist_setupGrid:  Grid lower corner = (-39.856, -37.893, -53.477)
Vclist_assignAtoms:  Have 3111488 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.085723
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.212055e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.499000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.726050e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 3.063761e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 9.060024e-02
Vprtstp: contraction number = 9.060024e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.137167e-02
Vprtstp: contraction number = 1.255148e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.578204e-03
Vprtstp: contraction number = 1.387838e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.303113e-04
Vprtstp: contraction number = 1.459326e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.459329e-05
Vprtstp: contraction number = 1.502023e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 5.423517e-06
Vprtstp: contraction number = 1.567794e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 9.688348e-07
Vprtstp: contraction number = 1.786359e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 6.491417e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.290795e+00
Vpmg_setPart:  lower corner = (-36.077, -44.7075, -43.1225)
Vpmg_setPart:  upper corner = (63.923, 55.2925, 56.8775)
Vpmg_setPart:  actual minima = (-36.077, -44.7075, -43.1225)
Vpmg_setPart:  actual maxima = (63.923, 55.2925, 56.8775)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.793325958923E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.649000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 69.9651
Vpbe_ctor2:  solute dimensions = 99.194 x 77.584 x 112.357
Vpbe_ctor2:  solute charge = -4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (107.558, 86.371, 120.709)
Vclist_setupGrid:  Grid lower corner = (-39.856, -37.893, -53.477)
Vclist_assignAtoms:  Have 3111488 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -26.077, -34.7075, -33.1225
VPMG::focusFillBound -- New mesh maxs = 53.923, 45.2925, 46.8775
VPMG::focusFillBound -- Old mesh mins = -36.077, -44.7075, -43.1225
VPMG::focusFillBound -- Old mesh maxs = 63.923, 55.2925, 56.8775
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-26.077, -34.7075, -33.1225)
Vpmg_setPart:  upper corner = (53.923, 45.2925, 46.8775)
Vpmg_setPart:  actual minima = (-36.077, -44.7075, -43.1225)
Vpmg_setPart:  actual maxima = (63.923, 55.2925, 56.8775)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-26.077, -34.7075, -33.1225)
VPMG::extEnergy    Disj part upper corner = (53.923, 45.2925, 46.8775)
VPMG::extEnergy    Old lower corner = (-36.077, -44.7075, -43.1225)
VPMG::extEnergy    Old upper corner = (63.923, 55.2925, 56.8775)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 15243.5 kT
VPMG::extEnergy: extDiEnergy = 7631.17 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.148333
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.606426e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.720300e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.605650e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.298439e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.079949e-01
Vprtstp: contraction number = 1.079949e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.371914e-02
Vprtstp: contraction number = 1.270350e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.942845e-03
Vprtstp: contraction number = 1.416157e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.878761e-04
Vprtstp: contraction number = 1.481724e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.301073e-05
Vprtstp: contraction number = 1.494071e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 6.713198e-06
Vprtstp: contraction number = 1.560819e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.277053e-06
Vprtstp: contraction number = 1.902301e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.215093e-07
Vprtstp: contraction number = 2.517588e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.514013e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.302384e+00
Vpmg_setPart:  lower corner = (-26.077, -34.7075, -33.1225)
Vpmg_setPart:  upper corner = (53.923, 45.2925, 46.8775)
Vpmg_setPart:  actual minima = (-26.077, -34.7075, -33.1225)
Vpmg_setPart:  actual maxima = (53.923, 45.2925, 46.8775)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.810517440244E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.770000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.182109e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:20:59 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_vverm_rhoGTPase.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 10115 atoms
Valist_getStatistics:  Max atom coordinate:  (64.884, 43.128, 36.533)
Valist_getStatistics:  Min atom coordinate:  (-39.094, -37.521, -52.749)
Valist_getStatistics:  Molecule center:  (12.895, 2.8035, -8.108)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 12.895 2.8035 -8.108
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 12.895 2.8035 -8.108
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 12.895 2.8035 -8.108
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 62.895 52.8035 41.892
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -37.105 -47.1965 -58.108
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 52.895 42.8035 31.892
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -27.105 -37.1965 -48.108
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 52.895 42.8035 31.892
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -27.105 -37.1965 -48.108
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 60.4119
Vpbe_ctor2:  solute dimensions = 106.744 x 83 x 91.311
Vpbe_ctor2:  solute charge = -10
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (115.054, 91.725, 100.358)
Vclist_setupGrid:  Grid lower corner = (-44.632, -43.059, -58.287)
Vclist_assignAtoms:  Have 2683752 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.968403
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.906431e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.623000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.578020e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.824292e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 9.450350e-02
Vprtstp: contraction number = 9.450350e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.147702e-02
Vprtstp: contraction number = 1.214454e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.488835e-03
Vprtstp: contraction number = 1.297232e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.068052e-04
Vprtstp: contraction number = 1.389040e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.348692e-05
Vprtstp: contraction number = 1.619249e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 7.031655e-06
Vprtstp: contraction number = 2.099821e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.917155e-06
Vprtstp: contraction number = 2.726464e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 5.872902e-07
Vprtstp: contraction number = 3.063342e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.534893e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.407787e+00
Vpmg_setPart:  lower corner = (-37.105, -47.1965, -58.108)
Vpmg_setPart:  upper corner = (62.895, 52.8035, 41.892)
Vpmg_setPart:  actual minima = (-37.105, -47.1965, -58.108)
Vpmg_setPart:  actual maxima = (62.895, 52.8035, 41.892)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.507417543677E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.828000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 60.4119
Vpbe_ctor2:  solute dimensions = 106.744 x 83 x 91.311
Vpbe_ctor2:  solute charge = -10
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (115.054, 91.725, 100.358)
Vclist_setupGrid:  Grid lower corner = (-44.632, -43.059, -58.287)
Vclist_assignAtoms:  Have 2683752 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -27.105, -37.1965, -48.108
VPMG::focusFillBound -- New mesh maxs = 52.895, 42.8035, 31.892
VPMG::focusFillBound -- Old mesh mins = -37.105, -47.1965, -58.108
VPMG::focusFillBound -- Old mesh maxs = 62.895, 52.8035, 41.892
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-27.105, -37.1965, -48.108)
Vpmg_setPart:  upper corner = (52.895, 42.8035, 31.892)
Vpmg_setPart:  actual minima = (-37.105, -47.1965, -58.108)
Vpmg_setPart:  actual maxima = (62.895, 52.8035, 41.892)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-27.105, -37.1965, -48.108)
VPMG::extEnergy    Disj part upper corner = (52.895, 42.8035, 31.892)
VPMG::extEnergy    Old lower corner = (-37.105, -47.1965, -58.108)
VPMG::extEnergy    Old upper corner = (62.895, 52.8035, 41.892)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 10565.9 kT
VPMG::extEnergy: extDiEnergy = 5305.45 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.913913
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.224261e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.658000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.597510e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.347113e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.078789e-01
Vprtstp: contraction number = 1.078789e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.339590e-02
Vprtstp: contraction number = 1.241754e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.796325e-03
Vprtstp: contraction number = 1.340951e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.512418e-04
Vprtstp: contraction number = 1.398643e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.941543e-05
Vprtstp: contraction number = 1.568825e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 8.891684e-06
Vprtstp: contraction number = 2.255889e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.680304e-06
Vprtstp: contraction number = 3.014394e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 9.055303e-07
Vprtstp: contraction number = 3.378461e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.492575e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.346247e+00
Vpmg_setPart:  lower corner = (-27.105, -37.1965, -48.108)
Vpmg_setPart:  upper corner = (52.895, 42.8035, 31.892)
Vpmg_setPart:  actual minima = (-27.105, -37.1965, -48.108)
Vpmg_setPart:  actual maxima = (52.895, 42.8035, 31.892)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.376923456104E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.744000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.226299e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:22:36 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_incindi_actin.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 5838 atoms
Valist_getStatistics:  Max atom coordinate:  (31.415, 26.179, 37.622)
Valist_getStatistics:  Min atom coordinate:  (-27.84, -31.876, -31.62)
Valist_getStatistics:  Molecule center:  (1.7875, -2.8485, 3.001)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 1.7875 -2.8485 3.001
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 1.7875 -2.8485 3.001
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 1.7875 -2.8485 3.001
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 51.7875 47.1515 53.001
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -48.2125 -52.8485 -46.999
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 41.7875 37.1515 43.001
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -38.2125 -42.8485 -36.999
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 41.7875 37.1515 43.001
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -38.2125 -42.8485 -36.999
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 41.9108
Vpbe_ctor2:  solute dimensions = 61.923 x 60.597 x 71.58
Vpbe_ctor2:  solute charge = -12
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (70.331, 69.131, 80.318)
Vclist_setupGrid:  Grid lower corner = (-33.378, -37.414, -37.158)
Vclist_assignAtoms:  Have 3220255 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.534160
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.178771e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.656400e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 5.861210e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.989863e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 9.718967e-02
Vprtstp: contraction number = 9.718967e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.119527e-02
Vprtstp: contraction number = 1.151900e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.356340e-03
Vprtstp: contraction number = 1.211529e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.724418e-04
Vprtstp: contraction number = 1.271376e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.439344e-05
Vprtstp: contraction number = 1.414590e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.192908e-06
Vprtstp: contraction number = 1.718867e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 9.909307e-07
Vprtstp: contraction number = 2.363349e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 6.389394e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.173739e+00
Vpmg_setPart:  lower corner = (-48.2125, -52.8485, -46.999)
Vpmg_setPart:  upper corner = (51.7875, 47.1515, 53.001)
Vpmg_setPart:  actual minima = (-48.2125, -52.8485, -46.999)
Vpmg_setPart:  actual maxima = (51.7875, 47.1515, 53.001)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.776708715117E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.236000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 41.9108
Vpbe_ctor2:  solute dimensions = 61.923 x 60.597 x 71.58
Vpbe_ctor2:  solute charge = -12
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (70.331, 69.131, 80.318)
Vclist_setupGrid:  Grid lower corner = (-33.378, -37.414, -37.158)
Vclist_assignAtoms:  Have 3220255 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -38.2125, -42.8485, -36.999
VPMG::focusFillBound -- New mesh maxs = 41.7875, 37.1515, 43.001
VPMG::focusFillBound -- Old mesh mins = -48.2125, -52.8485, -46.999
VPMG::focusFillBound -- Old mesh maxs = 51.7875, 47.1515, 53.001
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-38.2125, -42.8485, -36.999)
Vpmg_setPart:  upper corner = (41.7875, 37.1515, 43.001)
Vpmg_setPart:  actual minima = (-48.2125, -52.8485, -46.999)
Vpmg_setPart:  actual maxima = (51.7875, 47.1515, 53.001)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-38.2125, -42.8485, -36.999)
VPMG::extEnergy    Disj part upper corner = (41.7875, 37.1515, 43.001)
VPMG::extEnergy    Old lower corner = (-48.2125, -52.8485, -46.999)
VPMG::extEnergy    Old upper corner = (51.7875, 47.1515, 53.001)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 3.0507 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.518924
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.460993e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.521100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 5.929060e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.064262e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.088137e-01
Vprtstp: contraction number = 1.088137e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.292349e-02
Vprtstp: contraction number = 1.187671e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.626906e-03
Vprtstp: contraction number = 1.258876e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.070934e-04
Vprtstp: contraction number = 1.272928e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.041705e-05
Vprtstp: contraction number = 1.468760e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 5.668222e-06
Vprtstp: contraction number = 1.863501e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.916034e-06
Vprtstp: contraction number = 3.380309e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.596814e-07
Vprtstp: contraction number = 3.442952e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.364841e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.135102e+00
Vpmg_setPart:  lower corner = (-38.2125, -42.8485, -36.999)
Vpmg_setPart:  upper corner = (41.7875, 37.1515, 43.001)
Vpmg_setPart:  actual minima = (-38.2125, -42.8485, -36.999)
Vpmg_setPart:  actual maxima = (41.7875, 37.1515, 43.001)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.418341707872E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 8.044000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 1.892756e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:22:42 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_vverm_actin.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 5959 atoms
Valist_getStatistics:  Max atom coordinate:  (30.884, 25.844, 37.181)
Valist_getStatistics:  Min atom coordinate:  (-29.23, -31.794, -32.297)
Valist_getStatistics:  Molecule center:  (0.827, -2.975, 2.442)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 0.827 -2.975 2.442
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 0.827 -2.975 2.442
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 0.827 -2.975 2.442
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 50.827 47.025 52.442
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -49.173 -52.975 -47.558
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 40.827 37.025 42.442
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -39.173 -42.975 -37.558
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 40.827 37.025 42.442
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -39.173 -42.975 -37.558
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 42.1208
Vpbe_ctor2:  solute dimensions = 62.386 x 60.406 x 71.672
Vpbe_ctor2:  solute charge = -8
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (71.19, 68.714, 80.554)
Vclist_setupGrid:  Grid lower corner = (-34.768, -37.332, -37.835)
Vclist_assignAtoms:  Have 3277025 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.545568
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.211961e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.375800e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 5.866760e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.950253e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 9.378283e-02
Vprtstp: contraction number = 9.378283e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.111300e-02
Vprtstp: contraction number = 1.184972e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.429352e-03
Vprtstp: contraction number = 1.286198e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.923651e-04
Vprtstp: contraction number = 1.345821e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.961924e-05
Vprtstp: contraction number = 1.539741e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.804971e-06
Vprtstp: contraction number = 1.622247e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.350607e-06
Vprtstp: contraction number = 2.810853e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 5.149062e-07
Vprtstp: contraction number = 3.812406e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.274194e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.984669e+00
Vpmg_setPart:  lower corner = (-49.173, -52.975, -47.558)
Vpmg_setPart:  upper corner = (50.827, 47.025, 52.442)
Vpmg_setPart:  actual minima = (-49.173, -52.975, -47.558)
Vpmg_setPart:  actual maxima = (50.827, 47.025, 52.442)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.794191394036E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.005000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 42.1208
Vpbe_ctor2:  solute dimensions = 62.386 x 60.406 x 71.672
Vpbe_ctor2:  solute charge = -8
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (71.19, 68.714, 80.554)
Vclist_setupGrid:  Grid lower corner = (-34.768, -37.332, -37.835)
Vclist_assignAtoms:  Have 3277025 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -39.173, -42.975, -37.558
VPMG::focusFillBound -- New mesh maxs = 40.827, 37.025, 42.442
VPMG::focusFillBound -- Old mesh mins = -49.173, -52.975, -47.558
VPMG::focusFillBound -- Old mesh maxs = 50.827, 47.025, 52.442
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-39.173, -42.975, -37.558)
Vpmg_setPart:  upper corner = (40.827, 37.025, 42.442)
Vpmg_setPart:  actual minima = (-49.173, -52.975, -47.558)
Vpmg_setPart:  actual maxima = (50.827, 47.025, 52.442)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-39.173, -42.975, -37.558)
VPMG::extEnergy    Disj part upper corner = (40.827, 37.025, 42.442)
VPMG::extEnergy    Old lower corner = (-49.173, -52.975, -47.558)
VPMG::extEnergy    Old upper corner = (50.827, 47.025, 52.442)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 2.63326 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.548868
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.471824e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.360500e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 5.832850e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.148508e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.080802e-01
Vprtstp: contraction number = 1.080802e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.318944e-02
Vprtstp: contraction number = 1.220339e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.774600e-03
Vprtstp: contraction number = 1.345470e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.547052e-04
Vprtstp: contraction number = 1.435282e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.843148e-05
Vprtstp: contraction number = 1.508862e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 6.484619e-06
Vprtstp: contraction number = 1.687320e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.173567e-06
Vprtstp: contraction number = 3.351881e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 8.085111e-07
Vprtstp: contraction number = 3.719743e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.249407e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.006872e+00
Vpmg_setPart:  lower corner = (-39.173, -42.975, -37.558)
Vpmg_setPart:  upper corner = (40.827, 37.025, 42.442)
Vpmg_setPart:  actual minima = (-39.173, -42.975, -37.558)
Vpmg_setPart:  actual maxima = (40.827, 37.025, 42.442)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.430127421447E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.044000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 1.967173e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:26:15 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_incendi_Dihydropteridine.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 3589 atoms
Valist_getStatistics:  Max atom coordinate:  (28.156, 23.391, 22.771)
Valist_getStatistics:  Min atom coordinate:  (-28.438, -26.272, -30.207)
Valist_getStatistics:  Molecule center:  (-0.141, -1.4405, -3.718)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = -0.141 -1.4405 -3.718
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = -0.141 -1.4405 -3.718
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = -0.141 -1.4405 -3.718
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 49.859 48.5595 46.282
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -50.141 -51.4405 -53.718
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 39.859 38.5595 36.282
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -40.141 -41.4405 -43.718
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 39.859 38.5595 36.282
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -40.141 -41.4405 -43.718
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 35.2085
Vpbe_ctor2:  solute dimensions = 59.19 x 52.46 x 55.335
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (67.67, 60.739, 64.054)
Vclist_setupGrid:  Grid lower corner = (-33.976, -31.81, -35.745)
Vclist_assignAtoms:  Have 2652003 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.342987
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 9.677580e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.644100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.571090e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.833751e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 9.420735e-02
Vprtstp: contraction number = 9.420735e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.120265e-02
Vprtstp: contraction number = 1.189148e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.420888e-03
Vprtstp: contraction number = 1.268350e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.978857e-04
Vprtstp: contraction number = 1.392690e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.088547e-05
Vprtstp: contraction number = 1.560773e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 6.218505e-06
Vprtstp: contraction number = 2.013408e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.727763e-06
Vprtstp: contraction number = 2.778422e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 5.479822e-07
Vprtstp: contraction number = 3.171628e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.395775e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.242841e+00
Vpmg_setPart:  lower corner = (-50.141, -51.4405, -53.718)
Vpmg_setPart:  upper corner = (49.859, 48.5595, 46.282)
Vpmg_setPart:  actual minima = (-50.141, -51.4405, -53.718)
Vpmg_setPart:  actual maxima = (49.859, 48.5595, 46.282)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 5.164650135580E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.778000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 35.2085
Vpbe_ctor2:  solute dimensions = 59.19 x 52.46 x 55.335
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (67.67, 60.739, 64.054)
Vclist_setupGrid:  Grid lower corner = (-33.976, -31.81, -35.745)
Vclist_assignAtoms:  Have 2652003 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -40.141, -41.4405, -43.718
VPMG::focusFillBound -- New mesh maxs = 39.859, 38.5595, 36.282
VPMG::focusFillBound -- Old mesh mins = -50.141, -51.4405, -53.718
VPMG::focusFillBound -- Old mesh maxs = 49.859, 48.5595, 46.282
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-40.141, -41.4405, -43.718)
Vpmg_setPart:  upper corner = (39.859, 38.5595, 36.282)
Vpmg_setPart:  actual minima = (-50.141, -51.4405, -53.718)
Vpmg_setPart:  actual maxima = (49.859, 48.5595, 46.282)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-40.141, -41.4405, -43.718)
VPMG::extEnergy    Disj part upper corner = (39.859, 38.5595, 36.282)
VPMG::extEnergy    Old lower corner = (-50.141, -51.4405, -53.718)
VPMG::extEnergy    Old upper corner = (49.859, 48.5595, 46.282)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.456546 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.341634
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.211286e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.689900e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 7.157280e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.131910e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.104605e-01
Vprtstp: contraction number = 1.104605e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.359869e-02
Vprtstp: contraction number = 1.231090e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.777849e-03
Vprtstp: contraction number = 1.307369e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.443379e-04
Vprtstp: contraction number = 1.374345e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.824961e-05
Vprtstp: contraction number = 1.565439e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 9.911214e-06
Vprtstp: contraction number = 2.591193e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.328647e-06
Vprtstp: contraction number = 3.358466e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 1.126636e-06
Vprtstp: contraction number = 3.384667e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 4.160003e-07
Vprtstp: contraction number = 3.692409e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 8.604246e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 9.448099e+00
Vpmg_setPart:  lower corner = (-40.141, -41.4405, -43.718)
Vpmg_setPart:  upper corner = (39.859, 38.5595, 36.282)
Vpmg_setPart:  actual minima = (-40.141, -41.4405, -43.718)
Vpmg_setPart:  actual maxima = (39.859, 38.5595, 36.282)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.432322855990E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.591000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 2.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.119330e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Oct 10 10:26:21 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path aligned_vverm_Dihydropteridine.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 3487 atoms
Valist_getStatistics:  Max atom coordinate:  (28.934, 23.931, 21.813)
Valist_getStatistics:  Min atom coordinate:  (-24.98, -22.602, -30.141)
Valist_getStatistics:  Molecule center:  (1.977, 0.6645, -4.164)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 1.977 0.6645 -4.164
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 1.977 0.6645 -4.164
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.78125, 0.78125, 0.78125
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.625, 0.625, 0.625
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.8, 0.8, 0.8 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 1.977 0.6645 -4.164
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 51.977 50.6645 45.836
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -48.023 -49.3355 -54.164
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 41.977 40.6645 35.836
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -38.023 -39.3355 -44.164
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 41.977 40.6645 35.836
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -38.023 -39.3355 -44.164
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 32.8619
Vpbe_ctor2:  solute dimensions = 56.587 x 49.047 x 54.379
Vpbe_ctor2:  solute charge = -2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (64.99, 57.609, 63.03)
Vclist_setupGrid:  Grid lower corner = (-30.518, -28.14, -35.679)
Vclist_assignAtoms:  Have 2827175 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.352027
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.013465e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.700300e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.557800e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.870277e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 8.825348e-02
Vprtstp: contraction number = 8.825348e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.084868e-02
Vprtstp: contraction number = 1.229264e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.497893e-03
Vprtstp: contraction number = 1.380714e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.301429e-04
Vprtstp: contraction number = 1.536444e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.935779e-05
Vprtstp: contraction number = 1.710145e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 7.422858e-06
Vprtstp: contraction number = 1.885995e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.530984e-06
Vprtstp: contraction number = 2.062526e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.628044e-07
Vprtstp: contraction number = 2.369747e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.638909e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.476980e+00
Vpmg_setPart:  lower corner = (-48.023, -49.3355, -54.164)
Vpmg_setPart:  upper corner = (51.977, 50.6645, 45.836)
Vpmg_setPart:  actual minima = (-48.023, -49.3355, -54.164)
Vpmg_setPart:  actual maxima = (51.977, 50.6645, 45.836)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 5.026566502307E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.727000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 32.8619
Vpbe_ctor2:  solute dimensions = 56.587 x 49.047 x 54.379
Vpbe_ctor2:  solute charge = -2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (64.99, 57.609, 63.03)
Vclist_setupGrid:  Grid lower corner = (-30.518, -28.14, -35.679)
Vclist_assignAtoms:  Have 2827175 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -38.023, -39.3355, -44.164
VPMG::focusFillBound -- New mesh maxs = 41.977, 40.6645, 35.836
VPMG::focusFillBound -- Old mesh mins = -48.023, -49.3355, -54.164
VPMG::focusFillBound -- Old mesh maxs = 51.977, 50.6645, 45.836
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-38.023, -39.3355, -44.164)
Vpmg_setPart:  upper corner = (41.977, 40.6645, 35.836)
Vpmg_setPart:  actual minima = (-48.023, -49.3355, -54.164)
Vpmg_setPart:  actual maxima = (51.977, 50.6645, 45.836)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-38.023, -39.3355, -44.164)
VPMG::extEnergy    Disj part upper corner = (41.977, 40.6645, 35.836)
VPMG::extEnergy    Old lower corner = (-48.023, -49.3355, -54.164)
VPMG::extEnergy    Old upper corner = (51.977, 50.6645, 45.836)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.204377 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.337865
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.199143e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.586200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.599540e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.153155e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.090768e-01
Vprtstp: contraction number = 1.090768e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.400337e-02
Vprtstp: contraction number = 1.283808e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.030730e-03
Vprtstp: contraction number = 1.450172e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 3.211402e-04
Vprtstp: contraction number = 1.581403e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 5.400006e-05
Vprtstp: contraction number = 1.681510e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 9.452767e-06
Vprtstp: contraction number = 1.750511e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.697386e-06
Vprtstp: contraction number = 1.795650e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.123647e-07
Vprtstp: contraction number = 1.840269e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.563617e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.352822e+00
Vpmg_setPart:  lower corner = (-38.023, -39.3355, -44.164)
Vpmg_setPart:  upper corner = (41.977, 40.6645, 35.836)
Vpmg_setPart:  actual minima = (-38.023, -39.3355, -44.164)
Vpmg_setPart:  actual maxima = (41.977, 40.6645, 35.836)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.227778584172E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 6.502000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.038635e+01
