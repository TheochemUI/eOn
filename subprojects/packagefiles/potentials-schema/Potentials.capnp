# @brief RPC schema for distributed potential evaluations.
#
# This schema defines the binary communication contract between the light client
# and the RPC server components. The design derives from the C-style structures
# in the eOn [1] project (at the v4 writeup).
#
# # References
# [1] eOn Development Team. eOn Documentation. https://eondocs.org.

@0xbd1f89fa17369103;

# @struct ForceInput
# @brief Input configuration for a potential energy evaluation.
# @field lengthUnit Unit expression for positions/box (default "angstrom").
# @field energyUnit Unit expression for energy output (default "eV").
# Unit strings are parsed by rgpot::units::unit_conversion_factor().
# Examples: "angstrom", "bohr", "eV", "hartree", "kJ/mol", "kcal/mol".
struct ForceInput {
  pos        @0 :List(Float64); # @brief Flat array of atomic coordinates [natoms * 3].
  atmnrs     @1 :List(Int32);   # @brief Array of atomic numbers [natoms].
  box        @2 :List(Float64); # @brief Simulation cell vectors [9] (row-major 3x3).
  lengthUnit @3 :Text = "angstrom"; # @brief Unit for positions and box vectors.
  energyUnit @4 :Text = "eV";       # @brief Unit for energy and forces output.
  hasCharge @5 :Bool = false;       # @brief Override NWChemParams.charge for this step.
  charge @6 :Int32 = 0;             # @brief Per-step charge when hasCharge is true.
  hasMultiplicity @7 :Bool = false; # @brief Override NWChemParams.multiplicity for this step.
  multiplicity @8 :Int32 = 1;       # @brief Per-step spin multiplicity when hasMultiplicity is true.
}

# @struct PotentialResult
# @brief Results returned from a potential energy evaluation.
struct PotentialResult {
  energy  @0 :Float64;       # Calculated potential energy.
  forces  @1 :List(Float64); # Flat array of atomic forces [natoms * 3].
  hessian @2 :List(Float64); # Dense Cartesian Hessian [natoms * 3]^2.
  dipole  @3 :List(Float64); # Total dipole vector [3] in atomic units.
  quadrupole @4 :List(Float64); # Traceless quadrupole [6] in atomic units.
  optimizedPos @5 :List(Float64); # Optimized Cartesian coordinates [natoms * 3].
  frequencies @6 :List(Float64); # Harmonic vibrational frequencies [natoms * 3] in cm^-1.
  intensities @7 :List(Float64); # Harmonic IR intensities [natoms * 3] in atomic units.
  stress @8 :List(Float64); # Stress tensor [9] in energy/length^3 units.
  polarizability @9 :List(Float64); # Polarizability / aoresponse vector (backend-defined length).
  gradient @10 :List(Float64); # Flat nuclear gradient [natoms * 3] in energy/length units.
  normalModes @11 :List(Float64); # Dense Cartesian normal-mode matrix [natoms * 3]^2.
  zeroPointEnergy @12 :Float64; # Vibrational zero-point correction in energy units.
  thermalEnergy @13 :Float64; # Vibrational thermal energy correction in energy units.
  thermalEnthalpy @14 :Float64; # Vibrational thermal enthalpy correction in energy units.
  entropy @15 :Float64; # Vibrational total entropy (backend units).
  heatCapacityCv @16 :Float64; # Constant-volume heat capacity (backend units).
  projectedFrequencies @17 :List(Float64); # Frequencies with translations/rotations projected out.
  projectedIntensities @18 :List(Float64); # IR intensities for projected modes.
  # --- cpmdc OpenCPMD extensions (append-only after shared nwchemc surface) ---
  energyComponents @19 :List(Float64); # ener_com scalars etot..eefield (24) when componentsValid.
  componentsValid @20 :Bool = false;
  chargeIntegrals @21 :List(Float64); # chrg: csumg,csumr,csums,csumsabs
  chargeValid @22 :Bool = false;
  multiStateEnergies @23 :List(Float64); # ener_c + ener_d catalog
  multiStateValid @24 :Bool = false;
  mdTrajectoryRow @25 :List(Float64); # MD ENERGY row incl. EKINC when MD ran
  mdTrajectoryValid @26 :Bool = false;
  embedMdPropsSkipped @27 :Bool = true; # true when MD/PROP not harvested this eval (honest skip)
}

# @struct Capabilities
# @brief Backend self-description returned by `<p>_capabilities_result()`.
#
# Loaders (rgpot) negotiate against this message before dispatching work:
# which calculate operations the backend implements, which CommonMethodSpec
# fields its overlay lowers, and which PotentialConfig arms it accepts.
# Backends fill it from their compiled feature inventory, so a stub build
# reports available = false with the same operation surface as the embed.
struct Capabilities {
  backendName    @0 :Text;  # ABI prefix, e.g. "nwchemc" or "cpmdc".
  backendVersion @1 :Text;  # `<p>_version()` string.
  abiVersion     @2 :Int32; # `<p>_abi_version()` value.
  available      @3 :Bool;  # `<p>_available()`: embed shell linked in.
  operations     @4 :List(Operation); # Calculate operations the ABI serves.
  loweredCommonFields @5 :List(Text); # CommonMethodSpec field names the overlay lowers.
  configKinds    @6 :List(Text); # PotentialConfig arms accepted, e.g. "nwchem".
  schemaVersion  @7 :Text;  # Potentials.capnp release the backend compiled against.
  protocolFamily @8 :Text;  # Stable RPC family, e.g. "rgpot.potentials".
  protocolMajor  @9 :UInt16; # Wire-incompatible changes increment this value.
  protocolMinor  @10 :UInt16; # Additive wire-compatible changes increment this value.
  schemaId       @11 :Text; # Stable Cap'n Proto schema identity.
  bridgeAbiMajor @12 :UInt16; # eindir objective ABI major revision.
  bridgeAbiMinor @13 :UInt16; # eindir objective ABI minor revision.
  bridgeLayout   @14 :UInt32; # eindir objective layout revision.
  dlpackMajor    @15 :UInt16; # DLPack callback major revision.
  dlpackMinor    @16 :UInt16; # DLPack callback minor revision.
  bridgeFeatures @17 :UInt64; # eindir bridge feature bitset.

  enum Operation {
    energy         @0;
    forces         @1;
    gradient       @2;
    hessian        @3;
    dipole         @4;
    polarizability @5;
    quadrupole     @6;
    stress         @7;
    optimize       @8;
    frequencies    @9;
  }
}

# @struct NWChemParams
# @brief NWChem-specific knobs (one backend arm inside PotentialConfig / rgpot params).
#
# Not a standalone "rgpot config language" - only used when the active potential
# is NWChem (or configure targets that backend). Same fields in/out via Cap'n Proto.
struct NWChemDirective {
  keyword @0 :Text;       # Directive keyword inside a block, e.g. "convergence".
  args    @1 :List(Text); # Tokenized directive arguments.
}

struct NWChemGenericStanza {
  name       @0 :Text;                   # NWChem block name, e.g. "driver".
  directives @1 :List(NWChemDirective);  # Structured block body.
}

struct NWChemSetDirective {
  key       @0 :Text; # RTDB key, e.g. "dft:avg_fon".
  value     @1 :Text; # Backward-compatible single NWChem input literal.
  valueType @2 :ValueType = auto;
  values    @3 :List(Text); # Preferred structured value list.

  enum ValueType {
    auto    @0; # Infer through NWChem text parser in full decks; direct embed treats as text.
    text    @1; # NWChem "string" set value.
    double  @2;
    integer @3;
    logical @4;
  }
}

struct NWChemDftSmearing {
  sigmaHartree @0 :Float64 = 0.0; # Molecular DFT smear sigma in Hartree.
  mode         @1 :Mode = fixsz;

  enum Mode {
    fixsz   @0;
    nofixsz @1;
  }
}

# Logical DFT integration grid (replaces directive-only multi-token grid cards).
# Prefer quality preset; optional radial/angular tokens for custom grids (text on embed).
struct NWChemDftGridSpec {
  quality @0 :Quality = unspecified;
  radial  @1 :Text = ""; # Optional radial grid token (e.g. lebedev).
  angular @2 :Text = ""; # Optional angular grid token.

  enum Quality {
    unspecified @0; # Do not emit grid from this struct.
    xcoarse     @1;
    coarse      @2;
    medium      @3;
    fine        @4;
    xfine       @5;
    hugeGrid    @6; # Renders as NWChem token "huge"; named to dodge math.h HUGE in C++.
  }
}

struct NWChemDftStanza {
  xc         @0 :Text = "";                  # Exchange-correlation keyword.
  direct     @1 :Bool = false;               # Emit "direct".
  smearing   @2 :NWChemDftSmearing;          # Emit "smear ...".
  # Residual escape for unmodeled DFT cards only (not grid/convergence/odft/diis).
  directives @3 :List(NWChemDirective);
  iterations @4 :Int32 = 0;                  # Emit/promote dft iterations when >0.
  grid       @5 :Text = "";                  # Legacy single-token grid; prefer gridSpec.
  # Convergence thresholds (convergence energy/density/gradient); promote dft:e_conv / d_conv / g_conv.
  energyConv    @6 :Float64 = 0.0;
  densityConv   @7 :Float64 = 0.0;
  gradientConv  @8 :Float64 = 0.0;
  # Open-shell DFT (odft) -> dft:ipol=2 when true.
  odft          @9 :Bool = false;
  diis          @10 :NWChemToggle = unspecified; # Promote dft:diis; nodiis when disabled.
  nfock         @11 :Int32 = 0;              # Promote dft:nfock (DIIS Fock history).
  levelShift    @12 :Float64 = 0.0;          # Promote dft:rlshift + dft:levelshift=true when >0.
  vectorsInput  @13 :Text = "";              # Promote dft:input vectors on embed.
  vectorsOutput @14 :Text = "";              # Promote dft:output vectors on embed.
  gridSpec      @15 :NWChemDftGridSpec;      # Logical multi-token/preset grid (typed; text on embed).
  disp          @16 :NWChemDftDisp;          # Grimme dispersion: emit "disp vdw <n> ..." tokens.
  cdft          @17 :NWChemCdftSpec;         # Constrained DFT lines.
}

# @struct NWChemDftDisp
# @brief Grimme dispersion controls inside the dft block.
struct NWChemDftDisp {
  enabled @0 :Bool = false;
  vdw     @1 :Int32 = 0;      # 1=D1, 2=D2, 3=D3, 4=D3(BJ) -> "disp vdw <n>".
  s6      @2 :Float64 = 0.0;  # 0 => omit token.
  s8      @3 :Float64 = 0.0;
  sr6     @4 :Float64 = 0.0;
  alpha   @5 :Float64 = 0.0;
}

# @struct NWChemRelativisticStanza
# @brief "relativistic" block: ZORA / Douglas-Kroll / X2C controls.
struct NWChemRelativisticStanza {
  method     @0 :Method = none;
  spinOrbit  @1 :Bool = false;    # ZORA: "so" vs "scalar".
  atomicScf  @2 :Bool = false;    # ZORA atomic reference.
  cutoff     @3 :Float64 = 0.0;   # zora:cutoff; 0 => omit.
  dkOrder    @4 :Int32 = 0;       # "douglas-kroll <n>"; 0 => omit.
  dkModifier @5 :Text = "";       # fpp | dkh | dkfull | dk3 | dk3full.
  directives @6 :List(NWChemDirective);
  enum Method {
    none @0;
    zora @1;
    dk   @2;
    x2c  @3;
  }
}

# @struct NWChemCosmoStanza
# @brief COSMO implicit-solvation block.
struct NWChemCosmoStanza {
  dielec       @0 :Float64 = 0.0; # Solvent dielectric; 0 => omit.
  solvent      @1 :Text = "";     # Named solvent preset, e.g. "water".
  rsolv        @2 :Float64 = 0.0;
  radiusScheme @3 :Text = "";     # bondi | cosmo radii source.
  radii        @4 :List(Float64); # Per-atom radii override.
  iscren       @5 :Int32 = 0;
  minbem       @6 :Int32 = 0;
  lineq        @7 :Int32 = 0;
  doCosmoSmd   @8 :Bool = false;  # Gateway to the smd block.
  directives   @9 :List(NWChemDirective);
}

# @struct NWChemSmdStanza
# @brief SMD solvation parameters (rides on cosmo).
struct NWChemSmdStanza {
  solvent    @0 :Text = "";     # smd:solvent name.
  sola       @1 :Float64 = 0.0; # Abraham hydrogen-bond acidity.
  solb       @2 :Float64 = 0.0; # Basicity.
  soln       @3 :Float64 = 0.0; # Refractive index.
  solg       @4 :Float64 = 0.0; # Macroscopic surface tension.
  directives @5 :List(NWChemDirective);
}

# @struct NWChemConstraintsStanza
# @brief "constraints" block for constrained optimization.
struct NWChemConstraint {
  kind   @0 :Kind;
  atoms  @1 :List(UInt32);   # 1-based atom indices.
  value  @2 :Float64 = 0.0;  # Target value or spring constant.
  enum Kind {
    fixAtom @0;
    bond    @1;
    angle   @2;
    torsion @3;
    spring  @4;
  }
}

struct NWChemConstraintsStanza {
  clear       @0 :Bool = false;
  constraints @1 :List(NWChemConstraint);
  directives  @2 :List(NWChemDirective);
}

# @struct NWChemVibStanza
# @brief "vib" block: thermochemistry grid and mass overrides.
struct NWChemIsotope {
  atom @0 :UInt32; # 1-based atom index.
  mass @1 :Float64;
}

struct NWChemVibStanza {
  temperatures @0 :List(Float64); # "temp <n> t1 t2 ..." grid.
  pressure     @1 :Float64 = 0.0; # atm; 0 => omit.
  reuse        @2 :Text = "";     # "reuse <hessian file>".
  animate      @3 :Bool = false;
  masses       @4 :List(NWChemIsotope);
  directives   @5 :List(NWChemDirective);
}

# @struct NWChemBqStanza
# @brief "bq" external point charges.
struct NWChemBqCharge {
  x @0 :Float64;
  y @1 :Float64;
  z @2 :Float64;
  q @3 :Float64;
}

struct NWChemBqStanza {
  units      @0 :Text = "";  # angstroms | au.
  charges    @1 :List(NWChemBqCharge);
  loadFile   @2 :Text = "";  # "load <file> ..." alternative to inline charges.
  forces     @3 :Bool = false;
  directives @4 :List(NWChemDirective);
}

enum NWChemModuleName {
  custom   @0;
  basis    @1;
  bq       @2;
  ccsd     @3;
  cosmo    @4;
  dft      @5;
  dplot    @6;
  drdy     @7;
  driver   @8;
  esp      @9;
  etrans   @10;
  geometry @11;
  gw       @12;
  hessian  @13;
  mcscf    @14;
  md       @15;
  mm       @16;
  mp2      @17;
  ncc      @18;
  nwpw     @19;
  property @20;
  python   @21;
  qmd      @22;
  qmmm     @23;
  rimp2    @24;
  rism     @25;
  scf      @26;
  selci    @27;
  smd      @28;
  tce      @29;
  vib      @30;
  vscf     @31;
  xtb      @32;
  analysis @33;
  argos @34;
  argosDiana @35;
  argosPrep @36;
  argosPrepare @37;
  band @38;
  bandDplot @39;
  brillouinZone @40;
  bsemol @41;
  cckohn @42;
  cellOptimize @43;
  cgsd @44;
  constraints @45;
  cpmd @46;
  cpsd @47;
  ddscf @48;
  diana @49;
  dimpar @50;
  dimqm @51;
  dk @52;
  dmd @53;
  dntmc @54;
  fractionalOccupations @55;
  freeze @56;
  intgrl @57;
  mdXs @58;
  mepgs @59;
  metadynamics @60;
  modelpotential @61;
  neb @62;
  occup @63;
  prepare @64;
  pspFormatter @65;
  pspGenerator @66;
  pspw @67;
  pspwDplot @68;
  pspwQmmm @69;
  pspwWannier @70;
  qmdNamd @71;
  raman @72;
  rel @73;
  rtTddft @74;
  simulationCell @75;
  string @76;
  tamd @77;
  task @78;
  taskShell @79;
  tceMrcc @80;
  tddft @81;
  tddftGradient @82;
  tropt @83;
  vibZone @84;
  waterPseudopotential @85;
  x2c @86;
  zora @87;
}

struct NWChemModuleStanza {
  name       @0 :NWChemModuleName = custom;  # Known NWChem block name.
  customName @1 :Text = "";                  # Block name when name == custom.
  directives @2 :List(NWChemDirective);      # Structured block body.
}

struct NWChemPseudopotentialEntry {
  element     @0 :Text;             # Element symbol, e.g. "Si"; ignored when allElements is true.
  libraryType @1 :LibraryType = library;
  libraryName @2 :Text;             # NWPW library name or file token.
  allElements @3 :Bool = false;     # Use NWChem's "*" default entry for every element.

  enum LibraryType {
    library     @0; # NWChem "library" / PSPW library entry.
    pspwLibrary @1; # Explicit "pspw_library" entry.
    pawLibrary  @2; # PAW library entry.
    cpi         @3; # CPI pseudopotential file entry.
    teter       @4; # Teter pseudopotential file entry.
  }
}

enum NWChemPseudopotentialSpinMode {
  unspecified @0; # Do not emit pspspin.
  enabled     @1; # Promote nwpw:pspspin=true directly; spinRules render text rules.
  disabled    @2; # Render/promote pspspin off.
}

struct NWChemPseudopotentialSpinRule {
  channel                  @0 :Channel = up;
  angularMomentum          @1 :AngularMomentum = s;
  scale                    @2 :Float64 = 1.0;
  ionIndices               @3 :List(Int32);
  hasMagneticQuantumNumber @4 :Bool = false;
  magneticQuantumNumber    @5 :Int32 = 0;

  enum Channel {
    up   @0;
    down @1;
  }

  enum AngularMomentum {
    s @0;
    p @1;
    d @2;
    f @3;
  }
}

struct NWChemPseudopotentialUtermRule {
  angularMomentum @0 :AngularMomentum = s;
  uScale          @1 :Float64 = 0.0;
  jScale          @2 :Float64 = 0.0;
  ionIndices      @3 :List(Int32);

  enum AngularMomentum {
    s @0;
    p @1;
    d @2;
    f @3;
  }
}

enum NWChemPseudopotentialBlockName {
  pseudopotentials        @0; # Emit NWChem "pseudopotentials".
  pseudopotentialLibraries @1; # Emit NWChem "pseudopotential_libraries" alias.
}

struct NWChemPseudopotentialStanza {
  entries    @0 :List(NWChemPseudopotentialEntry);
  directives @1 :List(NWChemDirective); # Extra nwpw directives near the block.
  pspSpin    @2 :NWChemPseudopotentialSpinMode = unspecified; # Direct pspspin state, plus disabled text reset.
  spinRules  @3 :List(NWChemPseudopotentialSpinRule); # Emit/promote pspspin up/down scaling rules.
  semicoreSmall @4 :NWChemToggle = unspecified; # Promote nwpw:psp:semicore_small.
  uterm      @5 :NWChemToggle = unspecified; # Emit/promote NWPW uterm on/off.
  utermRules @6 :List(NWChemPseudopotentialUtermRule); # Emit/promote indexed uterm rules.
  blockName  @7 :NWChemPseudopotentialBlockName = pseudopotentials; # Select the NWChem pseudopotential block spelling.
}

struct NWChemKVector {
  x      @0 :Float64 = 0.0; # First reduced Brillouin-zone coordinate.
  y      @1 :Float64 = 0.0; # Second reduced Brillouin-zone coordinate.
  z      @2 :Float64 = 0.0; # Third reduced Brillouin-zone coordinate.
  weight @3 :Float64 = 0.0; # Optional NWChem kvector weight.
}

struct NWChemBrillouinZoneStanza {
  zoneName          @0 :Text = ""; # NWChem zone_name; empty uses zone_default.
  kVectors          @1 :List(NWChemKVector); # Explicit kvector entries.
  maxKpointsPrint   @2 :Int32 = 0; # max_kpoints_print value.
  monkhorstPackX    @3 :Int32 = 0; # First monkhorst-pack grid value.
  monkhorstPackY    @4 :Int32 = 0; # Second monkhorst-pack grid value.
  monkhorstPackZ    @5 :Int32 = 0; # Third monkhorst-pack grid value.
  directives        @6 :List(NWChemDirective);
  zoneStructureName @7 :Text = ""; # Emit/promote band_structure:zone_name.
  zoneFftName       @8 :Text = ""; # Emit/promote band_fft:zone_name.
  tetrahedronGridX  @9 :Int32 = 0; # First NWPW tetrahedron grid value.
  tetrahedronGridY  @10 :Int32 = 0; # Second NWPW tetrahedron grid value.
  tetrahedronGridZ  @11 :Int32 = 0; # Third NWPW tetrahedron grid value.
  tetrahedronZoneName @12 :Text = ""; # Optional trailing NWPW tetrahedron zone token.
  dosGridX          @13 :Int32 = 0; # First NWPW dos-grid value.
  dosGridY          @14 :Int32 = 0; # Second NWPW dos-grid value.
  dosGridZ          @15 :Int32 = 0; # Third NWPW dos-grid value.
  dosGridZoneName   @16 :Text = ""; # Optional trailing NWPW dos-grid zone token.
  dosFftGridX       @17 :Int32 = 0; # First NWPW dos-fft-grid value.
  dosFftGridY       @18 :Int32 = 0; # Second NWPW dos-fft-grid value.
  dosFftGridZ       @19 :Int32 = 0; # Third NWPW dos-fft-grid value.
  dosFftGridZoneName @20 :Text = ""; # Optional trailing NWPW dos-fft-grid zone token.
}

enum NWChemSimulationCellLatticeKind {
  unspecified @0; # Do not emit a compact cubic lattice directive.
  sc          @1; # Emit/promote NWChem "sc".
  fcc         @2; # Emit/promote NWChem "fcc".
  bcc         @3; # Emit/promote NWChem "bcc".
}

struct NWChemSimulationCellStanza {
  cellName            @0  :Text = ""; # NWChem cell_name; empty uses cell_default.
  boundaryConditions  @1  :Text = ""; # periodic, surface, free-space, ...
  latticeVectorsBohr  @2  :List(Float64); # Flat 3x3 simulation-cell vectors in Bohr.
  ngridX              @3  :Int32 = 0; # First ngrid value.
  ngridY              @4  :Int32 = 0; # Second ngrid value.
  ngridZ              @5  :Int32 = 0; # Third ngrid value.
  ngridSmallX         @6  :Int32 = 0; # First ngrid_small value.
  ngridSmallY         @7  :Int32 = 0; # Second ngrid_small value.
  ngridSmallZ         @8  :Int32 = 0; # Third ngrid_small value.
  boxDeltaBohr        @9  :Float64 = 0.0; # box_delta in Bohr.
  boxOrient           @10 :Bool = false; # Emit/promote box_orient.
  boxDifferentLengths @11 :Bool = false; # Emit/promote box_different_lengths.
  latticeKind         @12 :NWChemSimulationCellLatticeKind = unspecified;
  latticeLengthBohr   @13 :Float64 = 0.0; # sc/fcc/bcc compact lattice length in Bohr.
  directives          @14 :List(NWChemDirective);
}

enum NWChemNwpwBalanceMode {
  unspecified @0; # Do not emit balance/nobalance.
  balance     @1; # Emit/promote NWChem "balance".
  nobalance   @2; # Emit/promote NWChem "nobalance".
}

enum NWChemNwpwBoAlgorithm {
  unspecified    @0; # Do not emit bo_algorithm.
  verlet         @1; # NWChem bo_algorithm verlet.
  velocityVerlet @2; # NWChem bo_algorithm velocity-verlet.
  leapFrog       @3; # NWChem bo_algorithm leap-frog.
}

enum NWChemNwpwSpinMode {
  unspecified  @0; # Do not emit/promote dft/odft.
  restricted   @1; # Emit dft; promote CGSD/BAND/CPSD ispin = 1.
  unrestricted @2; # Emit odft; promote CGSD/BAND/CPSD ispin = 2.
}

enum NWChemNwpwMinimizer {
  unspecified   @0; # Do not emit/promote nwpw:minimizer.
  cgGrassman    @1; # Emit cg grassman; promote nwpw:minimizer = 1.
  cgStiefel     @2; # Emit cg stiefel; promote nwpw:minimizer = 4.
  cgStich       @3; # Emit cg stich; promote nwpw:minimizer = 9.
  lmbfgsGrassman @4; # Emit lmbfgs grassman; promote nwpw:minimizer = 2.
  lmbfgsStiefel @5; # Emit lmbfgs stiefel; promote nwpw:minimizer = 7.
  lmbfgsStich   @6; # Emit lmbfgs stich; promote nwpw:minimizer = 10.
  scfDensity    @7; # Emit scf density; promote nwpw:minimizer = 8.
  scfPotential  @8; # Emit scf potential; promote nwpw:minimizer = 5.
}

enum NWChemNwpwKsAlgorithm {
  unspecified @0; # Do not emit/promote nwpw:ks_algorithm.
  blockCg     @1; # Emit block-cg; promote nwpw:ks_algorithm = -1.
  cg          @2; # Emit cg; promote nwpw:ks_algorithm = 0.
  rmmDiis     @3; # Emit rmm-diis; promote nwpw:ks_algorithm = 1.
}

enum NWChemNwpwScfAlgorithm {
  unspecified @0; # Do not emit/promote nwpw:scf_algorithm.
  simple      @1; # Emit simple; promote nwpw:scf_algorithm = 0.
  broyden     @2; # Emit broyden; promote nwpw:scf_algorithm = 1.
  diis        @3; # Emit diis; promote nwpw:scf_algorithm = 2.
  anderson    @4; # Emit anderson; promote nwpw:scf_algorithm = 3.
}

enum NWChemNwpwToggle {
  unspecified @0; # Do not emit the option.
  enabled     @1; # Emit/promote the option as enabled.
  disabled    @2; # Emit/promote the option as disabled.
}

enum NWChemToggle {
  unspecified @0; # Do not emit the option.
  enabled     @1; # Emit/promote the option as enabled.
  disabled    @2; # Emit/promote the option as disabled.
}

# Closed choices for Gaussian basis library loads (embed + render).
# Open-ended library *names* stay Text on NWChemParams.basis / directives.
enum NWChemBasisLibraryRoot {
  unspecified @0; # Default classic libraries/ resolution (NWCHEM_BASIS_LIBRARY / nwchemrc / compiled).
  classic     @1; # Force basis/libraries/ (not libraries.bse).
  bse         @2; # Force basis/libraries.bse/ (spherical details from BSE files when present).
}

enum NWChemBasisAngularKind {
  unspecified @0; # Default cartesian on embed unless libraries.bse enables read_sphere.
  cartesian   @1; # Force cartesian shells.
  spherical   @2; # Force spherical shells (basis spherical).
}

enum NWChemBasisSegmentMode {
  unspecified @0; # Default segment contractions (NWChem library default).
  segment     @1; # Segment general contractions.
  nosegment   @2; # Do not segment (osegment = false).
}

enum NWChemPrintLevel {
  unspecified @0; # Do not emit a print directive.
  none        @1; # Emit "print none".
  low         @2; # Emit "print low".
  medium      @3; # Emit "print medium".
  high        @4; # Emit "print high".
  debug       @5; # Emit "print debug".
}

enum NWChemTceReference {
  unspecified @0; # Do not emit tce reference.
  dft         @1; # Emit DFT; promote tce:reference = 0.
  hf          @2; # Emit HF; promote tce:reference = 1.
  scf         @3; # Emit SCF; promote tce:reference = 1.
}

enum NWChemTceTwoElectronStorage {
  unspecified @0; # Do not emit/promote tce:model2e.
  default     @1; # Emit/promote tce:model2e = default.
  orbital     @2; # Emit 2eorb; promote tce:model2e = 2eorb.
  spin        @3; # Emit 2espin; promote tce:model2e = 2espin.
}

enum NWChemTceIoAlgorithm {
  unspecified @0; # Do not emit/promote tce:ioalg.
  fortran     @1; # Emit io fortran; promote tce:ioalg = 0.
  eaf         @2; # Emit io eaf; promote tce:ioalg = 1.
  ga          @3; # Emit io ga; promote tce:ioalg = 2.
  sf          @4; # Emit io sf; promote tce:ioalg = 3.
  replicated  @5; # Emit io replicated; promote tce:ioalg = 4.
  dra         @6; # Emit io dra; promote tce:ioalg = 5.
  gaEaf       @7; # Emit io ga_eaf; promote tce:ioalg = 6.
}

enum NWChemTceFreezeMode {
  unspecified @0; # Do not emit a symbolic TCE freeze directive.
  atomic      @1; # Emit "freeze atomic".
  core        @2; # Emit "freeze core".
  coreAtomic  @3; # Emit "freeze core atomic".
}

enum NWChemNwpwSmearType {
  unspecified        @0; # Do not emit fractional_smeartype.
  fixed              @1; # NWChem smear fixed.
  step               @2; # NWChem smear step.
  fermi              @3; # NWChem smear fermi.
  gaussian           @4; # NWChem smear gaussian.
  marzariVanderbilt  @5; # NWChem smear marzari-vanderbilt.
}

enum NWChemNwpwLcaoMode {
  unspecified @0; # Do not emit lcao/lcao_skip.
  lcao        @1; # Emit/promote NWChem "lcao".
  skip        @2; # Emit/promote NWChem "lcao_skip".
}

enum NWChemNwpwEfieldType {
  unspecified @0; # Do not emit/promote nwpw:efield_type.
  periodic    @1; # Emit/promote periodic electric-field mode.
  apc         @2; # Emit/promote APC electric-field mode.
  rgrid       @3; # Emit/promote rgrid electric-field mode.
}

enum NWChemNwpwMappingAlias {
  unspecified @0; # Do not emit a named NWPW mapping alias.
  slab1d      @1; # Emit/promote NWChem "1d-slab" mapping alias.
  hilbert2d   @2; # Emit/promote NWChem "2d-hilbert" mapping alias.
  hcurve2d    @3; # Emit/promote NWChem "2d-hcurve" mapping alias.
}

struct NWChemNwpwStanza {
  energyCutoff       @0 :Float64 = 0.0; # Emit/promote NWPW energy_cutoff.
  wavefunctionCutoff @1 :Float64 = 0.0; # Emit/promote NWPW wavefunction_cutoff.
  ewaldRcut          @2 :Float64 = 0.0; # Emit/promote NWPW ewald_rcut.
  ewaldNcut          @3 :Int32 = 0;     # Emit/promote NWPW ewald_ncut.
  directives         @4 :List(NWChemDirective);
  cellName                    @5  :Text = ""; # Emit/promote NWPW cell_name.
  inputWavefunctionFilename   @6  :Text = ""; # Emit/promote NWPW input_wavefunction_filename.
  outputWavefunctionFilename  @7  :Text = ""; # Emit/promote NWPW output_wavefunction_filename.
  fakeMass                    @8  :Float64 = 0.0; # Emit/promote NWPW fake_mass.
  timeStep                    @9  :Float64 = 0.0; # Emit/promote NWPW time_step.
  loopStart                   @10 :Int32 = 0;     # First value in NWPW loop.
  loopEnd                     @11 :Int32 = 0;     # Second value in NWPW loop.
  toleranceEnergy             @12 :Float64 = 0.0; # First NWPW tolerances value.
  toleranceDensity            @13 :Float64 = 0.0; # Second NWPW tolerances value.
  toleranceGradient           @14 :Float64 = 0.0; # Third NWPW tolerances value.
  exchangeCorrelation         @15 :Text = "";      # Emit/promote NWPW exchange_correlation.
  balanceMode                 @16 :NWChemNwpwBalanceMode = unspecified; # Emit/promote balance/nobalance.
  boStepStart                 @17 :Int32 = 0;       # First NWPW bo_steps value.
  boStepEnd                   @18 :Int32 = 0;       # Second NWPW bo_steps value.
  boTimeStep                  @19 :Float64 = 0.0;   # Emit/promote NWPW bo_time_step.
  boAlgorithm                 @20 :NWChemNwpwBoAlgorithm = unspecified; # Emit/promote NWPW bo_algorithm.
  boFakeMass                  @21 :Float64 = 0.0;   # Emit/promote NWPW bo_fake_mass.
  scalingFirst                @22 :Float64 = 0.0;   # First NWPW scaling value.
  scalingSecond               @23 :Float64 = 0.0;   # Second NWPW scaling value.
  npFftProcesses              @24 :Int32 = 0;       # First NWPW np_dimensions value.
  npOrbitalProcesses          @25 :Int32 = 0;       # Second NWPW np_dimensions value.
  npKspaceProcesses           @26 :Int32 = 0;       # Third NWPW np_dimensions value.
  spinOrbit                   @27 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW spin_orbit.
  parallelIo                  @28 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW parallel_io.
  xyzFilename                 @29 :Text = "";       # Emit/promote NWPW xyz_filename.
  ionMotionFilename           @30 :Text = "";       # Emit/promote NWPW ion_motion_filename.
  electronMotionFilename      @31 :Text = "";       # Emit/promote NWPW emotion_filename.
  hamiltonianMotionFilename   @32 :Text = "";       # Emit/promote NWPW hmotion_filename.
  orbitalMotionFilename       @33 :Text = "";       # Emit/promote NWPW omotion_filename.
  eigenvalueMotionFilename    @34 :Text = "";       # Emit/promote NWPW eigmotion_filename.
  fractionalOrbitalsStart     @35 :Int32 = 0;       # First NWPW fractional_orbitals value.
  fractionalOrbitalsEnd       @36 :Int32 = 0;       # Second NWPW fractional_orbitals value.
  smearTemperature            @37 :Float64 = 0.0;   # Emit/promote NWPW smear temperature RTDB value.
  smearAlpha                  @38 :Float64 = 0.0;   # Emit/promote NWPW smear alpha.
  smearType                   @39 :NWChemNwpwSmearType = unspecified; # Emit/promote NWPW smear type.
  virtualOrbitalsStart        @40 :Int32 = 0;       # First NWPW virtual_orbitals value.
  virtualOrbitalsEnd          @41 :Int32 = 0;       # Second NWPW virtual_orbitals value.
  lcaoMode                    @42 :NWChemNwpwLcaoMode = unspecified; # Emit/promote lcao/lcao_skip.
  ewaldGridX                  @43 :Int32 = 0;       # First NWPW ewald_ngrid value.
  ewaldGridY                  @44 :Int32 = 0;       # Second NWPW ewald_ngrid value.
  ewaldGridZ                  @45 :Int32 = 0;       # Third NWPW ewald_ngrid value.
  noseHoover                  @46 :NWChemNwpwToggle = unspecified; # Emit/promote Nose-Hoover or energy nose toggle.
  noseRestart                 @47 :NWChemNwpwToggle = unspecified; # Emit/promote Nose-Hoover restart/start.
  noseElectronPeriod          @48 :Float64 = 0.0;   # Nose-Hoover electron period Pe.
  noseElectronTemperature     @49 :Float64 = 0.0;   # Nose-Hoover electron temperature Te.
  noseIonPeriod               @50 :Float64 = 0.0;   # Nose-Hoover ion period Pr.
  noseIonTemperature          @51 :Float64 = 0.0;   # Nose-Hoover ion temperature Tr.
  noseElectronChainLength     @52 :Int32 = 0;       # Nose-Hoover Mchain value.
  noseIonChainLength          @53 :Int32 = 0;       # Nose-Hoover Nchain value.
  atomEfield                  @54 :NWChemNwpwToggle = unspecified; # Emit/promote atom_efield.
  atomEfieldGradient          @55 :NWChemNwpwToggle = unspecified; # Emit/promote atom_efield_grad.
  mulliken                    @56 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW mulliken/off.
  mullikenKawai               @57 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW mulliken kawai.
  periodicDipole              @58 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW periodic_dipole.
  electricField               @59 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW efield.
  electricFieldX              @60 :Float64 = 0.0;    # First NWPW efield vector value.
  electricFieldY              @61 :Float64 = 0.0;    # Second NWPW efield vector value.
  electricFieldZ              @62 :Float64 = 0.0;    # Third NWPW efield vector value.
  electricFieldCenterSet      @63 :Bool = false;     # Emit/promote NWPW efield center.
  electricFieldCenterX        @64 :Float64 = 0.0;    # First NWPW efield center value.
  electricFieldCenterY        @65 :Float64 = 0.0;    # Second NWPW efield center value.
  electricFieldCenterZ        @66 :Float64 = 0.0;    # Third NWPW efield center value.
  electricFieldType           @67 :NWChemNwpwEfieldType = unspecified; # Emit/promote NWPW efield type.
  smoothCutoff                @68 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW smooth_cutoff.
  smoothCutoffAfac            @69 :Float64 = 0.0;    # First NWPW smooth_cutoff value.
  smoothCutoffSigma           @70 :Float64 = 0.0;    # Second NWPW smooth_cutoff value.
  cutoffBootWavefunction      @71 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW cutoff_boot_wavefunction.
  fastErf                     @72 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW fast_erf.
  dipoleMotion                @73 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW dipole_motion.
  dipoleMotionFilename        @74 :Text = "";       # Emit/promote NWPW dipole_motion filename.
  rhoUseSymmetry              @75 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW symmetry.
  fmm                         @76 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW fmm.
  fmmLmax                     @77 :Int32 = 0;       # NWPW fmm lmax value.
  fmmLongRange                @78 :Int32 = 0;       # NWPW fmm lr value.
  born                        @79 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW Born solvation.
  bornDielectric              @80 :Float64 = 0.0;   # NWPW born dielectric value.
  bornRelax                   @81 :NWChemNwpwToggle = unspecified; # NWPW born relax flag.
  bornVRadiiAngstrom          @82 :List(Float64);    # NWPW born input radii in Angstrom.
  cpmdProperties              @83 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW cpmd_properties.
  useGridComparison           @84 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW use_grid_cmp.
  director                    @85 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW director.
  directorFilename            @86 :Text = "";       # Emit/promote NWPW director filename.
  cellExpandX                 @87 :Int32 = 0;       # First NWPW expand_cell value.
  cellExpandY                 @88 :Int32 = 0;       # Second NWPW expand_cell value.
  cellExpandZ                 @89 :Int32 = 0;       # Third NWPW expand_cell value.
  mapping                     @90 :Int32 = 0;       # NWPW mapping value.
  rotation                    @91 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW rotation.
  lmaxMultipole               @92 :Int32 = -1;      # NWPW integrate_mult_l / lmax_multipole value.
  fei                         @93 :Bool = false;    # Emit/promote NWPW Fei.
  feiFilename                 @94 :Text = "";       # Filename token in NWPW Fei.
  initialVelocitiesTemperature @95 :Float64 = 0.0;  # NWPW initial_velocities temperature in K.
  initialVelocitiesSeed       @96 :Int32 = -1;      # NWPW initial_velocities seed; negative uses NWChem default.
  makeHmass2                  @97 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW makehmass2.
  translateVectorSet          @98 :Bool = false;    # Emit/promote NWPW translate_vector.
  translateVectorX            @99 :Float64 = 0.0;   # First NWPW translate_vector value.
  translateVectorY            @100 :Float64 = 0.0;  # Second NWPW translate_vector value.
  translateVectorZ            @101 :Float64 = 0.0;  # Third NWPW translate_vector value.
  translateGeometryName       @102 :Text = "";      # Optional NWPW translate_vector geometry name.
  translateReorder            @103 :NWChemNwpwToggle = unspecified; # Optional NWPW translate_vector reorder flag.
  socketType                  @104 :Text = "";      # NWPW socket type token.
  socketIp                    @105 :Text = "";      # Optional NWPW socket IP or socket path token.
  apcSet                      @106 :Bool = false;   # Emit/promote NWPW APC values.
  apcGc                       @107 :Float64 = 0.0;  # NWPW APC Gc value.
  apcGamma                    @108 :List(Float64);  # NWPW APC gamma values.
  translation                 @109 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW translation.
  minimizer                   @110 :NWChemNwpwMinimizer = unspecified; # Emit/promote NWPW minimizer.
  ksAlgorithm                 @111 :NWChemNwpwKsAlgorithm = unspecified; # Emit/promote NWPW KS algorithm.
  scfAlgorithm                @112 :NWChemNwpwScfAlgorithm = unspecified; # Emit/promote NWPW SCF algorithm.
  precondition                @113 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW SCF precondition flag.
  kerkerG0Set                 @114 :Bool = false;   # Emit/promote NWPW SCF kerker g0.
  kerkerG0                    @115 :Float64 = 0.0;  # NWPW SCF kerker g0 value.
  ksAlphaSet                  @116 :Bool = false;   # Emit/promote NWPW SCF alpha.
  ksAlpha                     @117 :Float64 = 0.0;  # NWPW SCF alpha value.
  ksMaxitOrbSet               @118 :Bool = false;   # Emit/promote NWPW SCF inner iterations.
  ksMaxitOrb                  @119 :Int32 = 0;      # NWPW SCF inner iteration count.
  ksMaxitOrbsSet              @120 :Bool = false;   # Emit/promote NWPW SCF outer iterations.
  ksMaxitOrbs                 @121 :Int32 = 0;      # NWPW SCF outer iteration count.
  diisHistoriesSet            @122 :Bool = false;   # Emit/promote NWPW SCF DIIS histories.
  diisHistories               @123 :Int32 = 0;      # NWPW SCF DIIS history count.
  oneElectronGuessSet         @124 :Bool = false;   # Emit/promote NWPW one_electron_guess.
  oneElectronGuessItIn        @125 :Int32 = 0;      # NWPW one-electron guess inner iterations.
  oneElectronGuessItOut       @126 :Int32 = 0;      # NWPW one-electron guess outer iterations.
  oneElectronGuessItOrtho     @127 :Int32 = 0;      # NWPW one-electron guess orthogonalization iterations.
  vfieldFilenames             @128 :List(Text);     # Emit/promote NWPW vfield filenames.
  singlePrecisionHfx          @129 :Bool = false;   # Emit/promote NWPW single_precision_hfx.
  dosAlphaSet                 @130 :Bool = false;   # Emit/promote NWPW DOS alpha.
  dosAlpha                    @131 :Float64 = 0.0;  # NWPW DOS alpha value.
  dosNpointsSet               @132 :Bool = false;   # Emit/promote NWPW DOS point count.
  dosNpoints                  @133 :Int32 = 0;      # NWPW DOS point count.
  dosEminSet                  @134 :Bool = false;   # Emit/promote NWPW DOS minimum energy.
  dosEmin                     @135 :Float64 = 0.0;  # NWPW DOS minimum energy.
  dosEmaxSet                  @136 :Bool = false;   # Emit/promote NWPW DOS maximum energy.
  dosEmax                     @137 :Float64 = 0.0;  # NWPW DOS maximum energy.
  dosFilename                 @138 :Text = "";      # Emit/promote NWPW dos_filename.
  geometryOptimize            @139 :Bool = false;   # Emit/promote NWPW geometry_optimize.
  auxiliaryPotentials         @140 :Bool = false;   # Emit/promote NWPW auxiliary_potentials.
  multiplicity                @141 :Int32 = 0;      # Emit/promote NWPW mult.
  spinMode                    @142 :NWChemNwpwSpinMode = unspecified; # Emit/promote NWPW dft/odft.
  allowTranslation            @143 :Bool = false;   # Emit/promote NWPW allow_translation alias.
  cutoffWavefunction          @144 :Float64 = 0.0;  # Wavefunction value in NWPW cutoff alias.
  cutoffEnergy                @145 :Float64 = 0.0;  # Optional energy value in NWPW cutoff alias.
  mcStepStart                 @146 :Int32 = 0;      # First NWPW mc_steps alias value.
  mcStepEnd                   @147 :Int32 = 0;      # Second NWPW mc_steps alias value.
  etMovecsA                   @148 :Text = "";      # First movecs token in NWPW et.
  etMovecsB                   @149 :Text = "";      # Second movecs token in NWPW et.
  etIonA                      @150 :Text = "";      # First ion token in NWPW et.
  etIonB                      @151 :Text = "";      # Second ion token in NWPW et.
  temperatureIon              @152 :Float64 = 0.0;  # First NWPW temperature argument, Tion.
  temperatureIonPeriod        @153 :Float64 = 0.0;  # Second NWPW temperature argument, Pion.
  temperatureElectron         @154 :Float64 = 0.0;  # Third NWPW temperature argument, Telc.
  temperatureElectronPeriod   @155 :Float64 = 0.0;  # Fourth NWPW temperature argument, Pelc.
  temperatureRestart          @156 :NWChemNwpwToggle = unspecified; # restart/start token in NWPW temperature.
  temperatureIonChainLength   @157 :Int32 = 0;      # First chain argument in NWPW temperature.
  temperatureElectronChainLength @158 :Int32 = 0;   # Second chain argument in NWPW temperature.
  mappingAlias                @159 :NWChemNwpwMappingAlias = unspecified; # Named NWPW mapping alias.
  virtualAliasSet             @160 :Bool = false;   # Emit/promote NWPW virtual alias with defaults.
  virtualAliasStart           @161 :Int32 = 0;      # First NWPW virtual alias value.
  virtualAliasEnd             @162 :Int32 = 0;      # Second NWPW virtual alias value.
  boTimeStepSet               @163 :Bool = false;   # Emit/promote NWPW bo_time_step with the default value when boTimeStep is unset.
  boFakeMassSet               @164 :Bool = false;   # Emit/promote NWPW bo_fake_mass with the default value when boFakeMass is unset.
  scalingSet                  @165 :Bool = false;   # Emit/promote NWPW scaling with default values when scalingFirst is unset.
  npDimensionsSet             @166 :Bool = false;   # Emit/promote NWPW np_dimensions with default values when all dimensions are unset.
  tolerancesSet               @167 :Bool = false;   # Emit/promote NWPW tolerances with default values when all tolerance fields are unset.
  mappingSet                  @168 :Bool = false;   # Emit/promote NWPW mapping with the default value when mapping is unset.
  dosSet                      @169 :Bool = false;   # Emit/promote NWPW dos with the default alpha when DOS scalar fields are unset.
  lcaoMask                    @170 :NWChemNwpwToggle = unspecified; # Emit/promote NWPW lcao_mask.
  lcaoMaskUpOrbitals          @171 :List(Int32); # NWPW lcao_mask up/alpha orbital indexes.
  lcaoMaskDownOrbitals        @172 :List(Int32); # NWPW lcao_mask down/beta orbital indexes.
  occupations                 @173 :List(Float64); # Emit/promote NWPW occupations block occupation values.
  occupationStates            @174 :List(Int32);   # State indexes paired with occupations; unset values default to state 1.
  extraOrbitals               @175 :Int32 = 0;     # Emit/promote NWPW occupations extra_orbitals.
  scalingAtomIndices          @176 :List(Int32);   # One-based atom indexes appended to NWPW scaling.
}

# Logical SCF convergence card (was directive-only "convergence tight|loose|...").
struct NWChemScfConvergence {
  mode @0 :Mode = unspecified;

  enum Mode {
    unspecified @0; # Do not emit convergence line from this struct.
    tight       @1; # Emit "convergence tight".
    loose       @2; # Emit "convergence loose".
    default     @3; # Emit "convergence default" (NWChem spelling).
    energy      @4; # Emit "convergence energy" (energy-only criterion).
    density     @5; # Emit "convergence density".
    gradient    @6; # Emit "convergence gradient".
  }
}

# Logical SCF semidirect integral storage (was directive-only semidirect card).
# filesize/memsize in megawords when >0; promote int2e:* on embed when set.
# enabled alone (no sizes) emits "semidirect" on full and embed (no RTDB for the
# keyword). disabled is a no-op (leave NWChem default; no explicit off keyword).
struct NWChemScfSemidirect {
  enabled  @0 :NWChemToggle = unspecified; # enabled => emit semidirect; disabled => no-op.
  filesize @1 :Int32 = 0; # Promote int2e:filesize when >0 (omit size tokens on embed).
  memsize  @2 :Int32 = 0; # Promote int2e:memsize when >0 (omit size tokens on embed).
}

# @struct NWChemScfStanza
# @brief Typed SCF/HF block controls (vectors, convergence, thresh).
# directives is residual escape for unmodeled SCF cards only.
struct NWChemScfStanza {
  vectorsInput  @0 :Text = "";   # Emit "vectors input <path>"; embed promotes scf:input vectors.
  vectorsOutput @1 :Text = "";   # Emit "vectors output <path>"; embed promotes scf:output vectors.
  maxiter       @2 :Int32 = 0;   # SCF max iterations; embed writes RTDB directly.
  thresh        @3 :Float64 = 0; # SCF convergence threshold; embed writes RTDB directly.
  tol2e         @4 :Float64 = 0; # Two-electron tolerance; embed writes RTDB directly.
  noprint       @5 :Bool = false;# Emit "noprint" (text-only; no single RTDB flag).
  directives    @6 :List(NWChemDirective); # Residual unmodeled SCF cards only.
  # Wavefunction keyword (rhf/uhf/rohf/...). Emitted in the scf block and promoted
  # to RTDB key scf:scftype on embed so HF spin cases are C-ABI accessible without
  # relying only on NWChemParams.scfType / directives.
  wavefunctionType @7 :Text = "";
  # Open-shell count for UHF/ROHF; negative means unset. Promotes scf:nopen.
  nopen            @8 :Int32 = -1;
  diis             @9 :NWChemToggle = unspecified; # Promote scf:diis (enabled/disabled).
  diisBas          @10 :Int32 = 0; # Promote scf:diisbas when >0.
  maxsub           @11 :Int32 = 0; # Promote scf:maxsub when >0.
  lock             @12 :NWChemToggle = unspecified; # Promote scf:lock.
  adapt            @13 :NWChemToggle = unspecified; # Promote scf:adapt.
  noscf            @14 :NWChemToggle = unspecified; # Promote scf:noscf when enabled.
  convergence      @15 :NWChemScfConvergence; # Logical convergence card (not directives).
  semidirect       @16 :NWChemScfSemidirect;  # Logical semidirect / int2e storage.
}

# @struct NWChemCcsdStanza
# @brief Classic CCSD block controls backed by stable ccsd:* RTDB keys.
struct NWChemCcsdStanza {
  maxiter           @0 :Int32 = 0; # Emit/promote ccsd:maxiter.
  thresh            @1 :Float64 = 0.0; # Emit/promote ccsd:thresh.
  tol2e             @2 :Float64 = 0.0; # Emit/promote ccsd:tol2e.
  iprt              @3 :Int32 = 0; # Emit/promote ccsd:iprt.
  maxDiis           @4 :Int32 = 0; # Emit/promote ccsd:maxdiis / diisbas.
  frozenCore        @5 :Int32 = 0; # Emit/promote ccsd:frozen core.
  frozenVirtual     @6 :Int32 = 0; # Emit/promote ccsd:frozen virtual.
  useDisk           @7 :NWChemToggle = unspecified; # disabled emits/promotes nodisk.
  sameSpinScale     @8 :Float64 = 0.0; # Emit/promote ccsd:fss.
  oppositeSpinScale @9 :Float64 = 0.0; # Emit/promote ccsd:fos.
  directives        @10 :List(NWChemDirective);
  printLevel        @11 :NWChemPrintLevel = unspecified; # Emit CCSD print level.
  printItems        @12 :List(Text); # Extra CCSD print names after the level.
  noPrintItems      @13 :List(Text); # CCSD noprint names.
  doa               @14 :List(Int32); # Emit CCSD DOA debugging switches.
  dob               @15 :List(Int32); # Emit CCSD DOB debugging switches.
  dog               @16 :List(Int32); # Emit CCSD DOG debugging switches.
  doh               @17 :List(Int32); # Emit CCSD DOH debugging switches.
  dojk              @18 :List(Int32); # Emit CCSD DOJK debugging switches.
  dos               @19 :List(Int32); # Emit CCSD DOS debugging switches.
  dod               @20 :List(Int32); # Emit CCSD DOD debugging switches.
  useTriplesDriverNonblocking @21 :NWChemToggle = unspecified; # Emit/promote ccsd:use_trpdrv_nb.
  useCcsdOpenmp               @22 :NWChemToggle = unspecified; # Emit/promote ccsd:use_ccsd_omp.
  useTriplesDriverOpenmp      @23 :NWChemToggle = unspecified; # Emit/promote ccsd:use_trpdrv_omp.
  useTriplesDriverOffload     @24 :NWChemToggle = unspecified; # Emit/promote ccsd:use_trpdrv_offload.
}

# @struct NWChemTceStanza
# @brief Tensor Contraction Engine controls for CC, CI, MBPT, EOM-CC, and CR-EOMCC.
struct NWChemTceStanza {
  method                @0  :Text = ""; # Full-deck method keyword, e.g. ccsd(t), cr-ccsd(t), bwccsd.
  model                 @1  :Text = ""; # Emit/promote tce:model normalized value.
  model2e               @2  :NWChemTceTwoElectronStorage = unspecified; # Emit/promote tce:model2e.
  perturbative          @3  :Text = ""; # Promote tce:perturbative.
  ccsdVariant           @4  :Text = ""; # Promote tce:ccsdvar.
  noTriplesSingles      @5  :NWChemToggle = unspecified; # Promote tce:nts.
  reference             @6  :NWChemTceReference = unspecified; # Emit/promote tce:reference.
  frozenCore            @7  :Int32 = 0; # Emit/promote tce:frozen core.
  frozenVirtual         @8  :Int32 = 0; # Emit/promote tce:frozen virtual.
  thresh                @9  :Float64 = 0.0; # Emit/promote tce:thresh.
  levelShift            @10 :Float64 = 0.0; # Emit/promote tce:zlshift.
  leftLevelShift        @11 :Float64 = 0.0; # Emit/promote tce:zlshiftl.
  levelShift2Alpha      @12 :Float64 = 0.0; # First value in tce:zlshift2.
  levelShift2Beta       @13 :Float64 = 0.0; # Second value in tce:zlshift2.
  levelShift3Alpha      @14 :Float64 = 0.0; # First value in tce:zlshift3.
  levelShift3Beta       @15 :Float64 = 0.0; # Second value in tce:zlshift3.
  maxiter               @16 :Int32 = 0; # Emit/promote tce:maxiter.
  ioAlgorithm           @17 :NWChemTceIoAlgorithm = unspecified; # Emit/promote tce:ioalg.
  diis                  @18 :Int32 = 0; # Emit/promote tce:diis.
  diis2                 @19 :Int32 = 0; # Emit/promote tce:diis2.
  diis3                 @20 :Int32 = 0; # Emit/promote tce:diis3.
  eomSolver             @21 :Int32 = 0; # Emit eomsol; promote tce:eoms.
  hbarDimension         @22 :Int32 = 0; # Emit hbard; promote tce:hbard.
  nroots                @23 :Int32 = 0; # Emit/promote tce:nroots.
  target                @24 :Int32 = 0; # Emit/promote tce:target.
  targetSymmetry        @25 :Text = ""; # Emit targetsym; promote tce:targetsym.
  symmetry              @26 :NWChemToggle = unspecified; # Emit/promote tce:symmetry.
  densityMatrix         @27 :NWChemToggle = unspecified; # Emit densmat when enabled; promote tce:densmat.
  densityMatrixFile     @28 :Text = ""; # Emit/promote tce:file_densmat.
  left                  @29 :NWChemToggle = unspecified; # Promote tce:left.
  multipole             @30 :Int32 = 0; # Emit/promote tce:multipole.
  fragment              @31 :Int32 = -1; # Emit/promote tce:fragment when >= 0.
  recomputeFock         @32 :NWChemToggle = unspecified; # Emit fock/nofock; promote tce:recompf.
  activeOccupiedAlpha   @33 :Int32 = 0; # Emit/promote tce:active_oa.
  activeOccupiedBeta    @34 :Int32 = 0; # Emit/promote tce:active_ob.
  activeVirtualAlpha    @35 :Int32 = 0; # Emit/promote tce:active_va.
  activeVirtualBeta     @36 :Int32 = 0; # Emit/promote tce:active_vb.
  activeOccupied        @37 :Int32 = 0; # Emit oact; promote tce:oact.
  activeUnoccupied      @38 :Int32 = 0; # Emit uact; promote tce:uact.
  activeEnergyMin       @39 :Float64 = 0.0; # Emit emin_act; promote tce:eactmin.
  activeEnergyMax       @40 :Float64 = 0.0; # Emit emax_act; promote tce:eactmax.
  activeExcitationLevel @41 :Int32 = 0; # Emit t3a_lvl; promote tce:act_excit_lvl.
  maxDiff               @42 :Float64 = 0.0; # Emit/promote tce:maxdiff.
  atomicTileSize        @43 :Int32 = 0; # Emit attilesize; promote tce:maxs.
  split                 @44 :Int32 = 0; # Emit split; promote tce:ichopx.
  twoElectronMethod     @45 :Int32 = 0; # Emit 2emet; promote tce:i4im.
  diskBackend           @46 :Int32 = -1; # Emit idiskx; promote tce:idiskx when >= 0.
  tileSize              @47 :Int32 = 0; # Emit/promote tce:tilesize.
  cudaDevices           @48 :Int32 = 0; # Emit/promote tce:cuda.
  mrcc                  @49 :NWChemToggle = unspecified; # Promote tce:mrcc.
  tccSpaces             @50 :NWChemToggle = unspecified; # Emit tcc_spaces; promote tce:ltcc.
  eaCcsd                @51 :NWChemToggle = unspecified; # Emit/promote tce:eaccsd when enabled.
  ipCcsd                @52 :NWChemToggle = unspecified; # Emit/promote tce:ipccsd when enabled.
  directives            @53 :List(NWChemDirective);
  freezeMode            @54 :NWChemTceFreezeMode = unspecified; # Emit symbolic freeze directive; not directly promoted.
  dipole                @55 :Bool = false; # Emit/promote TCE dipole response keyword.
  printLevel            @56 :NWChemPrintLevel = unspecified; # Emit TCE print level.
  printItems            @57 :List(Text); # Extra TCE print names after the level.
  quadrupole            @58 :Bool = false; # Emit/promote TCE quadrupole response keyword.
  octupole              @59 :Bool = false; # Emit/promote TCE octupole response keyword.
}

# @struct NWChemMrccDataStanza
# @brief TCE MRCC data block for reference roots, active spaces, and occupation strings.
struct NWChemMrccDataStanza {
  root           @0 :Int32 = 0; # Emit "root <n>".
  casElectrons   @1 :Int32 = 0; # First value in "cas <electrons> <orbitals>".
  casOrbitals    @2 :Int32 = 0; # Second value in "cas <electrons> <orbitals>".
  nref           @3 :Int32 = 0; # Emit explicit nref; zero derives from references length.
  references     @4 :List(Text); # MRCC occupation strings following nref.
  se4t           @5 :Bool = false; # Emit "se4t".
  noAposteriori  @6 :Bool = false; # Emit "no_aposteriori".
  subgroupSize   @7 :Int32 = 0; # Emit "subgroupsize <n>".
  improveTiling  @8 :Bool = false; # Emit "improvetiling".
  usspt          @9 :Bool = false; # Emit "usspt".
  directives     @10 :List(NWChemDirective);
}

# @struct NWChemTaskStanza
# @brief Explicit NWChem "task <theory> <operation>" line.
# Prefer top-level theory/task for embed defaults; use this stanza when emitting
# a full input deck with multiple tasks or non-default theory/operation pairs.
struct NWChemTaskStanza {
  theory    @0 :Text = ""; # scf, dft, mp2, tce, ... (empty => caller omits theory token).
  operation @1 :Text = ""; # energy, gradient, hessian, optimize, property, ...
  ignore    @2 :Bool = false; # Emit "ignore" suffix when true.
}

# @struct NWChemDriverStanza
# @brief Geometry optimization / driver block.
struct NWChemDriverStanza {
  maxiter    @0 :Int32 = 0;     # Driver max steps; embed writes RTDB directly.
  tight      @1 :Bool = false;  # Tight convergence; embed writes RTDB directly.
  loose      @2 :Bool = false;  # Loose convergence; embed writes RTDB directly.
  xyz        @3 :Text = "";     # Emit "xyz <path>" when non-empty.
  directives @4 :List(NWChemDirective);
  gmaxTol    @5 :Float64 = 0;   # Emit "gmax"; embed writes driver:gmax_tol.
  grmsTol    @6 :Float64 = 0;   # Emit "grms"; embed writes driver:grms_tol.
  xmaxTol    @7 :Float64 = 0;   # Emit "xmax"; embed writes driver:xmax_tol.
  xrmsTol    @8 :Float64 = 0;   # Emit "xrms"; embed writes driver:xrms_tol.
}

# @struct NWChemPropertyStanza
# @brief Property evaluation block (dipole, mulliken, ...).
# Bool flags promote to prop:<name> integer 1 on embed (prop_input.F convention).
struct NWChemPropertyStanza {
  dipole     @0 :Bool = false;
  mulliken   @1 :Bool = false;
  quadrupol  @2 :Bool = false; # NWChem keyword "quadrupole" (typo preserved in field name only).
  directives @3 :List(NWChemDirective);
  octupole         @4 :Bool = false; # prop:octupole
  esp              @5 :Bool = false; # prop:esp
  efield           @6 :Bool = false; # prop:efield
  efieldGrad       @7 :Bool = false; # prop:efieldgrad
  electronDensity  @8 :Bool = false; # prop:electrondensity
  spinDensity      @9 :Bool = false; # prop:spindensity
  spinPopulation   @10 :Bool = false; # prop:spinpopulation
  shielding        @11 :Bool = false; # prop:shldopt path via shielding keyword
  hyperfine        @12 :Bool = false; # prop:hypopt
  polarizability   @13 :Bool = false; # prop:polarizability when set in prop_input
}

# @struct NWChemMp2Stanza
# @brief Classic MP2 block (mp2_input.F RTDB keys).
struct NWChemMp2Stanza {
  freezeCore       @0 :Int32 = 0;   # Promote mp2:number frozen core when >0 (also freeze atomic).
  freezeVirtual    @1 :Int32 = 0;   # Promote mp2:number frozen virtual when >0.
  tight            @2 :Bool = false; # Emit tight AO tolerances (aotol2e path).
  aotol2e          @3 :Float64 = 0.0; # Promote mp2:aotol2e when >0.
  aotol2eFock      @4 :Float64 = 0.0; # Promote mp2:aotol2e fock when >0.
  backtol          @5 :Float64 = 0.0; # Promote mp2:backtol when >0.
  sameSpinScale    @6 :Float64 = 0.0; # Promote mp2:fss (SCS-MP2).
  oppositeSpinScale @7 :Float64 = 0.0; # Promote mp2:fos.
  scs              @8 :NWChemToggle = unspecified; # Promote mp2:scs when enabled.
  scratchDisk      @9 :Float64 = 0.0; # Promote mp2:scratchdisk (GB) when >0.
  directives       @10 :List(NWChemDirective);
}

# @struct NWChemTddftStanza
# @brief Linear-response TDDFT block (tddft_input.F RTDB keys).
struct NWChemTddftStanza {
  nroots       @0 :Int32 = 0;    # Promote tddft:nroots.
  tda          @1 :NWChemToggle = unspecified; # Promote tddft:tda (Tamm-Dancoff).
  maxiter      @2 :Int32 = 0;    # Promote tddft:maxiter.
  thresh       @3 :Float64 = 0.0; # Promote tddft:thresh.
  maxvecs      @4 :Int32 = 0;    # Promote tddft:maxvecs.
  singlet      @5 :NWChemToggle = unspecified; # Promote tddft:singlet.
  triplet      @6 :NWChemToggle = unspecified; # Promote tddft:triplet.
  target       @7 :Int32 = 0;    # Promote tddft:target (1-based state index when >0).
  targetSym    @8 :Text = "";    # Promote tddft:targetsym.
  symmetry     @9 :NWChemToggle = unspecified; # Promote tddft:symmetry.
  algorithm    @10 :Int32 = 0;   # Promote tddft:algorithm when >0.
  energyCutoff @11 :Float64 = 0.0; # Promote tddft:ecut + tddft:lecut when >0.
  directives   @12 :List(NWChemDirective);
  # v1.8.0 deepening; keywords from tddft_input.F.
  cdspectrum   @13 :Bool = false;  # tddft:cdspectrum rotatory strengths.
  giao         @14 :Bool = false;  # tddft:giao gauge-including atomic orbitals.
  velocity     @15 :Bool = false;  # tddft:velocity velocity-gauge dipoles.
  civecs       @16 :Bool = false;  # tddft:lcivecs write CI vectors.
  transden     @17 :Bool = false;  # tddft:ltransden transition densities.
  ewinLow      @18 :Float64 = 0.0; # "ewin <low> <high>" energy window (Ha).
  ewinHigh     @19 :Float64 = 0.0;
  gradRoots    @20 :List(Int32);   # "grad ... root <n...> end" excited-state gradients.
}

# @struct NWChemBasisStanza
# @brief Structured Gaussian basis / ECP block (complements top-level basis name).
# Use when callers need spherical/cartesian, segment, library root, or per-element lines.
# Named enums are preferred for closed options; legacy spherical Bool remains for compatibility.
struct NWChemBasisStanza {
  spherical   @0 :Bool = false; # Legacy: true => spherical when angularKind is unspecified.
  segment     @1 :Text = "";    # Optional segment *label* for "* library <segment>" text render.
  ecp         @2 :Text = "";    # Optional ECP library name (embed stores ecp basis when non-empty).
  directives  @3 :List(NWChemDirective); # Per-element "H library 6-31g" etc.
  libraryRoot @4 :NWChemBasisLibraryRoot = unspecified; # classic vs BSE library tree.
  angularKind @5 :NWChemBasisAngularKind = unspecified; # cartesian vs spherical (overrides spherical bool).
  segmentMode @6 :NWChemBasisSegmentMode = unspecified; # segment vs nosegment for bas_tag_lib.
}

# @struct NWChemGeometryStanza
# @brief Geometry block metadata (units/symmetry/noautosym). Coordinates normally
# come from the C ABI positions/atomic_numbers arrays, not this stanza.
struct NWChemGeometryStanza {
  units      @0 :Text = "";    # angstrom, bohr, au, nm, ...
  symmetry   @1 :Text = "";    # c1, d2h, ... or empty.
  noautosym  @2 :Bool = false;
  noautoz    @3 :Bool = false;
  center     @4 :Bool = false; # Emit "center".
  directives @5 :List(NWChemDirective); # Extra geometry directives (not atom lines).
}

struct NWChemInputStanza {
  kind            @0 :Kind = generic;
  generic         @1 :NWChemGenericStanza;
  dft             @2 :NWChemDftStanza;
  set             @3 :NWChemSetDirective;
  raw             @4 :Text;
  module          @5 :NWChemModuleStanza;
  pseudopotential @6 :NWChemPseudopotentialStanza;
  scf             @7 :NWChemScfStanza;
  taskStanza      @8 :NWChemTaskStanza;
  driver          @9 :NWChemDriverStanza;
  property        @10 :NWChemPropertyStanza;
  basisStanza     @11 :NWChemBasisStanza;
  geometry        @12 :NWChemGeometryStanza;
  nwpw            @13 :NWChemNwpwStanza;
  ccsd            @14 :NWChemCcsdStanza;
  tce             @15 :NWChemTceStanza;
  mrccData        @16 :NWChemMrccDataStanza;
  brillouinZone   @17 :NWChemBrillouinZoneStanza;
  simulationCell  @18 :NWChemSimulationCellStanza;
  mp2             @19 :NWChemMp2Stanza;
  tddft           @20 :NWChemTddftStanza;
  relativistic    @21 :NWChemRelativisticStanza;
  cosmo           @22 :NWChemCosmoStanza;
  smd             @23 :NWChemSmdStanza;
  constraints     @24 :NWChemConstraintsStanza;
  vib             @25 :NWChemVibStanza;
  bq              @26 :NWChemBqStanza;
  dplot           @27 :NWChemDplotStanza;
  esp             @28 :NWChemEspStanza;
  qmd             @29 :NWChemQmdStanza;
  raman           @30 :NWChemRamanStanza;
  fon             @31 :NWChemFonStanza;
  neb             @32 :NWChemNebStanza;
  stringMethod    @33 :NWChemStringStanza;
  gw              @34 :NWChemGwStanza;
  etrans          @35 :NWChemEtransStanza;
  rism            @36 :NWChemRismStanza;
  dimQm           @37 :NWChemDimQmStanza;
  metadynamics    @38 :NWChemMetadynamicsStanza;
  cellOptimize    @39 :NWChemCellOptimizeStanza;
  mepgs           @40 :NWChemMepgsStanza;
  tropt           @41 :NWChemTroptStanza;

  enum Kind {
    generic         @0;
    dft             @1;
    set             @2;
    raw             @3;
    module          @4;
    pseudopotential @5;
    scf             @6;
    task            @7; # NWChemInputStanza.taskStanza
    driver          @8;
    property        @9;
    basis           @10; # NWChemInputStanza.basisStanza
    geometry        @11;
    nwpw            @12;
    ccsd            @13;
    tce             @14;
    mrccData        @15;
    brillouinZone   @16;
    simulationCell  @17;
    mp2             @18;
    tddft           @19;
    relativistic    @20;
    cosmo           @21;
    smd             @22;
    constraints     @23;
    vib             @24;
    bq              @25;
    dplot           @26;
    esp             @27;
    qmd             @28;
    raman           @29;
    fon             @30;
    neb             @31;
    stringMethod    @32;
    gw              @33;
    etrans          @34;
    rism            @35;
    dimQm           @36;
    metadynamics    @37;
    cellOptimize    @38;
    mepgs           @39;
    tropt           @40;
  }
}

struct NWChemParams {
  basis        @0 :Text = "sto-3g";  # Gaussian basis (sto-3g, 6-31g*, ...).
  theory       @1 :Text = "scf";     # scf | dft | blyp | b3lyp | ...
  scfType      @2 :Text = "rhf";     # HF: rhf/uhf; DFT: xc (blyp, b3lyp, ...).
  charge       @3 :Int32 = 0;
  multiplicity @4 :Int32 = 1;        # 2S+1.
  enginePath   @5 :Text = "";        # Reserved dynamic engine path; linked embed builds reject non-empty values.
  nwchemRoot   @6 :Text = "";        # NWCHEM_TOP for embed; empty => env.
  task         @7 :Text = "gradient"; # energy | gradient | property; frontend usually calls gradient.
  title        @8 :Text = "";         # Optional NWChem title/start prefix.
  memoryMb     @9 :UInt32 = 0;        # 0 => NWChem defaults / environment.
  scratchDir   @10 :Text = "";        # Optional NWChem scratch directory.
  permanentDir @11 :Text = "";        # Optional NWChem permanent directory.
  inputBlocks  @12 :List(Text);       # Raw NWChem directive blocks applied before task.
  inputStanzas @13 :List(NWChemInputStanza); # Structured NWChem input stanzas.
  # Long-tail / method-specific NWChem options not yet typed above: use
  # NWChemInputStanza.raw, inputBlocks, NWChemSetDirective, or NWChemModuleStanza.custom.
}

# Future backend option structs (extend here, then add a PotentialConfig union arm):
#   struct XTBParams { method @0 :Text = "GFN2-xTB"; ... }
#   struct TBLiteParams { method @0 :Text = "GFN2-xTB"; ... }
#   struct MetatomicParams { modelPath @0 :Text; device @1 :Text = "cpu"; ... }

# --- CPMD structured input (section-oriented, mirrors CPMD &SECTION decks) ---

struct CPMDDirective {
  keyword @0 :Text;       # Line keyword inside a section, e.g. "CONVERGENCE ORBITALS".
  args    @1 :List(Text); # Tokenized arguments on following lines or same line.
}

struct CPMDGenericSection {
  name       @0 :Text;                 # Section name without ampersands, e.g. "CPMD".
  directives @1 :List(CPMDDirective);  # Structured section body.
}

struct CPMDSetDirective {
  key   @0 :Text; # Dotted SECTION.KEYWORD, e.g. "CPMD.PRINT FORCES ON".
  value @1 :Text; # Optional value emitted on the following indented line.
}

struct CPMDKPoint {
  coordinates @0 :List(Float64); # Explicit KPOINTS coordinates kx, ky, kz.
  weight      @1 :Float64 = 1.0; # KPOINTS integration weight.
}

struct CPMDKPointBand {
  points @0 :Int32 = 0;      # Number of interpolated KPOINTS BANDS points.
  start  @1 :List(Float64);  # Initial KPOINTS BANDS vector.
  end    @2 :List(Float64);  # Final KPOINTS BANDS vector.
}

struct CPMDCouplingSurface {
  stateI      @0 :Int32 = 0;      # First coupled KS state index.
  stateJ      @1 :Int32 = 0;      # Second coupled KS state index.
  coefficient @2 :Float64 = 0.0;  # Coupling coefficient for NSURF.
}

struct CPMDSystemSection {
  symmetry        @0 :Int32 = 0;          # SYMMETRY code.
  angstrom        @1 :Bool = true;        # Emit ANGSTROM for CELL / ATOMS.
  cell            @2 :List(Float64);      # CELL vectors or a,b,c,alpha,beta,gamma (6 or 9).
  cutOffRy        @3 :Float64 = 70.0;     # CUTOFF in Rydberg (CPMD default unit).
  scale           @4 :Float64 = 0.0;      # Optional SCALE factor; 0 => omit.
  charge          @5 :Int32 = 0;          # CHARGE of the system.
  multiplicity    @6 :Int32 = 1;          # Spin multiplicity 2S+1 (maps to LSD/spin when >1).
  directives      @7 :List(CPMDDirective);
  densityCutOffRy @8 :Float64 = 0.0;      # DENSITY CUTOFF; 0 => omit.
  poissonSolver   @9 :Text;               # POISSON SOLVER argument, e.g. HOCKNEY.
  poissonParameter @10 :Float64 = 0.0;    # Optional POISSON ... PARAMETER value.
  surface         @11 :Text;              # SURFACE direction, e.g. XY/YZ/ZX.
  referenceCell   @12 :List(Float64);      # REFERENCE CELL values (6) or vectors (9).
  classicalCell   @13 :List(Float64);      # CLASSICAL CELL values (6).
  isotropicCell   @14 :Bool = false;       # ISOTROPIC CELL.
  zFlexibleCell   @15 :Bool = false;       # ZFLEXIBLE CELL.
  densityCutoffNumber @16 :Int32 = 0;      # DENSITY CUTOFF NUMBER; 0 => omit.
  dual            @17 :Float64 = 0.0;      # DUAL density cutoff factor; 0 => omit.
  constantCutoff  @18 :List(Float64);      # CONSTANT CUTOFF akin, skin, eckin (3).
  mesh            @19 :List(Int32);        # MESH nr1, nr2, nr3.
  scaleCartesian  @20 :Bool = false;       # SCALE CARTESIAN; combines with scale as S=.
  doubleGrid      @21 :Text;               # DOUBLE GRID argument, e.g. ON/OFF.
  symmetrizeCoordinates @22 :Bool = false; # SYMMETRIZE COORDINATES.
  tesr            @23 :Int32 = 0;          # TESR value; 0 => omit.
  polymer         @24 :Bool = false;       # POLYMER.
  cluster         @25 :Bool = false;       # CLUSTER.
  cutoffShape     @26 :Text;               # CUTOFF qualifier, e.g. SPHERICAL/NOSPHERICAL.
  hfxCutoff       @27 :List(Float64);       # HFX CUTOFF hfxwfe, hfxdee (2).
  boxWalls        @28 :Float64 = 0.0;       # BOX WALLS skin; 0 => omit.
  nSup            @29 :Int32 = 0;           # NSUP alpha spin states; 0 => omit.
  states          @30 :Int32 = 0;           # STATES count; 0 => omit.
  occupation      @31 :List(Float64);       # OCCUPATION values; length should match states.
  occupationFixed @32 :Bool = false;        # OCCUPATION FIXED.
  externalField   @33 :List(Float64);       # EXTERNAL FIELD vector (3).
  pressure        @34 :Float64 = 0.0;       # PRESSURE; 0 => omit.
  stressTensor    @35 :List(Float64);       # STRESS TENSOR values (9).
  shockVelocity   @36 :Float64 = 0.0;       # SHOCK VELOCITY; 0 => omit.
  checkSymmetryPrecision @37 :Float64 = 0.0; # CHECK SYMMETRY precision; 0 => omit.
  checkSymmetryOff @38 :Bool = false;       # CHECK SYMMETRY OFF.
  wCut            @39 :Float64 = 0.0;       # WCUT; 0 => omit.
  wGauss          @40 :List(Float64);       # WGAUSS sigma values.
  lowSpinExcitation @41 :Text;              # LOW SPIN EXCITATION options, e.g. ROKS.
  lowSpinExcitationLsets @42 :Bool = false; # LOW SPIN EXCITATION LSETS.
  lseParameters   @43 :List(Float64);       # LSE PARAMETERS lsea, lseb (2).
  modifiedGoedecker @44 :Bool = false;      # MODIFIED GOEDECKER.
  modifiedGoedeckerParameters @45 :List(Float64); # MODIFIED GOEDECKER PARAMETERS lambda_ab, lambda_ba (2).
  energyProfile   @46 :Bool = false;        # ENERGY PROFILE.
  pointGroup       @47 :Text;                # POINT GROUP selector line, e.g. AUTO, NAME=C2v, or 26.
  pointGroupDelta  @48 :Float64 = 0.0;       # POINT GROUP DELTA accuracy; 0 => omit.
  pointGroupMolecule @49 :Bool = false;      # POINT GROUP MOLECULE.
  scaleX         @50 :Float64 = 0.0;          # SCALE SX= value; 0 => omit.
  scaleY         @51 :Float64 = 0.0;          # SCALE SY= value; 0 => omit.
  scaleZ         @52 :Float64 = 0.0;          # SCALE SZ= value; 0 => omit.
  kpoints        @53 :List(CPMDKPoint);       # KPOINTS explicit weighted points.
  kpointsScaled  @54 :Bool = false;           # KPOINTS SCALED.
  kpointsOnlyDiagonal @55 :Bool = false;      # KPOINTS ONLYDIAG.
  kpointsMonkhorstPack @56 :List(Int32);      # KPOINTS MONKHORST-PACK nk1, nk2, nk3.
  kpointsMonkhorstSymmetrized @57 :Bool = false; # KPOINTS MONKHORST-PACK SYMMETRIZED.
  kpointsMonkhorstFull @58 :Bool = false;     # KPOINTS MONKHORST-PACK FULL.
  kpointsMonkhorstKdp @59 :Bool = false;      # KPOINTS MONKHORST-PACK KDP.
  kpointsMonkhorstShift @60 :List(Float64);   # KPOINTS MONKHORST-PACK mesh-line SHIFT vector.
  kpointBands    @61 :List(CPMDKPointBand);   # KPOINTS BANDS segments.
  kpointsBlock   @62 :Int32 = 0;              # KPOINTS BLOCK=n; 0 => omit.
  kpointsBlockAll @63 :Bool = false;          # KPOINTS BLOCK ALL.
  kpointsBlockCalculated @64 :Bool = false;   # KPOINTS BLOCK CALCULATED.
  kpointsBlockNoSwap @65 :Bool = false;       # KPOINTS BLOCK NOSWAP.
  lowSpinExcitationPenalty @66 :Float64 = 0.0; # LOW SPIN EXCITATION PENALTY; 0 => omit.
  cdftDonorAtoms @67 :List(Int32);             # DONOR atom indexes.
  cdftDonorWeights @68 :List(Int32);           # DONOR WMULT integer weights.
  cdftAcceptorAtoms @69 :List(Int32);          # ACCEPTOR atom indexes.
  cdftAcceptorHdasDonors @70 :List(Int32);     # ACCEPTOR HDAS donor atom indexes.
  cdftAcceptorWeights @71 :List(Int32);        # ACCEPTOR WMULT integer weights.
  couplingsFiniteDifference @72 :Bool = false; # COUPLINGS FD.
  couplingsFiniteDifferenceDisplacement @73 :Float64 = 0.0; # COUPLINGS FD=eps; 0 => omit value.
  couplingsProductDisplacement @74 :Float64 = 0.0; # COUPLINGS PROD=eps; 0 => omit.
  couplingsLinres @75 :Bool = false;           # COUPLINGS LINRES.
  couplingsLinresTolerance @76 :Float64 = 0.0; # COUPLINGS LINRES TOL=; 0 => omit.
  couplingsLinresNvects @77 :Int32 = 0;        # COUPLINGS LINRES NVECT=; 0 => omit.
  couplingsLinresSpecify @78 :Bool = false;    # COUPLINGS LINRES NVECT=... SPECIFY.
  couplingsLinresBruteForce @79 :Bool = false; # COUPLINGS LINRES BRUTE FORCE.
  couplingsLinresThresholds @80 :List(Float64); # COUPLINGS LINRES THRESHOLDS low/med/high pairs.
  couplingsSurfaces @81 :List(CPMDCouplingSurface); # COUPLINGS NSURF surface triples.
  couplingsFiniteDifferenceAtoms @82 :List(Int32); # COUPLINGS NAT atom indexes.
  cellAbsolute @83 :Bool = false;          # CELL ABSOLUTE.
  cellDegree @84 :Bool = false;            # CELL DEGREE.
  referenceCellAbsolute @85 :Bool = false; # REFERENCE CELL ABSOLUTE.
  referenceCellDegree @86 :Bool = false;   # REFERENCE CELL DEGREE.
  classicalCellAbsolute @87 :Bool = false; # CLASSICAL CELL ABSOLUTE.
  classicalCellDegree @88 :Bool = false;   # CLASSICAL CELL DEGREE.
  cellVectors @89 :Bool = false;            # CELL VECTORS.
  referenceCellVectors @90 :Bool = false;   # REFERENCE CELL VECTORS.
}

struct CPMDCpmdSection {
  optimizeWavefunction @0 :Bool = true;   # OPTIMIZE WAVEFUNCTION.
  molecularDynamics    @1 :Bool = false;  # MOLECULAR DYNAMICS.
  convergenceOrbitals  @2 :Float64 = 1.0e-6;
  maxStep              @3 :Int32 = 0;     # MAXSTEP; 0 => omit (CPMD default).
  timestep             @4 :Float64 = 0.0; # TIMESTEP; 0 => omit.
  restartWavefunction  @5 :Bool = false;  # RESTART WAVEFUNCTION.
  trajectory           @6 :Bool = false;  # TRAJECTORY.
  directives           @7 :List(CPMDDirective);
  optimizeGeometry     @8 :Bool = false;  # OPTIMIZE GEOMETRY.
  maxIter              @9 :Int32 = 0;     # MAXITER; 0 => omit.
  convergenceGeometry  @10 :Float64 = 0.0; # CONVERGENCE GEOMETRY; 0 => omit.
  electronMass         @11 :Float64 = 0.0; # EMASS; 0 => omit.
  molecularDynamicsCp        @12 :Bool = false; # MOLECULAR DYNAMICS CP.
  molecularDynamicsBo        @13 :Bool = false; # MOLECULAR DYNAMICS BO.
  molecularDynamicsEh        @14 :Bool = false; # MOLECULAR DYNAMICS EH.
  molecularDynamicsPt        @15 :Bool = false; # MOLECULAR DYNAMICS PT.
  molecularDynamicsClassical @16 :Bool = false; # MOLECULAR DYNAMICS CLASSICAL.
  molecularDynamicsFile      @17 :Text;         # MOLECULAR DYNAMICS FILE.
  nose                 @18 :Bool = false; # NOSE.
  noseIons             @19 :Bool = false; # NOSE IONS.
  noseElectrons        @20 :Bool = false; # NOSE ELECTRONS.
  berendsen            @21 :Text;         # BERENDSEN.
  langevin             @22 :Bool = false; # LANGEVIN.
  annealing            @23 :Text;         # ANNEALING.
  quench               @24 :Bool = false; # QUENCH.
  rattle               @25 :Bool = false; # RATTLE.
  shake                @26 :Bool = false; # SHAKE.
  constraint           @27 :Text;         # CONSTRAINT.
  trotter              @28 :Text;         # TROTTER.
  restart              @29 :Bool = false; # RESTART.
  printOptions         @30 :Text;         # PRINT.
  storeOptions         @31 :Text;         # STORE.
  centerMoleculeOff    @32 :Bool = false; # CENTER MOLECULE OFF.
  centerMoleculeOn     @33 :Bool = false; # CENTER MOLECULE ON.
  diis                 @34 :Bool = false; # DIIS.
  odiis                @35 :Bool = false; # ODIIS.
  pcg                  @36 :Bool = false; # PCG.
  diagonalization      @37 :Bool = false; # DIAGONALIZATION.
  freeEnergy           @38 :Bool = false; # FREE-ENERGY.
  interface            @39 :Bool = false; # INTERFACE.
  qmmm                 @40 :Bool = false; # QMMM.
  bicanonicalEnsemble  @41 :Bool = false; # BICANONICAL ENSEMBLE.
  cdft                 @42 :Bool = false; # CDFT.
  properties           @43 :Bool = false; # PROPERTIES.
  vdwCorrection        @44 :Text;         # VDW CORRECTION argument, e.g. ON/OFF.
  vdwWannier           @45 :Text;         # VDW WANNIER argument, e.g. ON/OFF.
  dcacp                @46 :Bool = false; # DCACP.
  isolatedMolecule     @47 :Bool = false; # ISOLATED MOLECULE.
  maxRuntime           @48 :Float64 = 0.0; # MAXRUNTIME; 0 => omit.
  timestepElectrons    @49 :Float64 = 0.0; # TIMESTEP ELECTRONS; 0 => omit.
  timestepIons         @50 :Float64 = 0.0; # TIMESTEP IONS; 0 => omit.
  cellMass             @51 :Float64 = 0.0; # CMASS; 0 => omit.
  temperatureElectron  @52 :Float64 = 0.0; # TEMPERATURE ELECTRON; 0 => omit.
  temperature          @53 :Float64 = 0.0; # TEMPERATURE; 0 => omit.
  temperatureRamp      @54 :Bool = false; # TEMPERATURE RAMP.
  temperatureRampTime  @55 :Float64 = 0.0; # TEMPERATURE RAMP second value.
  temperatureRampRate  @56 :Float64 = 0.0; # TEMPERATURE RAMP third value.
  rescaleOldVelocities @57 :Bool = false; # RESCALE OLD VELOCITIES.
  reverseVelocities    @58 :Bool = false; # REVERSE VELOCITIES.
  subtractComVelocity  @59 :Int32 = 0; # SUBTRACT COMVEL; 0 => omit.
  subtractRotVelocity  @60 :Int32 = 0; # SUBTRACT ROTVEL; 0 => omit.
  prngSeed             @61 :Int32 = 0; # PRNGSEED; 0 => omit.
  tempControlIons      @62 :Text; # TEMPCONTROL IONS payload.
  tempControlElectrons @63 :Text; # TEMPCONTROL ELECTRONS payload.
  tempControlCell      @64 :Text; # TEMPCONTROL CELL payload.
  berendsenIons        @65 :Text; # BERENDSEN IONS payload.
  berendsenElectrons   @66 :Text; # BERENDSEN ELECTRONS payload.
  berendsenCell        @67 :Text; # BERENDSEN CELL payload.
  noseIonsThermostat      @68 :Text; # NOSE IONS payload.
  noseElectronsThermostat @69 :Text; # NOSE ELECTRONS payload.
  noseCellThermostat      @70 :Text; # NOSE CELL payload.
  noseParameters          @71 :Text; # NOSE PARAMETERS payload.
  convergenceCell         @72 :Float64 = 0.0; # CONVERGENCE CELL; 0 => omit.
  convergenceAdapt        @73 :Float64 = 0.0; # CONVERGENCE ADAPT; 0 => omit.
  convergenceEnergy       @74 :Float64 = 0.0; # CONVERGENCE ENERGY; 0 => omit.
  convergenceCalfor       @75 :Float64 = 0.0; # CONVERGENCE CALFOR; 0 => omit.
  convergenceRelax        @76 :Int32 = 0; # CONVERGENCE RELAX; 0 => omit.
  convergenceRhofix       @77 :Float64 = 0.0; # CONVERGENCE RHOFIX; 0 => omit.
  convergenceInitial      @78 :Float64 = 0.0; # CONVERGENCE INITIAL; 0 => omit.
  convergenceConstraint   @79 :Text; # CONVERGENCE CONSTRAINT payload.
  annealingIons           @80 :Float64 = 0.0; # ANNEALING IONS; 0 => omit.
  annealingElectrons      @81 :Float64 = 0.0; # ANNEALING ELECTRONS; 0 => omit.
  annealingCell           @82 :Float64 = 0.0; # ANNEALING CELL; 0 => omit.
  dampingIons             @83 :Float64 = 0.0; # DAMPING IONS; 0 => omit.
  dampingElectrons        @84 :Float64 = 0.0; # DAMPING ELECTRONS; 0 => omit.
  dampingCell             @85 :Float64 = 0.0; # DAMPING CELL; 0 => omit.
  hessian                 @86 :Text; # HESSIAN inline options.
  project                 @87 :Text; # PROJECT inline option.
  stressTensorSample      @88 :Int32 = 0; # STRESS TENSOR sample; 0 => omit.
  stressTensorVirial      @89 :Bool = false; # STRESS TENSOR VIRIAL.
  classStressSample       @90 :Int32 = 0; # CLASSTRESS sample; 0 => omit.
  storeSelection          @91 :Text; # STORE inline selection.
  storeInterval           @92 :Int32 = 0; # STORE interval; 0 => omit.
  storeSelfConsistentInterval @93 :Int32 = 0; # STORE SC interval; 0 => omit.
  storeOffSelection       @94 :Text; # STORE OFF inline selection.
  restFileCount           @95 :Int32 = 0; # RESTFILE count; 0 => omit.
  restFileSample          @96 :Text; # RESTFILE SAMPLE payload.
  trajectoryOptions       @97 :Text; # TRAJECTORY inline options.
  trajectorySample        @98 :Int32 = 0; # TRAJECTORY SAMPLE value; 0 => omit.
  trajectoryRange         @99 :Text; # TRAJECTORY RANGE payload.
  movieSample             @100 :Int32 = 0; # MOVIE SAMPLE value; 0 => omit.
  movieOff                @101 :Bool = false; # MOVIE OFF.
  energyBands             @102 :Bool = false; # ENERGYBANDS.
  externalPotential       @103 :Bool = false; # EXTERNAL POTENTIAL.
  externalPotentialAdd    @104 :Bool = false; # EXTERNAL POTENTIAL ADD.
  electrostaticPotential  @105 :Bool = false; # ELECTROSTATIC POTENTIAL.
  electrostaticPotentialSample @106 :Int32 = 0; # ELECTROSTATIC POTENTIAL SAMPLE.
  dipoleDynamicsSample    @107 :Int32 = 0; # DIPOLE DYNAMICS SAMPLE.
  dipoleDynamicsWannier   @108 :Bool = false; # DIPOLE DYNAMICS WANNIER.
  rhoOut                  @109 :Bool = false; # RHOOUT.
  rhoOutSample            @110 :Int32 = 0; # RHOOUT SAMPLE.
  rhoOutBandsCount        @111 :Int32 = 0; # RHOOUT BANDS count.
  rhoOutBands             @112 :Text; # RHOOUT BANDS payload.
  elf                     @113 :Bool = false; # ELF.
  elfParameters           @114 :Text; # ELF PARAMETER payload.
  wannierParameters       @115 :Text; # WANNIER PARAMETER payload.
  wannierOptimization     @116 :Text; # WANNIER OPTIMIZATION option.
  wannierType             @117 :Text; # WANNIER TYPE option.
  wannierReference        @118 :Text; # WANNIER REFERENCE payload.
  wannierSerial           @119 :Bool = false; # WANNIER SERIAL.
  wannierDos              @120 :Bool = false; # WANNIER DOS.
  wannierMolecular        @121 :Bool = false; # WANNIER MOLECULAR.
  wannierWfnOutOptions    @122 :Text; # WANNIER WFNOUT inline options.
  wannierWfnOutPayload    @123 :Text; # WANNIER WFNOUT payload.
  compress                @124 :Text; # COMPRESS inline option.
  memory                  @125 :Text; # MEMORY inline options.
  realSpaceWfnKeep        @126 :Bool = false; # REAL SPACE WFN KEEP.
  realSpaceWfnSize        @127 :Float64 = 0.0; # REAL SPACE WFN SIZE; 0 => omit.
  splineOptions           @128 :Text; # SPLINE inline options.
  splinePoints            @129 :Int32 = 0; # SPLINE POINTS payload; 0 => omit.
  splineRange             @130 :Float64 = 0.0; # SPLINE RANGE payload; 0 => omit.
  finiteDifferences       @131 :Text; # FINITE DIFFERENCES payload.
  taskGroups              @132 :Text; # TASKGROUPS inline option.
  taskGroupsCount         @133 :Int32 = 0; # TASKGROUPS count; 0 => omit.
  distributeFnl           @134 :Text; # DISTRIBUTE FNL ON/OFF.
  filePath                @135 :Text; # FILEPATH payload.
  benchmark               @136 :Text; # BENCHMARK payload.
  mirror                  @137 :Bool = false; # MIRROR.
  shiftPotential          @138 :Text; # SHIFT POTENTIAL payload.
  glocalizationParameters @139 :Text; # GLOCALIZATION PARAMETERS payload.
  glocalizationOptimization @140 :Text; # GLOCALIZATION OPTIMIZATION option.
  gfunctionalType         @141 :Text; # GFUNCTIONAL TYPE option.
  spreadRspace            @142 :Text; # SPREAD RSPACE payload after '='.
  gUnitarityOptions       @143 :Text; # PIPPO unitarity options.
  stepTuning              @144 :Bool = false; # STEP TUNING.
  gAntisym                @145 :Bool = false; # G_ANTISYM.
  gAntisymPenalty         @146 :Bool = false; # G_ANTISYM PENALTY.
  gKick                   @147 :Bool = false; # G_KICK.
  gComplex                @148 :Bool = false; # G_COMPLEX.
  gReal                   @149 :Bool = false; # G_REAL.
  readMatrix              @150 :Bool = false; # READ MATRIX.
  gStepTune               @151 :Bool = false; # G_STEP TUNE.
  glocWfnOutOptions       @152 :Text; # GLOC WFNOUT inline options.
  glocWfnOutPayload       @153 :Text; # GLOC WFNOUT payload line.
  noGeoCheck              @154 :Bool = false; # NO_GEO_CHECK.
  brokenSymmetry          @155 :Bool = false; # BROKEN symmetry.
  distributedLinalg       @156 :Text; # DISTRIBUTED LINALG ON/OFF.
  linalgNewOrtho          @157 :Text; # LINALG NEWORTHO ON/OFF.
  disorthoBlockSize       @158 :Int32 = 0; # DISORTHO_BSIZE; 0 => omit.
  statesBlockSize         @159 :Int32 = 0; # BLOCKSIZE STATES; 0 => omit.
  allToAllPrecision       @160 :Text; # ALLTOALL precision option.
  gshell                  @161 :Bool = false; # GSHELL.
  localPotential          @162 :Bool = false; # LOCAL POTENTIAL.
  cdftOptions             @163 :Text; # CDFT inline options.
  cdftPayload             @164 :Text; # CDFT payload line.
  cdftHdaPayload          @165 :Text; # CDFT HDA payload line.
  vgfactor                @166 :Text; # VGFACTOR payload line.
  vMirror                 @167 :Bool = false; # VMIRROR.
  combineSystemsOptions   @168 :Text; # COMBINE SYSTEMS inline options.
  combineSystemsPayload   @169 :Text; # COMBINE SYSTEMS payload line.
  combineSystemsSabPayload @170 :Text; # COMBINE SYSTEMS SAB payload line.
  kshamOptions            @171 :Text; # KSHAM inline options.
  kshamPayload            @172 :Text; # KSHAM payload line.
  czonesSet               @173 :Text; # CZONES SET payload lines.
  woutOptions             @174 :Text; # WOUT inline options.
  woutPayload             @175 :Text; # WOUT payload line.
  xfmqcTrajectories       @176 :Int32 = 0; # XFMQC trajectory count; 0 => omit.
  molecularDynamicsFileOptions @177 :Text; # MOLECULAR DYNAMICS FILE inline options.
  molecularDynamicsBdTrajectories @178 :Int32 = 0; # MOLECULAR DYNAMICS BD payload.
  parrinelloRahmanOptions @179 :Text; # PARRINELLO-RAHMAN inline options.
  optimizeGeometryOptions @180 :Text; # OPTIMIZE GEOMETRY inline options.
  optimizeGeometrySample  @181 :Int32 = 0; # OPTIMIZE GEOMETRY SAMPLE payload.
  optimizeCombinedOptions @182 :Text; # OPTIMIZE COMBINED inline options.
  optimizeCombinedSample  @183 :Int32 = 0; # OPTIMIZE COMBINED SAMPLE payload.
  cheby                   @184 :Bool = false; # CHEBY.
  cayley                  @185 :Bool = false; # CAYLEY.
  rungeKutta              @186 :Bool = false; # RUNGE-KUTTA.
  forceMatch              @187 :Bool = false; # FORCEMATCH.
  debugOptions            @188 :Text; # DEBUG inline options.
  kohnShamEnergiesOptions @189 :Text; # KOHN-SHAM ENERGIES inline options.
  kohnShamEnergiesCount   @190 :Int32 = 0; # KOHN-SHAM ENERGIES payload count.
  surfaceHoppingOptions   @191 :Text; # SURFACE HOPPING inline options.
  roksOptions             @192 :Text; # ROKS inline options.
  roksExpertPayload       @193 :Text; # ROKS EXPERT payload line.
  pathSampling            @194 :Bool = false; # PATH SAMPLING.
  fixrhoUpwfnOptions      @195 :Text; # FIXRHO UPWFN inline options.
  fixrhoVectors           @196 :Int32 = 0; # FIXRHO VECT payload; 0 => omit.
  fixrhoLoop              @197 :Text; # FIXRHO LOOP min/max payload.
  fixrhoWftol             @198 :Float64 = 0.0; # FIXRHO WFTOL; 0 => omit.
  bogoliubovCorrection    @199 :Text; # BOGOLIUBOV CORRECTION inline option.
  vibrationalAnalysisOptions @200 :Text; # VIBRATIONAL ANALYSIS inline options.
  vibrationalAnalysisSample @201 :Int32 = 0; # VIBRATIONAL ANALYSIS SAMPLE payload.
  vibrationalAnalysisMode @202 :Int32 = 0; # VIBRATIONAL ANALYSIS MODE= value.
  electronicSpectra       @203 :Bool = false; # ELECTRONIC SPECTRA.
  spinOrbitCouplingStates @204 :List(Int32); # SPIN-ORBIT COUPLING state pair.
  propagationSpectra      @205 :Bool = false; # PROPAGATION SPECTRA.
  propagationDistrub      @206 :Bool = false; # PROPAGATION DISTRUB.
  gaugePulse              @207 :Bool = false; # GAUGEPULSE.
  gaugeFieldFrequency     @208 :Float64 = 0.0; # GAUGEFIELD payload; 0 => omit.
  nacv                    @209 :Bool = false; # NACV.
  orbitalHardnessOptions  @210 :Text; # ORBITAL HARDNESS inline options.
  pathIntegral            @211 :Bool = false; # PATH INTEGRAL.
  pathMinimization        @212 :Bool = false; # PATH MINIMIZATION.
  langevinOptions         @213 :Text; # LANGEVIN inline options.
  langevinParameter       @214 :Text; # LANGEVIN payload line.
  qmmmEasy                @215 :Bool = false; # QMMMEASY.
  interfaceOptions        @216 :Text; # INTERFACE inline options.
  trotterFactorCount      @217 :Int32 = 0; # TROTTER FACTOR= count; 0 => omit.
  trotterFactorPayload    @218 :Text; # TROTTER FACTOR payload lines.
  linearResponse          @219 :Bool = false; # LINEAR RESPONSE.
  harmonicReference       @220 :Text; # HARMONIC REFERENCE option.
  scaledMasses            @221 :Text; # SCALED MASSES option.
  tddft                   @222 :Bool = false; # TDDFT.
  ssic                    @223 :Float64 = 0.0; # SSIC payload; 0 => omit.
  nonorthogonalOrbitalsOptions @224 :Text; # NONORTHOGONAL ORBITALS option.
  nonorthogonalOrbitalsLimit @225 :Float64 = 0.0; # NONORTHOGONAL ORBITALS payload.
  lanczosDiagonalizationOptions @226 :Text; # LANCZOS DIAGONALIZATION inline options.
  lanczosParametersCount  @227 :Int32 = 0; # LANCZOS PARAMETERS N= count.
  lanczosParametersPayload @228 :Text; # LANCZOS PARAMETERS payload lines.
  davidsonDiagonalization @229 :Bool = false; # DAVIDSON DIAGONALIZATION.
  davidsonParameters      @230 :Text; # DAVIDSON PARAMETERS payload line.
  alexanderMixing         @231 :Float64 = 0.0; # ALEXANDER MIXING payload; 0 => omit.
  andersonMixingGspace    @232 :Bool = false; # ANDERSON MIXING G-SPACE.
  andersonMixingCount     @233 :Int32 = 0; # ANDERSON MIXING N= count.
  andersonMixingPayload   @234 :Text; # ANDERSON MIXING payload lines.
  broydenMixingOptions    @235 :Text; # BROYDEN MIXING inline options.
  broydenMixingPayload    @236 :Text; # BROYDEN MIXING payload line.
  diisMixingCount         @237 :Int32 = 0; # DIIS MIXING N= count.
  diisMixingPayload       @238 :Text; # DIIS MIXING payload lines.
  moverhoMixing           @239 :Float64 = 0.0; # MOVERHO payload; 0 => omit.
  extrapolateWfnOptions   @240 :Text; # EXTRAPOLATE WFN inline options.
  extrapolateWfnOrder     @241 :Int32 = 0; # EXTRAPOLATE WFN payload; 0 => omit.
  extrapolateConstraintOrder @242 :Int32 = 0; # EXTRAPOLATE CONSTRAINT payload.
  tsdeOptions             @243 :Text; # TSDE inline options.
  tsdpOptions             @244 :Text; # TSDP inline options.
  tcgp                    @245 :Bool = false; # TCGP.
  tsdc                    @246 :Bool = false; # TSDC.
  steepestDescentOptions  @247 :Text; # STEEPEST DESCENT inline options.
  conjugateGradientOptions @248 :Text; # CONJUGATE GRADIENT inline options.
  odiisOptions            @249 :Text; # ODIIS inline options.
  odiisVectors            @250 :Int32 = 0; # ODIIS payload; 0 => omit.
  hamiltonianCutoff       @251 :Float64 = 0.0; # HAMILTONIAN CUTOFF payload.
  gdiisVectors            @252 :Int32 = 0; # GDIIS payload; 0 => omit.
  lbfgsOptions            @253 :Text; # LBFGS inline options.
  lbfgsPayload            @254 :Text; # LBFGS payload line.
  prfoOptions             @255 :Text; # PRFO inline options.
  prfoPayload             @256 :Text; # PRFO payload line.
  hesscore                @257 :Bool = false; # HESSCORE.
  bfgs                    @258 :Bool = false; # BFGS.
  rfoOrder                @259 :Int32 = 0; # RFO ORDER= value; 0 => omit.
  inrParametersCount      @260 :Int32 = 0; # INR PARAMETERS N= count.
  inrParametersPayload    @261 :Text; # INR PARAMETERS payload lines.
  implicitNewtonOptions   @262 :Text; # IMPLICIT NEWTON inline options.
  implicitNewtonMaxIter   @263 :Int32 = 0; # IMPLICIT NEWTON payload; 0 => omit.
  mixsd                   @264 :Float64 = 0.0; # MIXSD payload; 0 => omit.
  mixdiis                 @265 :Float64 = 0.0; # MIXDIIS payload; 0 => omit.
  restartOptions          @266 :Text; # RESTART inline options.
  intFileOptions          @267 :Text; # INTFILE inline options.
  intFileName             @268 :Text; # INTFILE filename payload.
  initializeWavefunctionOptions @269 :Text; # INITIALIZE WAVEFUNCTION inline options.
  rattleParameters        @270 :Text; # RATTLE payload line.
  orthogonalizationOptions @271 :Text; # ORTHOGONALIZATION inline options.
  quenchOptions           @272 :Text; # QUENCH inline options.
  randomizeOptions        @273 :Text; # RANDOMIZE inline options.
  randomizeAmplitude      @274 :Float64 = 0.0; # RANDOMIZE payload; 0 => omit.
  useMts                  @275 :Bool = false; # USE_MTS.
  nabdyZmax               @276 :Int32 = 0; # NABDY_ZMAX payload; 0 => omit.
  nabdySoft               @277 :Float64 = 0.0; # NABDY_SOFT payload; 0 => omit.
  nabdyRedistributeAmplitude @278 :Bool = false; # NABDY_REDISTR_AMPLI.
  nabdyScaleP             @279 :Bool = false; # NABDY_SCALEP.
  nabdyThermo             @280 :Text; # NABDY_THERMO payload line.
  noseIonsOptions         @281 :Text; # NOSE IONS inline options.
  useInStream             @282 :Bool = false; # USE_IN_STREAM.
  useOutStream            @283 :Bool = false; # USE_OUT_STREAM.
  useCublas               @284 :Bool = false; # USE_CUBLAS.
  useCufft                @285 :Bool = false; # USE_CUFFT.
  blasNStreamsPerDevice   @286 :Int32 = 0; # BLAS_N_STREAMS_PER_DEVICE payload.
  blasNDevicesPerTask     @287 :Int32 = 0; # BLAS_N_DEVICES_PER_TASK payload.
  fftNStreamsPerDevice    @288 :Int32 = 0; # FFT_N_STREAMS_PER_DEVICE payload.
  fftNDevicesPerTask      @289 :Int32 = 0; # FFT_N_DEVICES_PER_TASK payload.
  useMpiIo                @290 :Bool = false; # USE_MPI_IO.
  traceOptions            @291 :Text; # TRACE inline options.
  traceProcedure          @292 :Text; # TRACE_PROCEDURE payload line.
  traceMaxDepth           @293 :Int32 = 0; # TRACE_MAX_DEPTH payload.
  traceMaxCalls           @294 :Int32 = 0; # TRACE_MAX_CALLS payload.
  structureOptions        @295 :Text; # STRUCTURE inline options.
  structureSelection      @296 :Text; # STRUCTURE SELECT payload lines.
  wannierNproc            @297 :Int32 = 0; # WANNIER NPROC payload.
  wannierRelocalizeInScf  @298 :Bool = false; # WANNIER RELOCALIZE_IN_SCF.
  wannierRecomputeDipoleMatricesEvery @299 :Int32 = 0; # WANNIER RECOMPUTE_DIPOLE_MATRICES_EVERY payload.
  wannierRelocalizeEvery  @300 :Int32 = 0; # WANNIER RELOCALIZE_EVERY payload.
  paraUseMpiInPlace       @301 :Bool = false; # PARA_USE_MPI_IN_PLACE.
  paraBuffSize            @302 :Int32 = 0; # PARA_BUFF_SIZE payload.
  paraStackBuffSize       @303 :Int32 = 0; # PARA_STACK_BUFF_SIZE payload.
  cpGroups                @304 :Text; # CP_GROUPS payload line.
  fileFusionPayload       @305 :Text; # FILE FUSION payload lines.
  fileMergePayload        @306 :Text; # FILE MERGE payload lines.
  fileSeparationPayload   @307 :Text; # FILE SEPARATION payload lines.
  noReset                 @308 :Int32 = 0; # NO_RESET payload.
  lbfgsNtrust             @309 :Int32 = 0; # LBFGS NTRUST payload; 0 => omit.
  lbfgsNrestt             @310 :Int32 = 0; # LBFGS NRESTT payload; 0 => omit.
  lbfgsNtrstr             @311 :Int32 = 0; # LBFGS NTRSTR payload; 0 => omit.
  lbfgsTrustr             @312 :Float64 = 0.0; # LBFGS TRUSTR payload; 0 => omit.
  prfoMode                @313 :Int32 = 0; # PRFO MODE payload; 0 => omit.
  prfoMdlock              @314 :Int32 = 0; # PRFO MDLOCK payload; 0 => omit.
  prfoTolenv              @315 :Float64 = 0.0; # PRFO TOLENV payload; 0 => omit.
  prfoTrustp              @316 :Float64 = 0.0; # PRFO TRUSTP payload; 0 => omit.
  prfoOmin                @317 :Float64 = 0.0; # PRFO OMIN payload; 0 => omit.
  prfoNsvib               @318 :Int32 = 0; # PRFO NSVIB payload; 0 => omit.
  prfoCoreAtoms           @319 :List(Int32); # PRFO CORE=<n> atom indexes.
  prfoNsmaxp              @320 :Int32 = 0; # PRFO NSMAXP payload; 0 => omit.
  prfoProjectedHessian    @321 :Bool = false; # PRFO PRJHES.
  prfoDisplacement        @322 :Float64 = 0.0; # PRFO DISPL payload; 0 => omit.
  prfoHessianType         @323 :Int32 = 0; # PRFO HESSTYPE payload; 0 => omit.
  lsd                     @324 :Bool = false; # LSD / LOCAL SPIN DENSITY.
  checkMemory             @325 :Bool = false; # CHECK MEMORY.
  trotterFactorizationOff @326 :Bool = false; # TROTTER FACTORIZATION OFF.
  damping                 @327 :Text; # DAMPING payload for ions/electrons/cell.
  dipoleDynamics          @328 :Bool = false; # DIPOLE DYNAMICS.
  realSpaceWfn            @329 :Bool = false; # REAL SPACE WFN.
  noseIonsLocalT0         @330 :Bool = false; # NOSE IONS LOCAL T0.
  noseIonsLocalThermostatCount @331 :Int32 = 0; # NOSE IONS LOCAL thermostat count.
  noseIonsLocalThermostats @332 :Text; # NOSE IONS LOCAL thermostat payload lines.
  noseIonsLocalRangeCount @333 :Int32 = 0; # NOSE IONS LOCAL range count.
  noseIonsLocalRanges     @334 :Text; # NOSE IONS LOCAL range payload lines.
  noseIonsCafesGroupCount @335 :Int32 = 0; # NOSE IONS CAFES group count.
  noseIonsCafesGroups     @336 :Text; # NOSE IONS CAFES group payload lines.
  trotterFactorization    @337 :Float64 = 0.0; # TROTTER FACTORIZATION payload; 0 => omit.
}

struct CPMDDftSection {
  functional    @0 :Text = "BLYP";        # FUNCTIONAL name (LDA, BLYP, PBE, ...).
  lsd           @1 :Bool = false;         # LSD (spin polarized).
  directives    @2 :List(CPMDDirective);
  gcCutoff      @3 :Float64 = 0.0;        # GC-CUTOFF; 0 => omit.
  xcDriver      @4 :Text;                 # XC_DRIVER.
  libxc         @5 :Text;                 # LIBXC.
  lrKernel      @6 :Text;                 # LR KERNEL.
  refunct       @7 :Text;                 # REFUNCT.
  mtsHighFunc   @8 :Text;                 # MTS_HIGH_FUNC.
  mtsLowFunc    @9 :Text;                 # MTS_LOW_FUNC.
  hfx           @10 :Bool = false;        # HFX.
  hfxScreening  @11 :Text;                # HFX-SCREENING.
  hubbard       @12 :Text;                # HUBBARD.
  alpha         @13 :Float64 = 0.0;       # ALPHA; 0 => omit.
  beta          @14 :Float64 = 0.0;       # BETA; 0 => omit.
  oldCode       @15 :Bool = false;        # OLDCODE.
  newCode       @16 :Bool = false;        # NEWCODE.
  correlation   @17 :Text;                # CORRELATION.
  exchange      @18 :Text;                # EXCHANGE.
  becke88       @19 :Bool = false;        # BECKE88.
  hubbardU      @20 :List(CPMDHubbardU);  # Structured HUBBARD U (preferred over hubbard Text).
  hfxWfcCutoff  @21 :Float64 = 0.0;       # HFX WFC cutoff; 0 => omit.
  hfxBlock      @22 :Int32 = 0;           # HFX BLOCK size; 0 => omit.
  hfxDistribution @23 :Text = "";         # HFX distribution scheme token.
}

struct CPMDAtomsPseudopotential {
  element @0 :Text;  # Element symbol, e.g. "H".
  path    @1 :Text;  # Pseudopotential file path or library token.
  lmax    @2 :Int32 = -1; # Optional LMAX; -1 => omit.
  loc     @3 :Int32 = -1; # LOC= local channel; -1 => omit.
  skip    @4 :Int32 = -1; # SKIP= channel; -1 => omit.
  kleinmanBylander @5 :Bool = false; # KLEINMAN-BYLANDER projector form.
  raggio  @6 :Float64 = 0.0; # RAGGIO Gaussian charge width; 0 => omit.
  nonlinearCore @7 :Bool = false; # Per-species NLCC flag.
}

# @struct CPMDAtomConstraint
# @brief One CONSTRAINTS-block entry inside &ATOMS.
struct CPMDAtomConstraint {
  kind   @0 :Text;          # FIX | DIST | BEND | TORSION | COORD | RIGID.
  atoms  @1 :List(Int32);   # 1-based indices, CPMD convention.
  target @2 :Float64 = 0.0; # Constraint value; 0 => current geometry.
  growth @3 :Float64 = 0.0; # GROWTH rate; 0 => omit.
}

struct CPMDIsotope {
  element @0 :Text;
  mass    @1 :Float64;
}

struct CPMDAtomVelocity {
  atom     @0 :Int32;         # 1-based index.
  velocity @1 :List(Float64); # 3 Cartesian components.
}

struct CPMDDummyAtom {
  type    @0 :Int32 = 1;     # DUMMY ATOMS type 1-4.
  atoms   @1 :List(Int32);
  weights @2 :List(Float64); # Weighted dummies (type 4).
}

struct CPMDAtomsSection {
  pseudopotentials @0 :List(CPMDAtomsPseudopotential);
  # Explicit Cartesian coordinates are supplied per step via ForceInput; this
  # section only carries PP and fixed &ATOMS directives for deck rendering.
  directives       @1 :List(CPMDDirective);
  constraints      @2 :List(CPMDAtomConstraint); # CONSTRAINTS block.
  isotopes         @3 :List(CPMDIsotope);        # ISOTOPE per-species masses.
  velocities       @4 :List(CPMDAtomVelocity);   # VELOCITIES block.
  dummyAtoms       @5 :List(CPMDDummyAtom);      # DUMMY ATOMS entries.
  changeBonds      @6 :List(Text);               # CHANGE BONDS lines.
  generate         @7 :Bool = false;             # GENERATE coordinates.
}

# @struct CPMDVdwSection
# @brief Typed &VDW empirical dispersion controls (new union arm; the
# directive-only vdw arm stays for wire compatibility).
struct CPMDVdwSection {
  empiricalCorrection @0 :Bool = false; # EMPIRICAL CORRECTION header.
  grimme     @1 :Text = "";             # D2 | D3 | D3BJ variant token.
  s6         @2 :Float64 = 0.0;         # Global scaling; 0 => omit.
  vdwCutoff  @3 :Float64 = 0.0;         # VDW-CUTOFF radius; 0 => omit.
  vdwCell    @4 :List(Int32);           # VDW-CELL replication n1 n2 n3.
  radius     @5 :Float64 = 0.0;         # RADIUS; 0 => omit.
  directives @6 :List(CPMDDirective);
}

# @struct CPMDHubbardU
# @brief Per-species DFT+U entry (structured replacement for the hubbard Text blob).
struct CPMDHubbardU {
  element @0 :Text;
  l       @1 :Int32 = 2;   # Angular channel.
  u       @2 :Float64;     # U value (eV).
}

struct CPMDDirectiveSection {
  directives @0 :List(CPMDDirective); # Keyword/value lines inside a named section.
  subsections @1 :List(CPMDGenericSection); # Nested NAME ... END NAME blocks.
}

enum CPMDSectionKind {
  generic @0;
  system  @1;
  cpmd    @2;
  dft     @3;
  atoms   @4;
  set     @5;
  raw     @6;
  atom    @7;
  basis   @8;
  clas    @9;
  eam     @10;
  exte    @11;
  hardness @12;
  info    @13;
  linres  @14;
  molstates @15;
  mts     @16;
  nlcc    @17;
  path    @18;
  pimd    @19;
  potential @20;
  prop    @21;
  ptddft  @22;
  resp    @23;
  tddft   @24;
  vdw     @25;
  vectors @26;
  wavefunction @27;
}

struct CPMDInputSection {
  union {
    generic @0 :CPMDGenericSection;
    system  @1 :CPMDSystemSection;
    cpmd    @2 :CPMDCpmdSection;
    dft     @3 :CPMDDftSection;
    atoms   @4 :CPMDAtomsSection;
    set     @5 :CPMDSetDirective;
    raw     @6 :Text; # Full &SECTION ... &END block text.
    atom    @7 :CPMDDirectiveSection;
    basis   @8 :CPMDDirectiveSection;
    clas    @9 :CPMDDirectiveSection;
    eam     @10 :CPMDDirectiveSection;
    exte    @11 :CPMDDirectiveSection;
    hardness @12 :CPMDDirectiveSection;
    info    @13 :CPMDDirectiveSection;
    linres  @14 :CPMDDirectiveSection;
    molstates @15 :CPMDDirectiveSection;
    mts     @16 :CPMDDirectiveSection;
    nlcc    @17 :CPMDDirectiveSection;
    path    @18 :CPMDDirectiveSection;
    pimd    @19 :CPMDDirectiveSection;
    potential @20 :CPMDDirectiveSection;
    prop    @21 :CPMDDirectiveSection;
    ptddft  @22 :CPMDDirectiveSection;
    resp    @23 :CPMDDirectiveSection;
    tddft   @24 :CPMDDirectiveSection;
    vdw     @25 :CPMDDirectiveSection;
    vectors @26 :CPMDDirectiveSection;
    wavefunction @27 :CPMDDirectiveSection;
    vdwParams @28 :CPMDVdwSection; # Typed &VDW; supersedes the directive-only vdw arm.
    propParams @29 :CPMDPropSection;     # Typed &PROP.
    linresParams @30 :CPMDLinresSection; # Typed &LINRES.
    pimdParams @31 :CPMDPimdSection;     # Typed &PIMD.
    pathParams @32 :CPMDPathSection;     # Typed &PATH.
    tddftParams @33 :CPMDTddftSection;   # Typed &TDDFT.
    respParams @34 :CPMDRespSection;     # Typed &RESP.
    exteParams @35 :CPMDExteSection;     # Typed &EXTE.
    vectorsParams @36 :CPMDVectorsSection; # Typed &VECTORS.
  }
}

# @struct CPMDParams
# @brief CPMD-specific knobs (backend arm inside PotentialConfig / rgpot params).
#
# Not a standalone config language — only used when the active potential is CPMD.
# Same fields in/out via Cap'n Proto. Geometry for each evaluation is ForceInput.
struct CPMDParams {
  functional   @0 :Text = "BLYP";     # Default DFT functional (also in dft section).
  cutOffRy     @1 :Float64 = 70.0;    # Plane-wave cutoff in Rydberg.
  charge       @2 :Int32 = 0;
  multiplicity @3 :Int32 = 1;         # 2S+1; >1 enables LSD when rendering defaults.
  task         @4 :Text = "gradient"; # energy | gradient | md | optimize (frontend hint).
  title        @5 :Text = "";         # Optional comment header in rendered deck.
  memoryMb     @6 :UInt32 = 0;        # Frontend memory hint; 0 => environment defaults.
  scratchDir   @7 :Text = "";         # Fallback &CPMD FILEPATH directory.
  permanentDir @8 :Text = "";         # Preferred &CPMD FILEPATH restart directory.
  cpmdRoot     @9 :Text = "";         # OpenCPMD source/build tree; empty => env CPMD_ROOT.
  enginePath   @10 :Text = "";        # Frontend engine selection hint.
  inputBlocks  @11 :List(Text);       # Raw &SECTION blocks applied before generated ones.
  inputSections @12 :List(CPMDInputSection); # Structured CPMD input sections.
  # Long-tail options: use CPMDInputSection.raw, inputBlocks, or set.
}

# @struct CommonMethodSpec
# @brief NOMAD-metainfo-aligned normalized method overlay (the thin Esperanto).
#
# Field names follow NOMAD's normalized method quantities (xc functional as
# libxc identifiers, Monkhorst-Pack k mesh, smearing kind/width, SCF
# thresholds). This struct lowers ONE-WAY into the native backend arm before
# it is applied; any knob the native arm also sets wins. It never replaces the
# arms -- the long tail of code-specific options stays code-native.
struct CommonMethodSpec {
  xcFunctionals @0 :List(Text);        # libxc names ("GGA_X_PBE", "GGA_C_PBE") or composites ("PBE", "B3LYP").
  basisSet      @1 :Text = "";         # Gaussian basis name, or "planewave".
  planewaveCutoffEv @2 :Float64 = 0.0; # >0 configures the plane-wave cutoff (eV).
  charge        @3 :Int32 = 0;
  spinMultiplicity @4 :Int32 = 0;      # 2S+1; 0 = unset (backend default).
  scfEnergyToleranceEv @5 :Float64 = 0.0; # SCF energy-change convergence (eV); 0 = unset.
  scfMaxIterations @6 :Int32 = 0;      # 0 = unset.
  kMesh @7 :List(Int32);               # Monkhorst-Pack divisions [3]; empty = molecular/Gamma.
  smearing @8 :Smearing;
  vanDerWaalsMethod @9 :Text = "";     # e.g. "DFT-D3", "DFT-D3(BJ)".
  relativityMethod  @10 :Text = "";    # e.g. "ZORA", "DKH".
  vanDerWaalsS6     @11 :Float64 = 0.0; # Global dispersion scaling s6; 0 = engine default.
                                        # Lowers to dft:vdw (NWChem) / vdwParams.s6 (CPMD).

  struct Smearing {
    kind    @0 :Kind = none;
    widthEv @1 :Float64 = 0.0;
    enum Kind {
      none             @0;
      fermi            @1;
      gaussian         @2;
      methfesselPaxton @3;
    }
  }
}


# @struct MetatomicParams
# @brief Metatomic (torch ML potential) backend arm, mirroring rgpot's
# MetatomicConfig so the model selection travels on the shared wire instead of
# rgpot-side ad-hoc config. lengthUnit/energyUnit stay on ForceInput.
struct MetatomicParams {
  modelPath            @0 :Text;            # Exported metatomic-torch model (.pt).
  device               @1 :Text = "cpu";    # torch device string, e.g. "cuda:0".
  extensionsDirectory  @2 :Text = "";       # TorchScript extensions dir.
  checkConsistency     @3 :Bool = false;    # metatomic consistency checks.
  uncertaintyThreshold @4 :Float64 = -1.0;  # <0 disables uncertainty gating (driver-side knob).
  dtypeOverride        @5 :Text = "";       # e.g. "float32"; empty => model dtype.
  outputs              @6 :List(Output);    # Requested model outputs; empty => ["energy"].

  # Mirrors upstream ModelOutput requests: standard names are "energy",
  # "energy_ensemble", "energy_uncertainty", "features",
  # "non_conservative_forces", "non_conservative_stress". Units ride on
  # ForceInput; selected_atoms is a per-step concern, not method config.
  struct Output {
    name              @0 :Text;
    perAtom           @1 :Bool = false;
    explicitGradients @2 :List(Text); # Gradient names, e.g. "positions", "strain".
  }
}

# @struct UmaParams
# @brief UMA / OMol AOTInductor backend arm, mirroring rgpot's UmaConfig so
# model selection travels on the shared wire. The package is a per-composition
# .pt2 export carrying its runtime contract (molecular_box / cutoff / z_set)
# as embedded AOTI metadata, which stays authoritative for values it holds.
struct UmaParams {
  modelPath    @0 :Text;          # AOTInductor .pt2 package path.
  taskName     @1 :Text = "omol"; # UMA task head.
  device       @2 :Text = "cpu";  # torch device string.
  charge       @3 :Int32 = 0;     # Total charge.
  spin         @4 :Int32 = 1;     # Spin multiplicity.
  cutoff       @5 :Float64 = 0.0; # Angstrom; <=0 keeps UmaConfig/sidecar value.
  maxNeighbors @6 :Int32 = 0;     # <=0 keeps UmaConfig/sidecar value.
}

# @struct PotentialSpec
# @brief One backend, named and configured, as a value that can be carried
# inside another backend's configuration.
#
# Backend selection is otherwise a property of which engine library a host
# dlopens, which works while a potential is a leaf. A learning surrogate is
# not a leaf: it wraps an oracle and has to say which one and how it is
# configured, so the choice has to travel on the wire rather than in the
# host's loader. This is that value, one arm per params struct already in
# this file.
#
# Purely additive: no existing field is renumbered, so a reader that does
# not know this struct is unaffected.
struct PotentialSpec {
  union {
    nwchem     @0 :NWChemParams;
    cpmd       @1 :CPMDParams;
    metatomic  @2 :MetatomicParams;
    uma        @3 :UmaParams;
    lammps     @4 :LammpsParams;
    # Backends whose configuration is entirely a library path plus the
    # geometry on ForceInput carry no params struct of their own.
    builtin    @5 :BuiltinParams;
  }
}

# @struct BuiltinParams
# @brief Arm for backends configured by name alone (analytic potentials and
# the dlopen-only engines), so PotentialSpec can name them without inventing
# a params struct that would hold one field.
struct BuiltinParams {
  name @0 :Text;   # e.g. "lj", "cuh2", "xtb", "morse".
  library @1 :Text = "";  # Optional explicit engine library path.
}

# @struct GprParams
# @brief Learning-surrogate arm: a Gaussian process that answers where its
# posterior is resolved and calls an oracle where it is not.
#
# The oracle is a nested PotentialSpec rather than a callback, so a host that
# selects this backend does not also have to supply ground truth through a
# side channel. That is the difference between an engine eOn can load and one
# eOn can use: its GenericEngineLoader resolves create/destroy/force and never
# calls set_callbacks, so an engine that expects its oracle by callback gets
# none.
#
# The gate is two-sided and stated in multiples of the model's own likelihood
# noise rather than in absolute energy, so it transfers across systems.
struct GprParams {
  oracle @0 :PotentialSpec;        # Ground truth. Required.
  sigmaLoFactor @1 :Float64 = 1.0;  # Answer from the model below this * sigma.
  sigmaHiFactor @2 :Float64 = 20.0; # Above this * sigma, do not learn the row.
  # A predicted force magnitude under this (eV/Angstrom) forces a
  # ground-truth call, so a driver's convergence claim is checked rather
  # than trusted. Non-positive disables the certificate, which leaves a
  # driver free to converge on the surrogate's own minimum.
  certifyBelowForce @3 :Float64 = 0.0;
  refitTrainingRevisions @4 :UInt64 = 5;  # Refit cadence in appends.
}

# @struct LammpsParams
# @brief LAMMPS backend arm driving the lammpc dynlib shim.
#
# lammpc calls the actual liblammps functions directly, in process: instance
# creation, pair setup and pair compute, position/force array access. No
# command-string parsing, no input scripts, no file traffic anywhere in the
# force loop. Geometry arrives per step via ForceInput; these fields fix the
# interaction model once per session.
struct LammpsParams {
  unitsStyle @0 :Text = "metal";  # LAMMPS units; metal => eV / Angstrom.
  atomStyle  @1 :Text = "atomic";
  pairStyle  @2 :Text = "";       # pair_style args, e.g. "eam/alloy" or "lj/cut 10.0".
  pairCoeffs @3 :List(Text);      # pair_coeff argument lines, one entry per line.
  typeToAtomicNumber @4 :List(Int32); # LAMMPS type i+1 -> Z; maps ForceInput.atmnrs.
  masses     @5 :List(Float64);   # Per-type mass (amu); empty => standard weight per Z.
  newtonPair @6 :Bool = true;
  boundary   @7 :Text = "p p p";  # boundary keyword args.
  extraSetup @8 :List(Text);      # Escape hatch: raw commands after the pair block.
  suffix     @9 :Text = "";       # Accelerator suffix (omp/gpu/kk) when compiled in.
}

# @struct NWChemDplotStanza
# @brief "dplot" block: density/orbital cube output.
struct NWChemDplotStanza {
  title       @0 :Text = "";
  gaussianCube @1 :Text = "";     # gaussian <file> output.
  civecs      @2 :Text = "";      # TDDFT transition-density source.
  limitXyz    @3 :List(Float64);  # LimitXYZ extents (6 or 9 values).
  spin        @4 :Text = "";      # total | alpha | beta | spindens.
  orbitals    @5 :List(Int32);    # orbitals view list.
  output      @6 :Text = "";      # output <file>.
  directives  @7 :List(NWChemDirective);
}

# @struct NWChemEspStanza
# @brief "esp" block: ESP / RESP atomic charges.
struct NWChemEspStanza {
  recalculate   @0 :Bool = false;
  rangeFactor   @1 :Float64 = 0.0; # range <f>; 0 => omit.
  probe         @2 :Float64 = 0.0; # probe <r>; 0 => omit.
  spacing       @3 :Float64 = 0.0; # spacing <s>; 0 => omit.
  restrain      @4 :Bool = false;  # RESP restraint.
  restrainHfree @5 :Bool = false;  # restrain hfree.
  directives    @6 :List(NWChemDirective);
}

# ---- Typed CPMD long-tail sections (upgrade from directive-only arms) ----

# @struct CPMDPropSection
# @brief Typed &PROP: post-SCF property requests.
struct CPMDPropSection {
  dipoleMoment   @0 :Bool = false; # DIPOLE MOMENT.
  localize       @1 :Bool = false; # LOCALIZE (Wannier).
  polarizability @2 :Bool = false; # POLARIZABILITY.
  cubecenter     @3 :List(Float64); # CUBECENTER x y z.
  cubefileDensity @4 :Bool = false; # CUBEFILE DENSITY.
  cubefileOrbitals @5 :List(Int32); # CUBEFILE ORBITALS list.
  chargesGaussian @6 :Bool = false; # CHARGES.
  directives     @7 :List(CPMDDirective);
}

# @struct CPMDLinresSection
# @brief Typed &LINRES: linear-response convergence controls.
struct CPMDLinresSection {
  convergence @0 :Float64 = 0.0; # CONVERGENCE; 0 => omit.
  maxSteps    @1 :Int32 = 0;     # MAXSTEP; 0 => omit.
  hthrs       @2 :Float64 = 0.0; # HTHRS preconditioner threshold; 0 => omit.
  directives  @3 :List(CPMDDirective);
}

# @struct CPMDPimdSection
# @brief Typed &PIMD: path-integral MD controls.
struct CPMDPimdSection {
  replicas    @0 :Int32 = 0;     # TROTTER DIMENSION / replica count; 0 => omit.
  facmass     @1 :Float64 = 0.0; # FACMASS; 0 => omit.
  centroidDynamics @2 :Bool = false; # CENTROID DYNAMICS.
  normalModes @3 :Bool = false;  # NORMAL MODES transform.
  staging     @4 :Bool = false;  # STAGING transform.
  directives  @5 :List(CPMDDirective);
}

# @struct CPMDPathSection
# @brief Typed &PATH: NEB / string minimum-energy-path controls.
struct CPMDPathSection {
  replicaNumber @0 :Int32 = 0;    # REPLICA NUMBER; 0 => omit.
  nebSpring     @1 :Float64 = 0.0; # NEB spring constant via directives; 0 => omit.
  factor        @2 :Float64 = 0.0; # FACTOR step scaling; 0 => omit.
  nloop         @3 :Int32 = 0;     # NLOOP path iterations; 0 => omit.
  alpha         @4 :Float64 = 0.0; # ALPHA mixing; 0 => omit.
  directives    @5 :List(CPMDDirective);
}

# @struct CPMDTddftSection
# @brief Typed &TDDFT: excited-state controls.
struct CPMDTddftSection {
  states       @0 :Int32 = 0;     # STATES n; 0 => omit.
  tammDancoff  @1 :Bool = false;  # TAMM-DANCOFF.
  diagonalizer @2 :Text = "";     # DIAGONALIZER token (DAVIDSON, ...).
  convergence  @3 :Float64 = 0.0; # CONVERGENCE; 0 => omit.
  directives   @4 :List(CPMDDirective);
}

# @struct NWChemQmdStanza
# @brief "qmd" block: Born-Oppenheimer molecular dynamics controls.
struct NWChemQmdStanza {
  nstepNucl   @0 :Int32 = 0;     # nstep_nucl; 0 => omit.
  dtNucl      @1 :Float64 = 0.0; # dt_nucl (au); 0 => omit.
  targTemp    @2 :Float64 = 0.0; # targ_temp (K); 0 => omit.
  thermostat  @3 :Text = "";     # none|berendsen|langevin|svr|rescale.
  thermostatParameter @4 :Float64 = 0.0; # tau / friction; 0 => omit.
  randomSeed  @5 :Int32 = 0;     # rand_seed; 0 => omit.
  comStep     @6 :Int32 = 0;     # com_step; 0 => omit.
  printXyz    @7 :Int32 = 0;     # print_xyz cadence; 0 => omit.
  linear      @8 :Bool = false;  # linear molecule flag.
  directives  @9 :List(NWChemDirective);
}

# @struct NWChemRamanStanza
# @brief "raman" block: Raman intensity controls (with task dft raman).
struct NWChemRamanStanza {
  normalMode  @0 :Bool = false;   # NORMAL vs resonance.
  dynamic     @1 :Bool = false;   # DYNAMIC polarizability mode.
  low         @2 :Float64 = 0.0;  # LOW frequency bound (cm-1); 0 => omit.
  high        @3 :Float64 = 0.0;  # HIGH bound; 0 => omit.
  steps       @4 :Int32 = 0;      # STEPS in the range; 0 => omit.
  fieldFrequency @5 :Float64 = 0.0; # Incident field frequency (au); 0 => omit.
  directives  @6 :List(NWChemDirective);
}

# @struct NWChemFonStanza
# @brief Fractional occupation (fon) controls promoted through dft.
struct NWChemFonStanza {
  partial     @0 :Int32 = 0;     # Number of partially occupied orbitals.
  electrons   @1 :Float64 = 0.0; # Electrons distributed over the partial set.
  filled      @2 :Int32 = 0;     # Fully occupied orbitals below the set.
  alphaOnly   @3 :Bool = false;  # Apply to alpha spin channel only.
  temperature @4 :Float64 = 0.0; # Electronic temperature (K); 0 => omit.
  directives  @5 :List(NWChemDirective);
}

# ---- CPMD typed long-tail batch 3 ----

# @struct CPMDRespSection
# @brief Typed &RESP: restrained electrostatic potential charges.
struct CPMDRespSection {
  restraintStrength @0 :Float64 = 0.0; # WEIGHT; 0 => omit.
  hyperbolic        @1 :Bool = false;  # HYPERBOLIC restraint form.
  tightness         @2 :Float64 = 0.0; # BETA tightness; 0 => omit.
  directives        @3 :List(CPMDDirective);
}

# @struct CPMDExteSection
# @brief Typed &EXTE: external field / interface controls.
struct CPMDExteSection {
  externalField @0 :List(Float64); # EFIELD vector components.
  directives    @1 :List(CPMDDirective);
}

# @struct CPMDVectorsSection
# @brief Typed &VECTORS: initial-guess and orthogonalization controls.
struct CPMDVectorsSection {
  newOrtho   @0 :Bool = false; # NEWORTHO.
  overlap    @1 :Bool = false; # OVERLAP diagnostics.
  directives @2 :List(CPMDDirective);
}

# @struct NWChemNebStanza
# @brief "neb" block; keywords mined from optim/neb/neb_input.F.
struct NWChemNebStanza {
  nbeads     @0 :Int32 = 0;     # neb:nbeads; 0 => omit.
  kbeads     @1 :Float64 = 0.0; # neb:kbeads spring constant; 0 => omit.
  maxiter    @2 :Int32 = 0;     # neb:steps; 0 => omit.
  stepsize   @3 :Float64 = 0.0; # neb:stepsize; 0 => omit.
  trust      @4 :Float64 = 0.0; # neb:trust radius; 0 => omit.
  nhist      @5 :Int32 = 0;     # neb:m history depth; 0 => omit.
  algorithm  @6 :Int32 = 0;     # neb:algorithm; 0 => omit.
  reset      @7 :Bool = false;  # neb:nebnew restart reset.
  convergence @8 :Convergence = unspecified; # loose | default | tight presets.
  gmax       @9 :Float64 = 0.0;  # Explicit neb:gmax; 0 => omit.
  grms       @10 :Float64 = 0.0; # neb:grms.
  xmax       @11 :Float64 = 0.0; # neb:xmax.
  xrms       @12 :Float64 = 0.0; # neb:xrms.
  directives @13 :List(NWChemDirective);
  enum Convergence {
    unspecified @0;
    loose       @1;
    default     @2;
    tight       @3;
  }
}

# @struct NWChemMepgsStanza
# @brief Minimum-energy-path Gonzalez-Schlegel IRC block; mepgs:*/gsopt:*
# RTDB keys from optim/mepgs/mepgs_input.F.
struct NWChemMepgsStanza {
  stride     @0 :Float64 = 0.0; # mepgs:stride step along the IRC; 0 => omit.
  maxmep     @1 :Int32 = 0;     # mepgs:nircopt IRC points; 0 => omit.
  maxiter    @2 :Int32 = 0;     # gsopt:nptopt per-point optimizer steps.
  inhess     @3 :Int32 = 0;     # mepgs:inhess initial hessian source.
  evib       @4 :Float64 = 0.0; # mepgs:evib vibrational energy step.
  direction  @5 :Direction = both; # forward/backward selection.
  opttol     @6 :Float64 = 0.0; # gsopt:opt_tol convergence.
  eprec      @7 :Float64 = 0.0; # gsopt:eprec energy precision.
  xyz        @8 :Bool = false;  # write xyz trajectory files.
  mswg       @9 :Bool = false;  # ircgs:mswg mass-weighted coordinates.
  directives @10 :List(NWChemDirective);
  enum Direction {
    both     @0;
    forward  @1;
    backward @2;
  }
}

# @struct NWChemTroptStanza
# @brief Trust-region optimizer block; tropt:* RTDB keys from
# optim/tropt/tropt_input.F.
struct NWChemTroptStanza {
  opttol     @0 :Float64 = 0.0;  # tropt:opt_tol.
  eprec      @1 :Float64 = 0.0;  # tropt:eprec.
  trust      @2 :Float64 = 0.0;  # tropt:trust radius.
  maxiter    @3 :Int32 = 0;      # tropt:nptopt.
  inhess     @4 :Int32 = 0;      # tropt:inhess.
  lbfgs      @5 :Bool = false;   # tropt:lbfgs limited-memory update.
  mh         @6 :Int32 = 0;      # tropt:mh lbfgs history depth.
  linopt     @7 :Int32 = 0;      # tropt:linopt line-search mode.
  qstep      @8 :NWChemToggle = unspecified; # qstep/noqstep quadratic step.
  moddir     @9 :Int32 = 0;      # tropt:moddir mode to follow (saddle).
  modsad     @10 :Int32 = 0;     # tropt:modsad saddle mode count.
  redoautoz  @11 :Bool = false;  # tropt:redoautoz rebuild internals.
  convergence @12 :Convergence = unspecified; # loose/default/tight presets.
  gmax       @13 :Float64 = 0.0; # explicit tropt:gmax_tol.
  grms       @14 :Float64 = 0.0; # tropt:grms_tol.
  xmax       @15 :Float64 = 0.0; # tropt:xmax_tol.
  xrms       @16 :Float64 = 0.0; # tropt:xrms_tol.
  xyz        @17 :Bool = false;  # write xyz per step.
  directives @18 :List(NWChemDirective);
  enum Convergence {
    unspecified @0;
    loose       @1;
    default     @2;
    tight       @3;
  }
}

# @struct NWChemMetadynamicsCv
# @brief One collective variable line inside the nwpw metadynamics block
# (metadynamics_input.F). bond/angle/dihedral take bare atom indices; the
# coordination number takes two index groups plus n/m/r0 switching function.
struct NWChemMetadynamicsCv {
  kind      @0 :Kind = bond;
  atoms     @1 :List(Int32);   # bond: 2, angle: 3, dihedral: 4 indices.
  index1    @2 :List(Int32);   # coordNumber group 1.
  index2    @3 :List(Int32);   # coordNumber group 2.
  weight    @4 :Float64 = 0.0; # w hill height (au); 0 => engine default 5e-5.
  sigma     @5 :Float64 = 0.0; # hill width (au); 0 => default 0.1.
  rangeLow  @6 :Float64 = 0.0; # "range <lo> <hi>"; both 0 => omit.
  rangeHigh @7 :Float64 = 0.0;
  nrange    @8 :Int32 = 0;     # grid points; 0 => default 501.
  n         @9 :Float64 = 0.0; # coordNumber exponent; 0 => default 6.
  m         @10 :Float64 = 0.0; # 0 => default 12.
  r0        @11 :Float64 = 0.0; # switching radius; 0 => default 3.0.
  enum Kind {
    bond        @0;
    angle       @1;
    dihedral    @2;
    coordNumber @3;
  }
}

# @struct NWChemMetadynamicsStanza
# @brief nwpw metadynamics block; renders nested inside "nwpw ... end".
struct NWChemMetadynamicsStanza {
  cvs               @0 :List(NWChemMetadynamicsCv);
  tempered          @1 :Bool = false; # well-tempered variant.
  boundary          @2 :Float64 = 0.0; # boundary potential; 0 => omit.
  potentialFilename @3 :Text = "";
  printShift        @4 :Int32 = 0;
  directives        @5 :List(NWChemDirective); # equation/bondings etc. stay literal.
}

# @struct NWChemCellOptimizeStanza
# @brief nwpw cell_optimize block (cell_optimize_input.F); renders nested
# inside "nwpw ... end".
struct NWChemCellOptimizeStanza {
  cellName          @0 :Text = "";   # cell_optimize:cell_name; "" => cell_default.
  strategy          @1 :Strategy = unspecified; # "optimize all|lattice|lattice_vectors".
  optimizeLattice   @2 :List(Text);  # lat_a lat_b lat_c alpha beta gamma subset.
  optimizeLatticeVectors @3 :List(Text); # a11..a33 subset.
  latticeTolerances @4 :List(Float64);   # up to 6 values.
  cycles            @5 :Int32 = 0;   # 0 => engine default 10.
  directives        @6 :List(NWChemDirective);
  enum Strategy {
    unspecified    @0;
    all            @1;
    lattice        @2;
    latticeVectors @3;
  }
}

# @struct NWChemRismStanza
# @brief 3D-RISM solvation block; rism:* RTDB keys from rism/rism_input.F.
struct NWChemRismStanza {
  permittivity @0 :Float64 = 0.0; # rism:permittivity.
  tau          @1 :Float64 = 0.0; # rism:tau.
  temp         @2 :Float64 = 0.0; # rism:temp (K).
  tol          @3 :Float64 = 0.0; # rism:tol convergence.
  lambda       @4 :Float64 = 0.0; # rism:lambda coupling.
  diis         @5 :Int32 = 0;     # rism:diis subspace size.
  ngrid        @6 :Int32 = 0;     # rism:ngrid radial grid points.
  closure      @7 :Text = "";     # rism:closure (hnc, kh, py).
  vdwRule      @8 :Text = "";     # rism:vdw:rule combination rule.
  vdwParameters @9 :Text = "";    # rism:vdw:parameters file.
  directives   @10 :List(NWChemDirective); # solute/solvent subblocks stay literal.
}

# @struct NWChemDimQmStanza
# @brief Discrete interaction model / QM block; dimqm:* keys from dimqm/dimqm_input.F.
struct NWChemDimQmStanza {
  frequency  @0 :Bool = false;    # frequency-dependent (complex) polarizabilities.
  algorithm  @1 :Int32 = 0;       # dimqm:alg; 0 => engine default.
  tolerance  @2 :Float64 = 0.0;   # dimqm:dimtol induced-dipole tolerance.
  efield     @3 :List(Float64);   # finite field [x, y, z]; empty => none.
  localfield @4 :Bool = false;    # local-field effects.
  noresp     @5 :Bool = false;    # disable DIM response (dimqm:lrsp=false).
  screen     @6 :Screen = unspecified; # dipole-interaction screening.
  screenFactor @7 :Float64 = 0.0; # factor for exp/erf screening.
  noseed     @8 :Bool = false;    # disable seeding.
  debug      @9 :Bool = false;
  printAtomicDipoles @10 :Bool = false;
  directives @11 :List(NWChemDirective);
  enum Screen {
    unspecified @0;
    none        @1;
    exp         @2;
    erf         @3;
  }
}

# @struct NWChemGwStanza
# @brief Molecular GW block; grammar and gw:* RTDB keys from gwmol/gw_input.F.
struct NWChemGwStanza {
  evgw       @0 :Bool = false;   # gw:evgw eigenvalue self-consistency.
  evgw0      @1 :Bool = false;   # gw:evgw0 (G updated, W0 frozen).
  eviter     @2 :Int32 = 0;      # gw:eviter; 0 => omit count.
  method     @3 :Method = unspecified; # gw:cdgw contour deformation vs analytic.
  ngl        @4 :Int32 = 0;      # gw:ngl Gauss-Legendre grid points (cdgw).
  solver     @5 :Solver = unspecified; # gw:graph vs gw:newton root search.
  qpiter     @6 :Int32 = 0;      # gw:qpiter Newton iterations.
  thresholdEv @7 :Float64 = 0.0; # gw:threshold; rendered "convergence <x> ev".
  occAlpha   @8 :Int32 = 0;      # gw:noqp_alpha states alpha occ.
  virAlpha   @9 :Int32 = 0;      # gw:nvqp_alpha states alpha vir.
  occBeta    @10 :Int32 = 0;     # gw:noqp_beta.
  virBeta    @11 :Int32 = 0;     # gw:nvqp_beta.
  first      @12 :Int32 = 0;     # gw:first orbital counting offset.
  ncap       @13 :Bool = false;  # gw:ncap modified eigenvalue start.
  eta        @14 :Float64 = 0.0; # gw:eta imaginary infinitesimal.
  rpa        @15 :Bool = false;  # gw:rpa correlation energy.
  diag       @16 :Bool = false;  # gw:diag invert by diagonalization.
  core       @17 :Bool = false;  # gw:core count from core upwards.
  directives @18 :List(NWChemDirective);
  enum Method {
    unspecified @0;
    cdgw        @1;
    analytic    @2;
  }
  enum Solver {
    unspecified @0;
    graph       @1;
    newton      @2;
  }
}

# @struct NWChemEtransStanza
# @brief Electron-transfer coupling block; et:* RTDB keys from etrans/et_input.F.
struct NWChemEtransStanza {
  tol2e      @0 :Float64 = 0.0;  # et:tol2e two-electron screening threshold.
  fock       @1 :NWChemToggle = unspecified; # fock vs nofock (et:method_2e).
  fmo        @2 :Bool = false;   # et:fmo fragment MO coupling.
  vectorsReactant @3 :Text = ""; # "vectors <reactant.movecs> <product.movecs>".
  vectorsProduct  @4 :Text = "";
  directives @5 :List(NWChemDirective);
}

# @struct NWChemStringStanza
# @brief Zero-temperature string method; keys from optim/string/string_input.F.
struct NWChemStringStanza {
  nbeads     @0 :Int32 = 0;     # string:nbeads; 0 => omit.
  maxiter    @1 :Int32 = 0;     # string:maxit.
  stepsize   @2 :Float64 = 0.0; # string:stepsize.
  tol        @3 :Float64 = 0.0; # string:tol gradient tolerance.
  interpol   @4 :Int32 = 0;     # string:interpol interpolation order.
  nhist      @5 :Int32 = 0;     # string:m history depth.
  algorithm  @6 :Int32 = 0;     # string:algorithm.
  printShift @7 :Int32 = 0;     # string:print_shift; 0 => omit.
  freeze1    @8 :Bool = false;  # string:freeze1 pin first bead.
  freezen    @9 :Bool = false;  # string:freezen pin last bead.
  hasmiddle  @10 :Bool = false; # string:hasmiddle middle bead supplied.
  impose     @11 :Bool = false; # string:impose re-impose geometry constraints.
  reset      @12 :Bool = false; # bead_list:new restart reset.
  xyzPath    @13 :Text = "";    # xyz_path initial path file.
  mode       @14 :Mode = unspecified; # string:mode execution mode.
  directives @15 :List(NWChemDirective);
  enum Mode {
    unspecified @0;
    serial      @1;
    parallel    @2;
  }
}

# @struct NWChemCdftConstraint
# @brief One dft-block cdft line: atom ranges, charge|spin, target value.
struct NWChemCdftConstraint {
  firstAtomStart  @0 :Int32;
  firstAtomEnd    @1 :Int32;
  secondAtomStart @2 :Int32 = 0; # 0 => single-range constraint.
  secondAtomEnd   @3 :Int32 = 0;
  kind            @4 :Kind = charge;
  value           @5 :Float64;
  enum Kind {
    charge @0;
    spin   @1;
  }
}

# @struct NWChemCdftSpec
# @brief Constrained-DFT controls inside the dft block (cdft_inp.F).
struct NWChemCdftSpec {
  constraints @0 :List(NWChemCdftConstraint);
  convergence @1 :Float64 = 0.0; # dft:dl_conv; 0 => default 1e-6.
  population  @2 :Population = unspecified; # dft:ipop scheme.
  enum Population {
    unspecified @0;
    becke       @1;
    mulliken    @2;
    lowdin      @3;
  }
}

# @struct ExprTerm
# @brief One named child in an ExprParams composition tree.
#
# `name` is the Lepton identifier used in ExprParams.expression.
# `pot` is a PotSpec so a child may itself be an expression.
struct ExprTerm {
  name @0 :Text;
  pot @1 :PotSpec;
}

# @struct ExprParams
# @brief Recursive PES algebra for RPC/config. Not a PotentialConfig arm.
#
# The C++ ExprPot constructor (expression string plus named
# unique_ptr<PotentialBase> children) is the first-slice API.
# This tree names the same composition on the wire (e.g. "0.5*lj + d3")
# without a flat scale list. Leaf PotSpec arms land when those pots
# grow config structs; do not invent dummy leaves.
struct ExprParams {
  expression @0 :Text;
  terms @1 :List(ExprTerm);
}

# @struct PotSpec
# @brief Recursive potential tree used only for composition.
#
# Not a PotentialConfig arm. PotentialConfig stays a single-backend
# configure union (none / nwchem / cpmd / metatomic).
struct PotSpec {
  union {
    expr @0 :ExprParams;
    none @1 :Void; # Empty node. Leaf arms land when those pots grow config structs.
  }
}

# @struct PotentialConfig
# @brief **rgpot user parameters (extensible, in/out via Cap'n Proto only).**
#
# This is the single user-facing options carrier for rgpot: pass in to configure
# a live Potential (RPC `configure` or in-process apply), and/or round-trip out
# when a backend supports get. One schema for wire + embed - no parallel
# TOML/JSON/YAML option files for backends.
#
# Tagged union: exactly one backend's options (or none). Add new arms as new
# potentials gain runtime knobs (e.g. metatomic @2 :MetatomicParams).
# `calculate` geometry stays on ForceInput; this struct is method/backend setup only.
struct PotentialConfig {
  union {
    none      @0 :Void;         # No backend-specific options (or no-op configure).
    nwchem    @1 :NWChemParams; # NWChemPot / potserv ... NWChem
    cpmd      @2 :CPMDParams;   # CPMDPot / potserv ... OpenCPMD
    metatomic @4 :MetatomicParams; # MetatomicPot (torch ML models)
    lammps    @5 :LammpsParams;    # lammpc dynlib shim (LAMMPS library API)
    # xtb       @6 :XTBParams;
    # tblite    @7 :TBLiteParams;
  }
  # Normalized overlay applied BEFORE the native arm; native settings win.
  common @3 :CommonMethodSpec;
}
# @interface Potential
# @brief The RPC interface for remote calculations.
interface Potential {
  # @brief Executes the potential and force calculation.
  # @param fip The input atomic configuration.
  # @return The resulting energy and force vector.
  calculate @0 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Apply rgpot user parameters (PotentialConfig) before calculate().
  # @param config Backend-tagged options (nwchem, future metatomic/xtb/...).
  # @return ok=false if the arm does not match the server backend or apply failed.
  configure @1 (config :PotentialConfig) -> (ok :Bool, message :Text);

  # @brief Compute only the scalar potential energy for this geometry.
  # @return PotentialResult.energy in ForceInput.energyUnit.
  calculateEnergy @2 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute potential energy and force vector for this geometry.
  # @return PotentialResult.energy and PotentialResult.forces.
  calculateForces @3 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute a dense Cartesian Hessian for this geometry.
  # @return PotentialResult.hessian in energyUnit / lengthUnit^2.
  calculateHessian @4 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute the total electric dipole vector for this geometry.
  # @return PotentialResult.dipole in atomic units.
  calculateDipole @5 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute the total traceless electric quadrupole tensor.
  # @return PotentialResult.quadrupole in atomic units as xx, xy, xz, yy, yz, zz.
  calculateQuadrupole @6 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute the periodic stress tensor for this geometry.
  # @return PotentialResult.stress in energyUnit / lengthUnit^3.
  calculateStress @7 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Optimize this geometry and return final coordinates.
  # @return PotentialResult.optimizedPos in ForceInput.lengthUnit.
  calculateOptimize @8 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute harmonic vibrational frequencies and IR intensities.
  # @return PotentialResult.frequencies in cm^-1 and intensities in atomic units.
  calculateFrequencies @9 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute the electric polarizability response vector.
  # @return PotentialResult.polarizability stores NWChem aoresponse:alpha.
  calculatePolarizability @10 (fip :ForceInput) -> (result :PotentialResult);

  # @brief Compute potential energy and nuclear gradient for this geometry.
  # @return PotentialResult.energy and PotentialResult.gradient.
  calculateGradient @11 (fip :ForceInput) -> (result :PotentialResult);
}
