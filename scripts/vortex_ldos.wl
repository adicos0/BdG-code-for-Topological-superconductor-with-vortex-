(* Adapted from original_code/vortex_ldos_original.m and added_vortex.wl.
   Loading this file defines functions only; no calculation starts automatically.
   See notebook 05 and LDOS_WORKFLOW.md for conventions and provenance. *)

(* SECTION: setup *)
$BdGProjectRoot = DirectoryName[DirectoryName[$InputFileName]];
ClearAll[bdgRequire, bdgNorm, SiteNumber, SiteCoordinates];
bdgRequire[test_, message_] := If[! TrueQ[test],
  Throw[Failure["BdGInput", <|"MessageTemplate" -> message|>]]];
bdgNorm[a_] := Max[Abs[Flatten[Normal[a]]]];
SiteNumber[{x_, y_}, ny_] := (x - 1) ny + y;
SiteCoordinates[i_, ny_] := {Quotient[i - 1, ny] + 1, Mod[i - 1, ny] + 1};

(* SECTION: import *)
ClearAll[LoadBdGModel];
LoadBdGModel[folder_, nx_Integer, ny_Integer] := Module[
  {n = nx ny, names, blocks, sc, zero, d, coords},
  bdgRequire[nx >= 3 && ny >= 3, "Each lattice dimension must be at least 3."];
  names = {"Huu", "Hud", "Hdu", "Hdd", "SCs"};
  bdgRequire[And @@ (FileExistsQ[FileNameJoin[{folder, # <> ".wdx"}]] & /@ names),
    "Five WDX files are required in the input folder."];
  blocks = AssociationThread[names, Import[FileNameJoin[{folder, # <> ".wdx"}]] & /@ names];
  bdgRequire[And @@ (MatrixQ[Normal[blocks[#]], NumericQ] & /@ names), "Non-numeric matrix."];
  bdgRequire[And @@ (Dimensions[blocks[#]] === {n, n} & /@ Take[names, 4]),
    "Normal blocks do not match Nx Ny; do not guess the lattice shape."];
  sc = Normal[blocks["SCs"]];
  bdgRequire[Dimensions[sc] === {2 n, 2 n}, "SCs must have size 2 NN by 2 NN."];
  zero = ConstantArray[0., {n, n}];
  d = sc[[1 ;; n, n + 1 ;; 2 n]];
  bdgRequire[bdgNorm[sc - ArrayFlatten[{{zero, d}, {d, zero}}]] < 10^-10,
    "Expected the symmetric off-diagonal pairing layout used by the original files."];
  bdgRequire[bdgNorm[d - DiagonalMatrix[Diagonal[d]]] < 10^-10,
    "This workflow supports onsite pairing only."];
  bdgRequire[bdgNorm[blocks["Huu"] - ConjugateTranspose[blocks["Huu"]]] < 10^-10 &&
    bdgNorm[blocks["Hdd"] - ConjugateTranspose[blocks["Hdd"]]] < 10^-10 &&
    bdgNorm[blocks["Hdu"] - ConjugateTranspose[blocks["Hud"]]] < 10^-10,
    "Imported normal blocks are not Hermitian partners."];
  coords = SiteCoordinates[#, ny] & /@ Range[n];
  Join[blocks, <|"Nx" -> nx, "Ny" -> ny, "NN" -> n,
    "Coordinates" -> coords, "Delta" -> Diagonal[d], "InputFolder" -> folder|>]
];

(* SECTION: impurities *)
ClearAll[AddImpurity];
AddImpurity[model_, kind_, strength_, site_Integer] := Module[{n = model["NN"], p},
  bdgRequire[MemberQ[{"none", "scalar", "magnetic"}, kind], "Unknown impurity kind."];
  bdgRequire[1 <= site <= n && NumericQ[strength] && TrueQ[Im[strength] == 0], "Invalid impurity."];
  p = SparseArray[{{site, site} -> N[strength]}, {n, n}];
  Join[model, <|"Huu" -> model["Huu"] + If[kind === "none", 0 p, p],
    "Hdd" -> model["Hdd"] + Switch[kind, "scalar", p, "magnetic", -p, _, 0 p],
    "ImpurityKind" -> kind, "ImpurityStrength" -> strength, "ImpuritySite" -> site|>]
];

(* SECTION: boundary *)
ClearAll[BoundaryPhases];
BoundaryPhases[model_] := Module[{c = model["Coordinates"], n = model["NN"], lengths},
  lengths = {model["Nx"], model["Ny"]};
  Table[Table[With[{d = c[[j, axis]] - c[[i, axis]]},
    If[Abs[d] > lengths[[axis]]/2, Sign[d], 0]], {i, n}, {j, n}], {axis, 2}]
];

(* SECTION: matrix *)
ClearAll[BdGMatrix];
BdGMatrix[model_, delta_, {kx_, ky_}, phases_, shift_: 0.] := Module[
  {n = model["NN"], phase, normal, minusNormal, s, zero, d, pair, h},
  phase = Exp[I (kx model["Nx"] phases[[1]] + ky model["Ny"] phases[[2]])];
  normal[p_] := ArrayFlatten[{
    {(Normal[model["Huu"]] + shift IdentityMatrix[n]) p, Normal[model["Hud"]] p},
    {Normal[model["Hdu"]] p, (Normal[model["Hdd"]] + shift IdentityMatrix[n]) p}}];
  minusNormal = normal[Conjugate[phase]];
  s = DiagonalMatrix[Join[ConstantArray[1., n], ConstantArray[-1., n]]];
  zero = ConstantArray[0., {n, n}];
  d = DiagonalMatrix[delta];
  pair = ArrayFlatten[{{zero, d}, {d, zero}}];
  h = ArrayFlatten[{{normal[phase], pair},
    {ConjugateTranspose[pair], -s . Conjugate[minusNormal] . s}}];
  bdgRequire[bdgNorm[h - ConjugateTranspose[h]] < 10^-9, "BdG matrix is not Hermitian."];
  N[h]
];

(* SECTION: vortex *)
ClearAll[VortexFactor];
VortexFactor[model_, core_Integer, coherenceSites_, phaseMode_: "Winding"] := Module[
  {centre, offsets},
  bdgRequire[1 <= core <= model["NN"] && coherenceSites > 0, "Invalid vortex centre or length."];
  bdgRequire[MemberQ[{"Winding", "OriginalAbsoluteAngle"}, phaseMode], "Unknown phase mode."];
  centre = SiteCoordinates[core, model["Ny"]];
  offsets = (# - centre) & /@ model["Coordinates"];
  Map[Function[r, If[r === {0, 0}, 0.,
    Tanh[Norm[r]/coherenceSites] Exp[I If[phaseMode === "Winding",
      ArcTan[r[[1]], r[[2]]],
      If[r[[1]] == 0, Pi/2, ArcTan[Abs[r[[2]]]/Abs[r[[1]]]]]]]]], offsets]
];

(* SECTION: relaxation *)
ClearAll[RelaxPairing];
RelaxPairing[model_, initial_, factor_, phases_, settings_] := Module[
  {delta = N[initial], shift = 0., h, energies, vectors, f, proposal,
   next, nextShift, density, residual, history = {}, converged = False,
   n = model["NN"], temperature = settings["kT"], j},
  bdgRequire[temperature > 0 && settings["MaxIterations"] > 0, "Invalid iteration settings."];
  Do[
    h = BdGMatrix[model, delta, {0., 0.}, phases, shift];
    {energies, vectors} = Eigensystem[h];
    (* Same clipped Fermi function, component pairing sum and mixing as the saved solver. *)
    f = (If[# > .02, 0., If[# < -.02, 1., 1./(Exp[#/temperature] + 1.)]] &) /@ energies;
    density = Mean[Table[Abs[vectors[[All, j]]]^2 . f +
      Abs[vectors[[All, j + 3 n]]]^2 . (1. - f), {j, n}]];
    proposal = settings["Coupling"] Table[
      (vectors[[All, j]] Conjugate[vectors[[All, j + 3 n]]] +
       vectors[[All, j + n]] Conjugate[vectors[[All, j + 2 n]]]) . f, {j, n}] factor;
    residual = Max[Abs[proposal - delta]];
    next = (proposal + delta)/2.;
    nextShift = shift - .1 (settings["TargetDensity"] - density);
    AppendTo[history, {iteration, residual, density, shift}];
    converged = residual < settings["Tolerance"] &&
      Abs[density - settings["TargetDensity"]] < settings["DensityTolerance"];
    delta = next; shift = nextShift;
    If[converged, Break[]], {iteration, settings["MaxIterations"]}];
  <|"Delta" -> delta, "Shift" -> shift, "Converged" -> converged, "History" -> history|>
];

(* SECTION: ldos *)
ClearAll[CalculateLDOS];
CalculateLDOS[model_, delta_, phases_, sites_, energies_, eta_, ns_Integer, shift_: 0.] := Module[
  {n = model["NN"], rho, h, ev, vec, weights, spectral, kx, ky},
  bdgRequire[eta > 0 && ns >= 1 && Length[sites] > 0 &&
    And @@ (IntegerQ[#] && 1 <= # <= n & /@ sites), "Invalid LDOS grid or sites."];
  rho = ConstantArray[0., {Length[energies], Length[sites]}];
  Do[
    kx = 2 Pi a/(ns model["Nx"]); ky = 2 Pi b/(ns model["Ny"]);
    h = BdGMatrix[model, delta, {kx, ky}, phases, shift];
    {ev, vec} = Eigensystem[h];
    (* Full spectrum + electron components only: one state per site per spin on average. *)
    weights = (Abs[vec[[All, sites]]]^2 + Abs[vec[[All, sites + n]]]^2)/2.;
    spectral = Table[eta/(Pi ((energy - ev)^2 + eta^2)), {energy, energies}];
    rho += spectral . weights,
    {a, 0, ns - 1}, {b, 0, ns - 1}];
  rho/ns^2
];

(* SECTION: plots *)
ClearAll[LDOSCurves, PairingViews, SpectrumView, LDOSMap];
LDOSCurves[energies_, rho_, sites_, references_: {}, unit_: "model energy"] :=
 ListLinePlot[Table[Transpose[{energies, rho[[All, j]]}], {j, Length[sites]}],
  Frame -> True, FrameLabel -> {"Energy (" <> unit <> ")", "Spin-averaged LDOS (1/" <> unit <> ")"},
  PlotLegends -> ("site " <> ToString[#] & /@ sites), PlotStyle -> {Green, Blue, Orange},
  GridLines -> {references, None}, GridLinesStyle -> Directive[Red, Dotted], PlotRange -> All];
PairingViews[model_, delta_, core_Integer] := Module[{xyz, y0, line},
  xyz = MapThread[Append, {model["Coordinates"], Abs[delta]}];
  y0 = SiteCoordinates[core, model["Ny"]][[2]];
  line = Table[{x, Abs[delta[[SiteNumber[{x, y0}, model["Ny"]]]]]}, {x, model["Nx"]}];
  <|"Line" -> ListLinePlot[line, Frame -> True, FrameLabel -> {"x (lattice sites)", "|Delta|"}, PlotRange -> All],
    "Map" -> ListDensityPlot[xyz, InterpolationOrder -> 0, PlotLegends -> Automatic,
      FrameLabel -> {"x (lattice sites)", "y (lattice sites)"}, PlotRange -> All],
    "Surface" -> ListPlot3D[xyz, InterpolationOrder -> 1, Mesh -> All,
      AxesLabel -> {"x (lattice sites)", "y (lattice sites)", "|Delta|"}, PlotRange -> All]|>
];
SpectrumView[h_, references_: {}] := ListPlot[Sort[Eigenvalues[h]], Frame -> True,
  FrameLabel -> {"Sorted-state index", "Energy (model units)"}, PlotRange -> All,
  GridLines -> {None, references}, GridLinesStyle -> Directive[Green, Dashed]];
LDOSMap[model_, sites_, rho_, energyIndex_Integer] := Module[{xyz},
  bdgRequire[Sort[sites] === Range[model["NN"]], "A full spatial LDOS map needs all lattice sites."];
  xyz = MapThread[Append, {model["Coordinates"][[sites]], rho[[energyIndex]]}];
  ListDensityPlot[xyz, InterpolationOrder -> 0, PlotLegends -> Automatic,
    FrameLabel -> {"x (lattice sites)", "y (lattice sites)"}, PlotRange -> All]
];

(* SECTION: export *)
ClearAll[ExportLDOS];
ExportLDOS[folder_, result_] := Module[{rows, history},
  bdgRequire[! DirectoryQ[folder], "Choose a new output folder to preserve existing results."];
  CreateDirectory[folder, CreateIntermediateDirectories -> True];
  rows = Flatten[Table[{result["Energy"][[e]], result["Sites"][[j]], result["LDOS"][[e, j]]},
    {e, Length[result["Energy"]]}, {j, Length[result["Sites"]]}], 1];
  Export[FileNameJoin[{folder, "ldos.csv"}], Prepend[rows, {"energy", "site", "spin_averaged_ldos"}], "CSV"];
  Export[FileNameJoin[{folder, "case.wxf"}], result, "WXF"];
  If[KeyExistsQ[result, "IterationHistory"],
    history = Prepend[result["IterationHistory"], {"iteration", "gap_residual", "density", "onsite_shift"}];
    Export[FileNameJoin[{folder, "iterations.csv"}], history, "CSV"]];
  folder
];

(* SECTION: run *)
ClearAll[RunLDOSCase];
Options[RunLDOSCase] = {"Nx" -> 20, "Ny" -> 20, "Ns" -> 10,
  "Wmax" -> .10, "Intervals" -> 400, "Eta" -> .005,
  "Sites" -> {190, 180, 1}, "ImpurityKind" -> "none", "Strength" -> 0., "CoreSite" -> 190,
  "CleanBaselineConfirmed" -> False, "PairingMode" -> "Imported",
  "SeedVortex" -> False, "CoherenceSites" -> .3/Pi, "PhaseMode" -> "Winding",
  "Coupling" -> -1.25, "kT" -> .001, "TargetDensity" -> 1., "MaxIterations" -> 100,
  "Tolerance" -> 10^-7, "DensityTolerance" -> 10^-6};
RunLDOSCase[folder_, OptionsPattern[]] := Catch[Module[
  {model, phases, delta, factor, mode, settings, relaxed, shift = 0., energy, sites,
   rho, result, kind = OptionValue["ImpurityKind"], core = OptionValue["CoreSite"]},
  bdgRequire[OptionValue["Wmax"] > 0 && IntegerQ[OptionValue["Intervals"]] &&
    OptionValue["Intervals"] > 0 && IntegerQ[OptionValue["Ns"]], "Invalid sampling parameters."];
  bdgRequire[kind === "none" || TrueQ[OptionValue["CleanBaselineConfirmed"]],
    "Confirm that input blocks are an impurity-free baseline before adding a potential."];
  model = LoadBdGModel[folder, OptionValue["Nx"], OptionValue["Ny"]];
  bdgRequire[1 <= core <= model["NN"], "Core site is outside the lattice."];
  model = AddImpurity[model, kind, OptionValue["Strength"], core];
  phases = BoundaryPhases[model]; delta = model["Delta"];
  factor = ConstantArray[1., model["NN"]];
  If[TrueQ[OptionValue["SeedVortex"]],
    factor = VortexFactor[model, core, OptionValue["CoherenceSites"], OptionValue["PhaseMode"]];
    bdgRequire[Max[Abs[delta]] > 0, "A nonzero input pairing amplitude is required to seed a vortex."];
    delta = Max[Abs[delta]] factor];
  mode = OptionValue["PairingMode"];
  bdgRequire[MemberQ[{"Imported", "Relax"}, mode], "PairingMode must be Imported or Relax."];
  bdgRequire[mode =!= "Relax" || TrueQ[OptionValue["SeedVortex"]],
    "Relax mode requires an explicit vortex seed; it applies the saved constrained-vortex update."];
  settings = <|"Coupling" -> OptionValue["Coupling"], "kT" -> OptionValue["kT"],
    "TargetDensity" -> OptionValue["TargetDensity"], "MaxIterations" -> OptionValue["MaxIterations"],
    "Tolerance" -> OptionValue["Tolerance"], "DensityTolerance" -> OptionValue["DensityTolerance"]|>;
  If[mode === "Relax",
    relaxed = RelaxPairing[model, delta, factor, phases, settings];
    bdgRequire[TrueQ[relaxed["Converged"]],
      "Pairing did not converge. Inspect RelaxPairing history; no converged LDOS result was produced."];
    delta = relaxed["Delta"]; shift = relaxed["Shift"]];
  sites = Replace[OptionValue["Sites"], All -> Range[model["NN"]]];
  energy = N[Subdivide[-OptionValue["Wmax"], OptionValue["Wmax"], OptionValue["Intervals"]]];
  rho = CalculateLDOS[model, delta, phases, sites, energy, OptionValue["Eta"], OptionValue["Ns"], shift];
  result = <|"Model" -> model, "Delta" -> delta, "Shift" -> shift, "Sites" -> sites,
    "Energy" -> energy, "LDOS" -> rho, "Ns" -> OptionValue["Ns"], "Eta" -> OptionValue["Eta"],
    "PairingMode" -> mode, "SeedVortex" -> OptionValue["SeedVortex"],
    "PhaseMode" -> OptionValue["PhaseMode"], "CoherenceSites" -> OptionValue["CoherenceSites"],
    "IterationSettings" -> settings, "CoreSite" -> core,
    "EnergyUnit" -> "model energy", "Normalization" -> "electron spin average; all eigenvalues"|>;
  If[mode === "Relax", result = Join[result, <|"IterationHistory" -> relaxed["History"]|>]];
  result
]];
