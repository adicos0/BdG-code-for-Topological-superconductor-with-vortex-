(* Load functions only. Evaluate the examples below individually after editing paths.
   No 20x20 eigensystem or parameter sweep runs on Get of this file. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "vortex_ldos.wl"}]];

(* Example A: use saved pairing exactly as imported.
result = RunLDOSCase[$BdGProjectRoot, "Ns" -> 10, "Sites" -> {190, 180, 1}];
If[! FailureQ[result],
  LDOSCurves[result["Energy"], result["LDOS"], result["Sites"]]
]
*)

(* Example B: only after identifying an impurity-free input folder.
   Ns=20 is the documented vortex/impurity sampling choice, not a default proof
   that all other historical inputs are known. Repeat FROM THIS BASELINE for
   scalar/magnetic V=0.5, 5, 50; line profiles also use V=1.
result = RunLDOSCase["/absolute/path/to/verified/clean/input",
  "CleanBaselineConfirmed" -> True, "ImpurityKind" -> "scalar", "Strength" -> .5,
  "CoreSite" -> 190, "Sites" -> {190, 180, 1}, "Ns" -> 20,
  "PairingMode" -> "Relax", "SeedVortex" -> True,
  "PhaseMode" -> "OriginalAbsoluteAngle", "TargetDensity" -> 1.1];
If[! FailureQ[result],
  curves = LDOSCurves[result["Energy"], result["LDOS"], result["Sites"], {-.05, .05}];
  views = PairingViews[result["Model"], result["Delta"], result["CoreSite"]];
  phases = BoundaryPhases[result["Model"]];
  h = Catch[BdGMatrix[result["Model"], result["Delta"], {0., 0.}, phases, result["Shift"]]];
  spectrum = SpectrumView[h, {-.05, .05}];
  ExportLDOS["/absolute/path/to/new/output/case_scalar_0p5", result]
]
*)

(* Example C: full spatial LDOS is different from a pairing-amplitude map.
result = RunLDOSCase[$BdGProjectRoot, "Sites" -> All];
If[! FailureQ[result],
  zeroRow = First[Ordering[Abs[result["Energy"]], 1]];
  LDOSMap[result["Model"], result["Sites"], result["LDOS"], zeroRow]
]
*)

(* Example D: after separately converging three cases with strengths 0.5, 1, 5,
   collect their pairing cuts. Each result must use the same clean baseline,
   core position and conventions; do not reuse an already perturbed normal block.
   results = {case05, case1, case5};
If[And @@ (! FailureQ[#] & /@ results),
  GraphicsRow[(PairingViews[#["Model"], #["Delta"], #["CoreSite"]]["Line"] &) /@ results]
]
*)
