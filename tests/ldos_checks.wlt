(* Run in Mathematica: TestReport["/absolute/path/to/repo/tests/ldos_checks.wlt"] *)
Get[FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]], "scripts", "vortex_ldos.wl"}]];
nx = 3; ny = 3; n = nx ny;
z = ConstantArray[0., {n, n}];
a = DiagonalMatrix[ConstantArray[.2, n]];
a[[1, 7]] = -.3; a[[7, 1]] = -.3;
soc = .04 (1. - I) IdentityMatrix[n];
model = <|"Nx" -> nx, "Ny" -> ny, "NN" -> n,
 "Coordinates" -> (SiteCoordinates[#, ny] & /@ Range[n]),
 "Huu" -> a, "Hdd" -> a, "Hud" -> soc, "Hdu" -> ConjugateTranspose[soc],
 "Delta" -> ConstantArray[.07 + .02 I, n]|>;
p = BoundaryPhases[model];
k = {.13, -.17};
h = Catch[BdGMatrix[model, model["Delta"], k, p]];
hm = Catch[BdGMatrix[model, model["Delta"], -k, p]];
s = DiagonalMatrix[Join[ConstantArray[1., n], ConstantArray[-1., n]]];
z2 = ConstantArray[0., {2 n, 2 n}];
c = ArrayFlatten[{{z2, s}, {s, z2}}];
VerificationTest[SiteNumber[SiteCoordinates[190, 20], 20], 190, TestID -> "site-index round trip"]
VerificationTest[p[[1, 7, 1]], -1, TestID -> "original boundary orientation"]
VerificationTest[bdgNorm[h - ConjugateTranspose[h]] < 10^-12, True, TestID -> "Hermiticity"]
VerificationTest[bdgNorm[c . Conjugate[h] . c + hm] < 10^-12, True, TestID -> "opposite-k particle-hole relation"]
VerificationTest[Catch[AddImpurity[model, "scalar", 0., 5]]["Huu"] == model["Huu"], True, TestID -> "zero impurity"]
mag = Catch[AddImpurity[model, "magnetic", .5, 5]];
VerificationTest[{mag["Huu"][[5, 5]] - a[[5, 5]], mag["Hdd"][[5, 5]] - a[[5, 5]]}, {.5, -.5}, TestID -> "magnetic signs"]
factor = Catch[VortexFactor[model, 5, .3]];
VerificationTest[factor[[5]], 0., TestID -> "zero amplitude at vortex core"]
VerificationTest[Max[Abs[factor]] <= 1., True, TestID -> "bounded vortex envelope"]
{ev, vec} = Eigensystem[h];
weights = (Abs[vec[[All, {1, 5}]]]^2 + Abs[vec[[All, {1, 5} + n]]]^2)/2.;
VerificationTest[Max[Abs[Total[weights] - {1., 1.}]] < 10^-12, True, TestID -> "all-energy spectral sum rule"]
rho = Catch[CalculateLDOS[model, model["Delta"], p, {1, 5}, {-.1, 0., .1}, .005, 1]];
VerificationTest[Dimensions[rho], {3, 2}, TestID -> "energy-by-site output"]
VerificationTest[Min[Flatten[rho]] >= 0., True, TestID -> "nonnegative LDOS"]
