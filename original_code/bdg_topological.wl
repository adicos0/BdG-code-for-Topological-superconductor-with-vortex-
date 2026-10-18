(* Original input 1; BdG topological.nb, line 23 *)
(*H1x*)
SetDirectory["\</Users/aditya/Downloads/vortex\>"]
m = 0.50; 
lam = 1;(* next nearest neighbour hopping (diagonal)*)
t[0][0] = 3/2*m;
t[-1][0] = -m;
t[0][-1] = -m;
t[1][0] = -m;
t[0][1] = -m;
t[-1][-1] = m/4;
t[1][-1] = m/4;
t[-1][1] = m/4;
t[1][1] = m/4;
(*SOC spinorbit*)
l[-1][0] = lam(1-I);
l[0][-1] = lam(1-I);
l[1][0] = lam(1-I);
l[0][1] = lam(1-I);
Nx = 20;
Ny = 20;
i = 0;

For[ix = 1, ix <= Nx, ix++,
 For[iy = 1, iy <= Ny, iy++, i++;
  xpos[i] = ix;
  ypos[i] = iy;]]
NN = i;

Vimp = 0;
empty= SparseArray[{x_, y_} -> 0, {NN, NN}];
(*Hopping Matrix*)
tijmn = Table[
   If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)]], t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)],
    If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j])]], t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j])],
     If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)]], t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)],
      If[NumericQ[t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] - Ny)]], t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] - Ny)],
       If[NumericQ[t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j])]], t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j])],
        If[NumericQ[t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] + Ny)]], t[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] + Ny)],
         If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)]], t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)],
          If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j])]], t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j])],
           If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]], t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]]]]]]]]]],
   {i, 1, NN}, {j, 1, NN}];
For[i = 1, i <= NN, i++,
  For[j = 1, j <= NN, j++,
   tij[i][j] = tijmn[[j, i]]
   ]];
Clear[tijmn]
H1= SparseArray[{x_, y_} /; NumericQ[tij[x][y]] -> tij[x][y], {NN, NN}];
tijmn = Table[
   If[NumericQ[l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)]], l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)],
    If[NumericQ[l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j])]], l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j])],
     If[NumericQ[l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)]], l[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)],
      If[NumericQ[l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] - Ny)]], l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] - Ny)],
       If[NumericQ[l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j])]], l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j])],
        If[NumericQ[l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] + Ny)]], l[xpos[i] - (xpos[j])][ypos[i] - (ypos[j] + Ny)],
         If[NumericQ[l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)]], l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)],
          If[NumericQ[l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j])]], l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j])],
           If[NumericQ[l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]], l[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]]]]]]]]]],
   {i, 1, NN}, {j, 1, NN}];
For[i = 1, i <= NN, i++,
  For[j = 1, j <= NN, j++,
   tij[i][j] = tijmn[[j, i]]
   ]];
Clear[tijmn]
soc= SparseArray[{x_, y_} /; NumericQ[tij[x][y]] -> tij[x][y], {NN, NN}];
Clear[mu]
Chempot1[mu_]:= SparseArray[{x_, x_} ->  -mu, {NN, NN}];
(*Chempot2[mu_]:= SparseArray[{y_, y_} ->  -mu, {NN, NN}];*)
(*Chempot[mu]=ArrayFlatten[{{Chempot1[mu], empty}, {empty,Chempot1[mu]}}];*)
Hop1x[mu_]:= ArrayFlatten[{{H1-Chempot1[mu],soc}, {ConjugateTranspose[soc], -H1-Chempot1[mu]}}];(*adding pauli matrix \[Sigma]x*)
(*HermitianMatrixQ[Hop1x[0]];
MatrixForm[Hop1x[]];*)


(* Original input 2; BdG topological.nb, line 922 *)
(*------onsite------*)
Clear[s,st,stab,sij,SCs,dmn,dmm]

js={};
For[i=1,i<=NN,i++,jisel={};
For[j=1,j<=NN,j++,
If[xpos[i]==xpos[j],
If[ypos[i]==ypos[j],AppendTo[jisel,j]]]];AppendTo[js,jisel]]



SCs1=SparseArray[{x_,y_}/;Not[TrueQ[sij[x][y]==Null]]\[Rule]0,{NN,NN}];
SCs2[st_]:= SparseArray[{Band[{1,1}]\[Rule] Flatten[st]}, {NN, NN}];
SCs3=SparseArray[{x_,y_}/;Not[TrueQ[sij[x][y]==Null]]\[Rule]0,{NN,NN}];
SCs[st_]:= ArrayFlatten[{{SCs1, SCs2[st]}, {SCs2[st], SCs1}}];
Clear[sij]
s[i_,j_]:=Take[st[[i,j]]];

(* Original input 3; BdG topological.nb, line 1049 *)


f[x_]:=1./(Exp[x/kT]+1.);
kT = 0.001 (*eV*);
f[x_] :=(*UnitStep[x]*)(*1./(Exp[x/kT]+1.);*)If[Re[x] > 0.02, 0.0,
   If[Re[x] < -0.02, 1.0,
    1./(Exp[Re[x]/kT] + 1.)]];

(*mtuss={};*)
(*For[k=0,k\[LessEqual]0.2,k=k+0.1,
Print[StringJoin["\<calculte..\>",ToString[k]]];*)
(*flag1 = 1;*)
(*m = 0.5*0;*)
(*Us = 0.0;*)
vs = -1.45;(*s-wave pairing interaction*)
(*vd = -0.2*0; d-wave pairing interaction*)
mu = 0;

nexp = 1.2; (*expected electro density in your system*)
(*ntot = 1;
nut=Table[(ntot/2.)+((m*(-1)^(xpos[i]+ypos[i]))/2),{i,1,NN}];
ndt=Table[(ntot/2.)-((m*(-1)^(xpos[i]+ypos[i]))/2),{i,1,NN}];*)

Clear[v,as,b,bm];

(*dsignn={}; (*initialisation of d-wave sign*)
For[i=1,i<=NN,i++,jisel={};
For[j=1,j<=NN,j++,
If[xpos[i]==xpos[j],
If[ypos[i]==ypos[j]+1,AppendTo[jisel,-1],
If[ypos[i]==ypos[j]-1,AppendTo[jisel,-1],
If[ypos[i]==ypos[j]+(Ny-1),AppendTo[jisel,-1],
If[ypos[i]==ypos[j]-(Ny-1),AppendTo[jisel,-1]]]]],
If[ypos[i]==ypos[j],
If[xpos[i]==xpos[j]+1,AppendTo[jisel,1],
If[xpos[i]==xpos[j]-1,AppendTo[jisel,1],
If[xpos[i]==xpos[j]+(Nx-1),AppendTo[jisel,1],
If[xpos[i]==xpos[j]-(Nx-1),AppendTo[jisel,1]]]]]]]];
AppendTo[dsignn,jisel]];
spmtnn=Table[dsignn[[i,j]]*dijnn,{i,1,NN},{j,1,4}];
dijnn=0.3;*)

st=Table[ds,{i,1,NN},{j,1}]; (*initialisation of s-wave*)
ds=-0.01; (*Gap size*)
(*AbsoluteTiming[*)
(*MatrixForm[SCs[st]]*)

(* Original input 4; BdG topological.nb, line 1398 *)
(*
 Huu = Ho (*+imp*)(**+Hartree*);
  Hvv = -ConjugateTranspose[Ho (*-imp*)(*+Hartree2*)];
  Huv= SCs[st]; (*Superconducting part \[CapitalDelta]*)
  Hvu=ConjugateTranspose[Huv];
  H = ArrayFlatten[{{Huu, Huv}, {Hvu, Hvv}}]; 
MatrixForm[H];*)

(* Original input 5; BdG topological.nb, line 1441 *)
(*calculation mean fields via iteration*)(*Main Hamiltonian BDG*)
iter = 100;
  (*em =  SparseArray[{x_, x_} -> 0, {NN, NN}];*)
For[a = 1, a <= iter, a++,
  (*mtus={};*)
  Huu = H1-Chempot1[mu] (*+imp*)(**+Hartree*);
  Hdd = H1-Chempot1[mu];
Hud = soc;
Hdu = Conjugate[soc];
hh1 = ArrayFlatten[{{Huu,Hud}, {Hdu, Hdd}}];
hh2 = ArrayFlatten[{{-Conjugate[Huu],Conjugate[Hud]}, {Conjugate[Hdu], -Conjugate[Hdd]}}];
    Huv= SCs[st]; (*Superconducting part \[CapitalDelta]*)
  Hvu=Conjugate[Huv];
  H = ArrayFlatten[{{hh1, Huv}, {Hvu, hh2}}]; (*Joins all the matrices*)
 (* Print[HermitianMatrixQ[H]];*)
  (*Print[MatrixForm[H]];*)
  
  ES = Eigensystem[H];
 
  fermi = Map[f, ES[[1]]];
  (*tab = Table[1,{i,Length[fermi]}];*)
  nutcalc = Table[((Abs[ES[[2]][[;; , i]]])^2.) . fermi, {i, 1, NN}];
  
  
  ndtcalc = Table[((Abs[ES[[2]][[;; , i + 3NN]]])^2.) . (1. - fermi), {i, 1, NN}];
  
  
  ni = nutcalc + ndtcalc;
  mi = nutcalc - ndtcalc;
  nt = Total[ni]/NN;
  mt = Total[Abs[mi]]/NN;
  
  
  
  (*Print["\<till here\>"];*)
  stcalc=Table[vs*((Inner[Times,ES[[2]][[;;,i]],Conjugate[ES[[2]][[;;,j+
  3NN]]],List])+(Inner[Times,ES[[2]][[;;,i+NN]],Conjugate[ES[[2]][[;;,i+
  2NN]]],List])).fermi,{i,1,NN},{j,js[[i]]}];
  (*
vsite = 378;
  l=0.3;(*penetration depth*)
  zetanot=l/Pi;(*coherence length*)
    stcalc=Table[
  If[xpos[i]==xpos[vsite],
  phi=Pi/2,
  phi=ArcTan[(Abs[ypos[i]-ypos[vsite]])/(Abs[xpos[i]-xpos[vsite]])];
  ];
  vs*((Inner[Times,ES[[2]][[;;,i]],Conjugate[ES[[2]][[;;,j+
  3NN]]],List])+(Inner[Times,ES[[2]][[;;,i+NN]],Conjugate[ES[[2]][[;;,i+
  2NN]]],List])).fermi*
  Tanh[Sqrt[(xpos[i]-xpos[vsite])^2+(ypos[i]-ypos[vsite])^2]/zetanot]*
  Exp[I*phi]
,{i,1,NN},{j,js[[i]]}];
*)

  (*spmtcalcnn=Table[vd*((Inner[Times,ES[[2]][[;;,i]],ES[[2]][[;;,j+
  NN]],List])+(Inner[Times,ES[[2]][[;;,j]],ES[[2]][[;;,i+
  NN]],List])).fermi,{i,1,NN},{j,jspmnn[[i]]}];*)
  
  
  sta=Total[Flatten[stcalc]]/NN;(*
  stann=Total[Total[Abs[spmtcalcnn]]]/NN*0.25;*)
  
  
  
  Print[a , "\< mu=\>", mu, "\<  mt=\>", mt,(*"\< Us=\>", Us,*)"\<  n=\>", nt ,"\<  sup=\>", Abs[ sta(*"\<  supnn=\>",  stann*)]];
  mu = mu - 0.1*(nexp - nt);
  (*nut=0.5*nutcalc+0.5*nut;
  ndt=0.5*ndtcalc+0.5*ndt;*)
  st=0.5*stcalc+0.5*st;
  (*spmtnn=0.5*spmtcalcnn+0.5*spmtnn;*)
  ];
(*AppendTo[mtuss,{vd,stann}];
];*)


Export["\<Huu.wdx\>",Huu];
Export["\<Hud.wdx\>",Hud];
Export["\<Hdu.wdx\>",Hdu];
Export["\<Hdd.wdx\>",Hdd];
Export["\<SCs.wdx\>",SCs[st]];

(* Original input 6; BdG topological.nb, line 4032 *)
SortBy[ES[[1]]];

(* Original input 7; BdG topological.nb, line 4073 *)
ListPlot[Sort[ES[[1]]],PlotRange\[Rule]{-0.6,0.6} ,Frame\[Rule]True,PlotLabel\[Rule]"\<\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(o\)]\)=\[PlusMinus]0.7, M=0, m=0.5, \[Lambda]=1, <n>=1, \!\(\*SubscriptBox[\(\[CapitalDelta]\), \(con\)]\)=0.064, \!\(\*SubscriptBox[\(\[Mu]\), \(con\)]\)=-0.46 \>",FrameLabel\[Rule]{n,Energy},LabelStyle\[Rule]"\<DisplayFormula\>",GridLines\[Rule]{None,{-0.05,0.05}},GridLinesStyle\[Rule]Directive[Green, Dashed], Frame\[Rule]True, FrameStyle\[Rule]Gray,ImageSize\[Rule]Large]
