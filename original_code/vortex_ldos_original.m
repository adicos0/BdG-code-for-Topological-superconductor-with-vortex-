(* ::Package:: *)

SetDirectory["/Users/aditya/Downloads/vortex 2"]
kernum = 16;

Ns = 10; (*supercell*)
kkmx = Ns*Nx;
kkmy = Ns*Ny;
wmax = 0.10;
eta = 0.005; (**)
nw = 400; 
ival = Table[i,{i,1,200}];

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


2


Print["reading phasex3"]


pxijmn3 = Table[
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)]], -1., 
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - ypos[j]]], -1., 
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)]], -1., 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)]],Plus[1.], 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - ypos[j]]], Plus[1.], 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]], Plus[1.]]]]]]], 
     
{i,1,NN},{j,1,NN}]; 

For[i = 1, i <= NN, i++,
 For[j = 1, j <= NN, j++, 
   pxij3[i][j] = pxijmn3[[j,i]]]]
   
Clear[pxijmn3]
Phasex = SparseArray[{x_, y_} /; NumericQ[pxij3[x][y]] -> pxij3[x][y], {NN, NN}]
Clear[pxij3]



Print["reading phasey3"]


pyijmn3 = Table[
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)]], -1.,
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)]], Plus[1.],
If[NumericQ[t[xpos[i] - xpos[j]][ypos[i] - (ypos[j] - Ny)]], -1., 
If[NumericQ[t[xpos[i] - xpos[j]][ypos[i] - (ypos[j] + Ny)]], Plus[1.],
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)]], -1.,
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)]], Plus[1.]]]]]]], 
  
   {i,1,NN},{j,1,NN}]; 
    
For[i = 1, i <= NN, i++,
 For[j = 1, j <= NN, j++, 
 
   pyij3[i][j] = pyijmn3[[j,i]]]]
   
Clear[pyijmn3]
Phasey = SparseArray[{x_, y_} /; NumericQ[pyij3[x][y]] -> pyij3[x][y], {NN, NN}]
Clear[pyij3]


Print["Reading H matrix now"]
Huu = Import["Huu.wdx"]; 
Hud = Import["Hud.wdx"]; 
Hdu = Import["Hdu.wdx"]; 
Hdd = Import["Hdd.wdx"]; 
SCs = Import["SCs.wdx"];
Print["Reading H matrix done"]


(*Dimensions[Phasex1]*)


empty=SparseArray[{x_,y_}/;NumericQ[2==2]->2,{NN,NN}];



Clear[kx, ky]
Huuc = Huu*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey]; 
Hudc = Hud*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey]; 
Hduc = Hdu*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey]; 
Hddc = Hdd*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey]; 

hh1=ArrayFlatten[{{Huuc,Hudc},{Hduc,Hddc}}];
SCs1=Take[SCs,{1,NN},{1,NN}];
SCs2=Take[SCs,{NN+1,2NN},{NN+1,2NN}];
SCs1c=SCs1*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey];
SCs2c=SCs1*Exp[I*Re[kx]*Nx*Phasex + I*Re[ky]*Ny*Phasey];
hh2 = ArrayFlatten[{{-Conjugate[Huuc],Conjugate[Hudc]}, {Conjugate[Hduc], -Conjugate[Hddc]}}];

huv=ArrayFlatten[{{empty,SCs1c},{SCs2c,empty}}];
H=ArrayFlatten[{{hh1,huv},{Conjugate[huv],hh2}}];

H == ConjugateTranspose[H];
(*Clear[Huu,SCs,Hopping12];*)



(*NN=Nx*Ny
NN=2*NN*)


Print["starting ldos calc"];
 CloseKernels[];
 LaunchKernels[kernum];
ldos=Total[Flatten[ParallelTable[
(*For[kx = 0, kx <= (Ns - 1.)*2*(Pi/kkmx), kx += 2*(Pi/kkmx),
 For[ky = 0, ky <= (Ns - 1.)*2*(Pi/kkmy), ky += 2*(Pi/kkmy),*)
  ES = Eigensystem[H];
   elist = Table[ES[[1,l]], {l, 1, 4*NN}];
 
   ulist=Table[Abs[ES[[2,l,ival]]]^2,{l,1,4*NN}(*,{i,1,(nsite-1)*Nband+Nband}*)];
vlist=Table[Abs[ES[[2,l,ival+2*NN]]]^2,{l,1,4*NN}(*,{i,1,(nsite-1)*Nband+Nband}*)];


ldossc=Table[Im[Total[Table[(ulist[[l,;;]]) (1/(-wmax+wmax*2*w/nw-elist[[l]]+I*eta))+
(vlist[[l,;;]]) (1/(-wmax+2*wmax*w/nw+elist[[l]]+I*eta)),{l,1,4*NN}]]],{w,0,nw}]; ldossc,
{kx,0,(Ns - 1.)*2*(Pi/kkmx),2*(Pi/kkmx)},{ky,0,(Ns - 1.)*2*(Pi/kkmy),2*(Pi/kkmy)}],1]];
      (*ldos = ldos + ldossc;*) 
      Print["done"];
      
      
      
      (*Print["kx=", kx, "   ky=", ky, "  nsx=", (kx/(2*Pi))*kkmx + 1.];*) 
    
  (*Export["ldos_sup.dat",(* ldos*)Table[{-wmax + wmax*2*((n-1)/nw),(-Pi^(-1))*ldos[[n,i]]/Ns^2}, {n, 1, nw+1}, {i, 1, Length[ival]}]]

Print["done"];*)
fl = StringJoin["wmax",ToString[wmax],"_Ny",ToString[Ny],"_Nx",ToString[Nx],
"_eta",ToString[\[Eta]],"_Ns",ToString[Ns]];
DeleteFile[StringJoin[fl,"_ldos_sup_cdw.dat"]];
str1 = OpenAppend[StringJoin[fl,"_ldos_sup_cdw.dat"], FormatType -> StandardForm];
For[n=1,n<=nw+1,n++,
For[i=1,i<=Length[ival],i++,
Write[str1,-wmax + wmax*2*((n-1)/nw),"  ", (-1/Pi)*ldos[[n,i]]/Ns^2];
]
];
Close[str1];

(*Exit[];*)





Lds= Import["wmax0.1_Ny20_Nx20_eta\[Eta]_Ns10_ldos_sup_cdw.dat","Table"];
ListPlot[Lds, Joined->True, PlotRange->All,ImageSize->Large, Frame->True, FrameLabel->
{Energy[eV],LDOS},GridLines -> {{{0.04, Red}, {-0.04, Red}},None }, GridLinesStyle-> Dotted,  Epilog -> {Text["(\[CapitalDelta]=0.04)", {0.05, 8}],Text["(\[CapitalDelta]=-0.04)", {-0.05, 8}]}]




(*ListPlot[Lds[[12500;;27500]],Joined->True, PlotRange->[0.7,0.2]]*)


(*m = 2;
Table[Lds[[i]],{i,m,Length[Lds],2}];*)


(*LDos1 = Import["wmax0.15_Ny20_Nx20_eta\[Eta]_Ns20_ldos_sup_cdw.dat","Table"];
En = Import["Energy.dat","Table"];
ListDensityPlot[LDos1,  ColorFunction -> "TemperatureMap", 
  PlotRange -> Full, PlotLegends -> Automatic, 
  FrameLabel -> {"Position", "Energy"}, ImageSize -> Large]*)


ListPlot[ {Table[(-1/Pi)*ldos[[n,184]],{n,1,nw}](*,Table[(-1/Pi)*ldos[[n,1]],{n,1,nw}]*)}, Joined->True]
