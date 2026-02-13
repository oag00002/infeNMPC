#compared with orig, M_0 and x_0 are defined directly as unbounded vars.


param nfe;  
param ncp;
param time; # horizon length

set fe := 1..nfe ;  
set cp := 1..ncp ; 

param omega{cp,cp};   param h{fe}; 

param inf:=10000000;
param Mup:=500;  param Mlow:=20;    #bounds on M
param Tup:=335;  param Tlow:=280;   #bounds on T,Tbottom,Ttop
param Vup:=12000; param Vlow:=0;    #bounds on V
#param Lup:=12000; param Llow:=0;    #bounds on L
param L1up:=380; param L1low:=80;   #bounds on L[1,q,c]
param Recup:=18; param Reclow:=8;   #bounds on Rec
param xup:=1; param xlow:=0;        #bounds on x, x_0, x_nom, xout, yout
param Qrup:=140; param Qrlow:=40;   #bounds on Qr
param Qclow:=40; param Qcup:=inf;                 #bound on Qc
param Dup:=500; param Dlow:=300;
param Bup:=380; param Blow:=200;
##param Bup:=380; param Blow:=50;
##param Dup:=500; param Dlow:=200;
##param Bup:=500; param Blow:=50;
param Psteamup:=400; param Psteamlow:=300;
##param Psteamup:=365; param Psteamlow:=300;

param Tsup:=360; param Tslow:=310;
param Fsup:=100; param Fslow:=0;
param Pchestup:=30; param Pchestlow:=10;
param Toutup:=430; param Toutlow:=290;
param Wup:=500000; param Wlow:=0;        #bounds on Wrdb, Wrebl, Wrebv
param rorebup:=15; param roreblow:=0.1;  #bounds on roreb, rowreb
param rorebmup:=550; param rorebmlow:=0; #bounds on rorebm
param rebup:=100; param reblow:=0;       #bounds on Mvreb, Mlreb, Vvreb, Vlreb

param Refluxup:=inf; param Refluxlow:=-inf;

param Ntray          ;
param feedtray;
set tray := 1..Ntray  ;
#let Ntray := 158;
#let feedtray := 43;
#let nfe       := 18; 
#let ncp       := 3;
let ncp:=2;

param kQr>=0,<=3;     #cte for Efficiency of Reboiler
param Eff;
param kwreb >=1E-6,<=1E6;
param kFlow >=1E-6,<=1E6;
param Uheat >=0,<=10;		#Coefficient heat transfer


# define the hold in each tray
##var   M{tray,fe,cp} >=Mlow,<=Mup ;
##var   Mdot{tray, fe,cp};
##var M_0{tray, fe}>=Mlow,<=Mup;
##var M_nom{tray}>=Mlow,<=Mup;

var   M{tray,fe,cp};
var Meps{tray,fe,cp}>=0;

Mupper{i in tray, j in fe, k in cp}: M[i,j,k]<=Mup+Meps[i,j,k];
Mlower{i in tray, j in fe, k in cp}: M[i,j,k]>=Mlow-Meps[i,j,k];

var M_0{tray, fe}>=Mlow,<=Mup;
#var M0eps{tray,fe}>=0;
#M0upper{i in tray, j in fe}: M_0[i,j]<=Mup+M0eps[i,j];
#M0lower{i in tray, j in fe}: M_0[i,j]>=Mlow-M0eps[i,j];

var M_nom{tray}>=Mlow,<=Mup;
#Mnomlower{i in tray}: M_nom[i]>=Mlow-M0eps[i,1];
#Mnomupper{i in tray}: M_nom[i]<=Mup+M0eps[i,1];

# define the temperature in each tray
var T{tray, fe, cp};
var Teps{tray,fe,cp}>=0;
Tupper{i in tray, j in fe, k in cp}:T[i,j,k]<=Tup+Teps[i,j,k];
Tlower{i in tray, j in fe, k in cp}:T[i,j,k]>=Tlow-Teps[i,j,k];
##var Tbottom{fe, cp} >=Tlow, <=Tup;
##var Ttop{fe, cp} >=Tlow, <=Tup;

###var Tdot{tray,fe,cp};


param Tref;

param Ttop_des; param Tbottom_des;

param feednmpc{fe, cp};

param feedup:=750;
param feedlow:=750;

#parameters to MHE

param eD{fe}:=0; param eReflux{fe}:=0; param eB{fe}:=0;
param eEff{fe,cp}:=0; param eU{fe,cp}:=0; 
param eL{fe,cp}:=0; param eWeir{fe,cp}:=0; 

#parameters to simulation
param xtop{fe,cp}; param xbottom{fe,cp}; param Mtop{fe,cp}; param Mbottom{fe,cp}; param Mset{tray}; param xset{tray}; param Refluxset{fe}; param feedset{fe,cp}; param Psteamset{fe}; param Dset{fe}; param Qcset{fe,cp}; param Bset{fe};
param Ttopset{fe,cp}; param Tbottomset{fe,cp}; 

# define the pressure
param p{tray}>=15, <= 20;
#this max restriction to P was added to avoid negative sqrt
#param p{tray}>=0, <= 3000000;
param prect; param pstrip; 

param CapAm; param CapBm; 
param CapAn; param CapBn; 


# define the liquid and vapour flowrate in each tray
var V{tray,fe,cp};
var Veps{tray,fe,cp}>=0;

Vupper{i in tray, j in fe, k in cp}:V[i,j,k]<=Vup+Veps[i,j,k];
Vlower{i in tray, j in fe, k in cp}:V[i,j,k]>=Vlow-Veps[i,j,k];

param Lup{tray};#:=12000;
param Llow{tray};#:=0;
var L{tray,fe,cp};
var Leps{tray,fe,cp}>=0;
Lupper{i in tray, j in fe, k in cp}:L[i,j,k]<=Lup[i]+Leps[i,j,k];
Llower{i in tray, j in fe, k in cp}:L[i,j,k]>=Llow[i]-Leps[i,j,k];

var Rec{fe,cp};
var Receps{fe,cp}>=0;
Recupper{j in fe, k in cp}:Rec[j,k]<=Recup+Receps[j,k];
Reclower{j in fe, k in cp}:Rec[j,k]>=Reclow-Receps[j,k];

var Reflux{fe};


# define the molar fraction in liquid and vapour
##var x{tray,fe,cp}>=xlow ,<=xup;
##var xdot{tray, fe, cp};
##var x_0{tray, fe} >=xlow , <=xup;
##var x_nom{tray} >=xlow , <=xup;

var x{tray,fe,cp};
var xeps{tray,fe,cp}>=0;

xupper{i in tray, j in fe, k in cp}: x[i,j,k]<=xup+xeps[i,j,k];
xlower{i in tray, j in fe, k in cp}: x[i,j,k]>=xlow-xeps[i,j,k];

var x_0{tray, fe} >=xlow , <=xup;
#var x0eps{tray,fe}>=0;
#x0upper{i in tray, j in fe}: x_0[i,j]<=1+x0eps[i,j];
#x0lower{i in tray, j in fe}: x_0[i,j]>=-x0eps[i,j];

var x_nom{tray} >=xlow , <=xup;
#xnomupper{i in tray}: x_nom[i]<=1+x0eps[i,1];
#xnomlower{i in tray}: x_nom[i]>=-x0eps[i,1];

var y{tray,fe,cp};
var yeps{tray,fe,cp}>=0;
yupper{i in tray,j in fe, k in cp}:y[i,j,k]<=xup+yeps[i,j,k];
ylower{i in tray,j in fe, k in cp}:y[i,j,k]>=xlow-yeps[i,j,k];
#>=xlow, <=xup;

# define the enthalpy 
#var hl{tray, fe,cp}>=5000;
#var hv{tray,fe,cp}>=15000 ;

# define the 
var Qr{fe,cp};
var Qreps{fe,cp}>=0;
Qrupper{j in fe, k in cp}:Qr[j,k]<=Qrup+Qreps[j,k];
Qrlower{j in fe, k in cp}:Qr[j,k]>=Qrlow-Qreps[j,k];

var Qc{fe,cp};
var Qceps{fe,cp}>=0;
Qclower{j in fe, k in cp}:Qc[j,k]>=Qclow-Qceps[j,k];

var D{fe};
var Deps{fe}>=0;
Dupper{j in fe}:D[j]<=Dup+Deps[j];
Dlower{j in fe}:D[j]>=Dlow-Deps[j];

var B{fe}>=Blow,<=Bup;
#var Beps{fe}>=0;
#Bupper{j in fe}:B[j]<=Bup+Beps[j];
#Blower{j in fe}:B[j]>=Blow-Beps[j];

# define the feed
#var feed{fe,cp}>=feedlow, <=feedup;
#cfeed0{q in fe,c in cp}: feed[q,c]=feednmpc[q,c];
param feed{q in fe, c in cp}=feednmpc[q,c];

param xf; # feed molar fraction
param hf; # feed enthalpy

#Reboiler 
var Psteam{fe}>=Psteamlow,<=Psteamup;
#var Psteameps{fe}>=0;
#Psteamupper{j in fe}:Psteam[j]<=Psteamup+Psteameps[j];
#Psteamlower{j in fe}:Psteam[j]>=Psteamlow-Psteameps[j];#>=Psteamlow,<=Psteamup;	#Vapor pressure of Steam

var Ts{fe,cp};
var Tseps{fe,cp}>=0;
Tsupper{j in fe, k in cp}:Ts[j,k]<=Tsup+Tseps[j,k];
Tslower{j in fe, k in cp}:Ts[j,k]>=Tslow-Tseps[j,k];# >=Tslow , <=Tsup; 	#Temperature steam

var Fs{fe,cp};
var Fseps{fe,cp}>=0;
Fsupper{j in fe, k in cp}:Fs[j,k]<=Fsup+Fseps[j,k];
Fslower{j in fe, k in cp}:Fs[j,k]>=Fslow-Fseps[j,k];#>=Fslow,<=Fsup;		#Flow vapor reboiler
param pout>=15,<=60;			#Pressure reboiler

var Pchest{fe,cp};
var Pchesteps{fe,cp}>=0;
Pchestupper{j in fe, k in cp}:Pchest[j,k]<=Pchestup+Pchesteps[j,k];
Pchestlower{j in fe, k in cp}:Pchest[j,k]>=Pchestlow-Pchesteps[j,k];#>=Pchestlow,<=Pchestup;		#Pressure leaving reboiler

var Tout{fe,cp} >=Toutlow , <=Toutup; 	#Out Temperature of reboiler
var xout{fe,cp}>=xlow,<=xup;		#Pressure reboiler
var yout{fe,cp}>=xlow,<=xup; 		#Out Temperature of reboiler
var Wreb{fe,cp};
var Wrebeps{fe,cp}>=0;
Wrebupper{j in fe, k in cp}:Wreb[j,k]<=Wup+Wrebeps[j,k];
Wreblower{j in fe, k in cp}:Wreb[j,k]>=Wlow-Wrebeps[j,k];#>=Wlow,<=Wup;		#Flow reboiler

var Wrebl{fe,cp};
var Wrebleps{fe,cp}>=0;
Wreblupper{j in fe, k in cp}:Wrebl[j,k]<=Wup+Wrebleps[j,k];
Wrebllower{j in fe, k in cp}:Wrebl[j,k]>=Wlow-Wrebleps[j,k];#>=Wlow,<=Wup;	#Flow reboiler liquid

var Wrebv{fe,cp};
var Wrebveps{fe,cp}>=0;
Wrebvupper{j in fe, k in cp}:Wrebv[j,k]<=Wup+Wrebveps[j,k];
Wrebvlower{j in fe, k in cp}:Wrebv[j,k]>=Wlow-Wrebveps[j,k];#>=Wlow,<=Wup;	#Flow reboiler vapor

var roreb{fe,cp};
var rorebeps{fe,cp}>=0;
rorebupper{j in fe, k in cp}:roreb[j,k]<=rorebup+rorebeps[j,k];
roreblower{j in fe, k in cp}:roreb[j,k]>=roreblow-rorebeps[j,k];#>=roreblow, <=rorebup;		#ro leaving reboiler

var rorebm{fe,cp};
var rorebmeps{fe,cp}>=0;
rorebmupper{j in fe, k in cp}:rorebm[j,k]<=rorebmup+rorebmeps[j,k];
rorebmlower{j in fe, k in cp}:rorebm[j,k]>=rorebmlow-rorebmeps[j,k];#>=rorebmlow, <=rorebmup;		#ro leaving reboiler

var rowreb{fe,cp};
var rowrebeps{fe,cp}>=0;
rowrebupper{j in fe, k in cp}:rowreb[j,k]<=rorebup+rowrebeps[j,k];
rowreblower{j in fe, k in cp}:rowreb[j,k]>=roreblow-rowrebeps[j,k];#>=roreblow, <=rorebup;	#ro reboiler total

var Mvreb{fe,cp};
var Mvrebeps{fe,cp}>=0;
Mvrebupper{j in fe, k in cp}:Mvreb[j,k]<=rebup+Mvrebeps[j,k];
Mvreblower{j in fe, k in cp}:Mvreb[j,k]>=reblow-Mvrebeps[j,k];#>=reblow, <=rebup;		#Hold up vapor

var Mlreb{fe,cp};
var Mlrebeps{fe,cp}>=0;
Mlrebupper{j in fe, k in cp}:Mlreb[j,k]<=rebup+Mlrebeps[j,k];
Mlreblower{j in fe, k in cp}:Mlreb[j,k]>=reblow-Mlrebeps[j,k];#>=reblow, <=rebup;		#Hold up liquid

var Vvreb{fe,cp};
var Vvrebeps{fe,cp}>=0;
Vvrebupper{j in fe, k in cp}:Vvreb[j,k]<=rebup+Vvrebeps[j,k];
Vvreblower{j in fe, k in cp}:Vvreb[j,k]>=reblow-Vvrebeps[j,k];#>=reblow, <=rebup;		#Volume vapor on reboiler

var Vlreb{fe,cp};
var Vlrebeps{fe,cp}>=0;
Vlrebupper{j in fe, k in cp}:Vlreb[j,k]<=rebup+Vlrebeps[j,k];
Vlreblower{j in fe, k in cp}:Vlreb[j,k]>=reblow-Vlrebeps[j,k];#>=reblow, <=rebup;		#Volume liquid on reboiler

param Atotal;				#Total area reboiler

param x1ref{fe};
param xntrayref{fe};
param M1ref{fe};
param Mntrayref{fe};
param Dref{fe};
param Bref{fe};
param Psteamref{fe};
param Refluxref{fe};


# physical property
param Avm; param Bvm; param Cvm; param Dvm; param Evm;
param Avn; param Bvn; param Cvn; param Dvn; param Evn;
param Alm; param Blm; param Clm; param Dlm; 
param Aln; param Bln; param Cln; param Dln;

param h0m;
param h0n;

param MWm; param MWn;


param Adm;param Bdm; param Ndm;
param Adn;param Bdn; param Ndn;

param r; 
param gm; param Tkm; param Pkm;
param gn; param Tkn; param Pkn;

#Temperature
#hrt1{q in fe, c in cp}: 0 = (Tbottom[q,c]+eTbottom[q,c]) - T[1,q,c];
#hrt2{q in fe, c in cp}: 0 = (Ttop[q,c]+eTtop[q,c]) - T[Ntray,q,c];


# pressure
#partial pressure p, see Diehl's P116

var pm{i in tray,q in fe,c in cp}=	 exp(CapAm-CapBm/(T[i,q,c]));

var pn{i in tray,q in fe,c in cp}= exp(CapAn-CapBn/(T[i,q,c]));   

dp{i in tray, q in fe, c in cp} : p[i] = pm[i,q,c]*x[i,q,c]+(1-x[i,q,c])*pn[i,q,c];


# define the equilibrium coefficients

yr{i in tray, q in fe ,c in cp}: y[i,q,c] = ((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* (Eff+eEff[q,c]) + x[i,q,c] ;

#var y{i in tray,q in fe, c in cp}=((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* (Eff+eEff[q,c]) + x[i,q,c] ;
#yr{i in tray, q in fe, c in cp} : xlow <= y[i,q,c] <=xup ;


#restriction monotomic Temperature

#mT{i in tray diff{1,Ntray}, q in fe, c in cp} : T[i+1,q,c]<= T[i,q,c];

# calculate the hkl and hkv

var hl{i in tray, q in fe ,c in cp}= x[i,q,c]*(Alm*(T[i,q,c]-Tref)  + Blm/2*(T[i,q,c]^2-Tref^2) + Clm/3*(T[i,q,c]^3-Tref^3) + Dlm/4*(T[i,q,c]^4-Tref^4) ) + (1-x[i,q,c])*(Aln*(T[i,q,c]-Tref) + Bln/2*(T[i,q,c]^2-Tref^2) + Cln/3*(T[i,q,c]^3-Tref^3) + Dln/4*(T[i,q,c]^4-Tref^4) );

var hlout{q in fe ,c in cp}= xout[q,c]*(Alm*(Tout[q,c]-Tref)  + Blm/2*(Tout[q,c]^2-Tref^2) + Clm/3*(Tout[q,c]^3-Tref^3) + Dlm/4*(Tout[q,c]^4-Tref^4) ) + (1-xout[q,c])*(Aln*(Tout[q,c]-Tref) + Bln/2*(Tout[q,c]^2-Tref^2) + Cln/3*(Tout[q,c]^3-Tref^3) + Dln/4*(Tout[q,c]^4-Tref^4) );


var hv{i in tray, q in fe ,c in cp}= y[i,q,c]*((h0m+Avm*(T[i,q,c]-Tref)  + Bvm/2*(T[i,q,c]^2-Tref^2) + Cvm/3*(T[i,q,c]^3-Tref^3) + Dvm/4*(T[i,q,c]^4-Tref^4) + Evm/5*(T[i,q,c]^5-Tref^5))+(4.18*(r*Tkm/Pkm*(0.083-0.422/(T[i,q,c]/Tkm)^(8/5)+gm*(0.139-0.172/(T[i,q,c]/Tkm)^(21/5)))-T[i,q,c]*(0.672/(T[i,q,c]/Tkm)^(13/5)/Tkm+0.714*gm/(T[i,q,c]/Tkm)^(26/5)/Tkm))*p[i]))+(1-y[i,q,c])*((h0n+Avn*(T[i,q,c]-Tref)  + Bvn/2*(T[i,q,c]^2-Tref^2) + Cvn/3*(T[i,q,c]^3-Tref^3) + Dvn/4*(T[i,q,c]^4-Tref^4) + Evn/5*(T[i,q,c]^5-Tref^5))+(4.18*(r*Tkn/Pkn*(0.083-0.422/(T[i,q,c]/Tkn)^(8/5)+gn*(0.139-0.172/(T[i,q,c]/Tkn)^(21/5)))-T[i,q,c]*(0.672/(T[i,q,c]/Tkn)^(13/5)/Tkn+0.714*gn/(T[i,q,c]/Tkn)^(26/5)/Tkn))*p[i]));

var hvout{q in fe ,c in cp}= yout[q,c]*((h0m+Avm*(Tout[q,c]-Tref)  + Bvm/2*(Tout[q,c]^2-Tref^2) + Cvm/3*(Tout[q,c]^3-Tref^3) + Dvm/4*(Tout[q,c]^4-Tref^4) + Evm/5*(Tout[q,c]^5-Tref^5))+(4.18*(r*Tkm/Pkm*(0.083-0.422/(Tout[q,c]/Tkm)^(8/5)+gm*(0.139-0.172/(Tout[q,c]/Tkm)^(21/5)))-Tout[q,c]*(0.672/(Tout[q,c]/Tkm)^(13/5)/Tkm+0.714*gm/(Tout[q,c]/Tkm)^(26/5)/Tkm))*p[1]))+(1-yout[q,c])*((h0n+Avn*(Tout[q,c]-Tref)  + Bvn/2*(Tout[q,c]^2-Tref^2) + Cvn/3*(Tout[q,c]^3-Tref^3) + Dvn/4*(Tout[q,c]^4-Tref^4) + Evn/5*(Tout[q,c]^5-Tref^5))+(4.18*(r*Tkn/Pkn*(0.083-0.422/(Tout[q,c]/Tkn)^(8/5)+gn*(0.139-0.172/(Tout[q,c]/Tkn)^(21/5)))-Tout[q,c]*(0.672/(Tout[q,c]/Tkn)^(13/5)/Tkn+0.714*gn/(Tout[q,c]/Tkn)^(26/5)/Tkn))*pout));


# mass balance
cm1_init_constr{i in tray}:M_nom[i]=Mset[i];
cm1{i in tray}: M_0[i,1] = M_nom[i];

#var M_nom{i in tray}=Mset[i];

cm{i in tray, q in fe: q > 1}: M_0[i,q] = M[i,q-1,ncp];
#var M_0{i in tray, q in fe}=(if q==1 then (M_nom[i])
#                             else (M[i,q-1,ncp]));


#hmm{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : Mdot[i,q,c] = V[i-1,q,c] - V[i,q,c] + L[i+1,q,c]-L[i,q,c];

#hmb{q in fe, c in cp}: Mdot[1,q,c] = L[2,q,c] - L[1,q,c] - V[1,q,c]-Wreb[q,c]+ Wrebl[q,c] + Wrebv[q,c];

#hmfeed{q in fe, c in cp} : Mdot[feedtray,q,c] = V[feedtray-1,q,c] - V[feedtray,q,c] + L[feedtray+1,q,c]-L[feedtray,q,c]+feed[q,c];

#hmc{q in fe, c in cp}: Mdot[Ntray,q,c] = V[Ntray-1,q,c] - L[Ntray,q,c]-( D[q] + eD[q] ) ;

var Mdot{i in tray,q in fe, c in cp}= (if i==1 then (L[2,q,c] - L[1,q,c] - V[1,q,c]-Wreb[q,c]+ Wrebl[q,c] + Wrebv[q,c])
                                        else if i==feedtray then (V[feedtray-1,q,c] - V[feedtray,q,c] + L[feedtray+1,q,c]-L[feedtray,q,c]+feed[q,c])
                                        else if i==Ntray then (V[Ntray-1,q,c] - L[Ntray,q,c]-( D[q] + eD[q] ))
                                        else (V[i-1,q,c] - V[i,q,c] + L[i+1,q,c]-L[i,q,c]));

#hmb{q in fe, c in cp}: 0= L[2,q,c] - L[1,q,c] - V[1,q,c];

#hmb1{q in fe, c in cp} : 80 <= L[1,q,c]<= 380;



#hmc{q in fe, c in cp}: 0= V[Ntray-1,q,c] - L[Ntray,q,c] - D[q,c] ;

fm{i in tray, q in fe, c in cp}: M[i,q,c] = M_0[i,q] + h[q]* sum{k in cp} omega[k,c]*Mdot[i,q,k];

hrc{q in fe, c in cp}: 0 = (D[q]+eD[q])*Rec[q,c] - L[Ntray,q,c];

hrc1{q in fe, c in cp}: 0 = (Reflux[q]+eReflux[q]) - L[Ntray,q,c];

hrc2{q in fe, c in cp}: 0 = (B[q]+eB[q]) - L[1,q,c];

# component balance

#initial constraints
cx1_init_constr{i in tray}:x_nom[i]=xset[i];
cx1{i in tray}: x_0[i,1] = x_nom[i];

#var x_nom{i in tray}=xset[i];

cx{i in tray, q in fe: q > 1} : x_0[i,q] = x[i,q-1,ncp];
#var x_0{i in tray, q in fe}=(if q==1 then (x_nom[i])
#                             else (x[i,q-1,ncp]));

var xdot{i in tray, q in fe, c in cp}=(if i==1 then ((L[2,q,c]*(x[2,q,c]-x[1,q,c])-V[1,q,c]*(y[1,q,c]-x[1,q,c])+ Wrebl[q,c]*(xout[q,c]-x[1,q,c]) + Wrebv[q,c]*(yout[q,c]-x[1,q,c]))/M[1,q,c])
                                       else if i==feedtray then ((V[feedtray-1,q,c]*(y[feedtray-1,q,c]-x[feedtray,q,c])+L[feedtray+1,q,c]*(x[feedtray+1,q,c]-x[feedtray,q,c])-V[feedtray,q,c]*(y[feedtray,q,c]-x[feedtray,q,c]) + feed[q,c]*(xf-x[feedtray,q,c]))/M[feedtray,q,c])
                                       else if i==Ntray then ((V[Ntray-1,q,c]*(y[Ntray-1,q,c]-x[Ntray,q,c]))/M[Ntray,q,c])
                                       else ((V[i-1,q,c]*(y[i-1,q,c]-x[i,q,c])+L[i+1,q,c]*(x[i+1,q,c]-x[i,q,c])-V[i,q,c]*(y[i,q,c]-x[i,q,c]))/M[i,q,c]));

##gx{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : M[i,q,c]*xdot[i,q,c] = V[i-1,q,c]*(y[i-1,q,c]-x[i,q,c])+L[i+1,q,c]*(x[i+1,q,c]-x[i,q,c])-V[i,q,c]*(y[i,q,c]-x[i,q,c]);



##gxfeed{q in fe, c in cp} : M[feedtray,q,c]*xdot[feedtray,q,c] = V[feedtray-1,q,c]*(y[feedtray-1,q,c]-x[feedtray,q,c])+L[feedtray+1,q,c]*(x[feedtray+1,q,c]-x[feedtray,q,c])-V[feedtray,q,c]*(y[feedtray,q,c]-x[feedtray,q,c]) + feed[q,c]*(xf-x[feedtray,q,c]);

##gxb{q in fe, c in cp} : M[1,q,c]*xdot[1,q,c] = L[2,q,c]*(x[2,q,c]-x[1,q,c])-V[1,q,c]*(y[1,q,c]-x[1,q,c])+ Wrebl[q,c]*(xout[q,c]-x[1,q,c]) + Wrebv[q,c]*(yout[q,c]-x[1,q,c]);

##gxc{q in fe, c in cp} : M[Ntray,q,c]*xdot[Ntray,q,c] = V[Ntray-1,q,c]*(y[Ntray-1,q,c]-x[Ntray,q,c]);

fx{i in tray, q in fe, c in cp} : x[i,q,c] = x_0[i,q] + h[q]*sum{k in cp} omega[k,c]*xdot[i,q,k];

# define Tdot
 
##lTdot{i in tray, q in fe, c in cp} : Tdot[i,q,c] = -(pm[i,q,c] - pn[i,q,c])*xdot[i,q,c]/(x[i,q,c]*exp(CapAm-CapBm/(T[i,q,c]))*CapBm/(T[i,q,c])^2+(1-x[i,q,c])*exp(CapAn-CapBn/(T[i,q,c]))*CapBn/(T[i,q,c])^2);
var Tdot{i in tray, q in fe, c in cp} = -(pm[i,q,c] - pn[i,q,c])*xdot[i,q,c]/(x[i,q,c]*exp(CapAm-CapBm/(T[i,q,c]))*CapBm/(T[i,q,c])^2+(1-x[i,q,c])*exp(CapAn-CapBn/(T[i,q,c]))*CapBn/(T[i,q,c])^2);

Clb{q in fe, c in cp} : M[1,q,c] = 300;
Clt{q in fe, c in cp} : M[Ntray,q,c] = 300;

# energy balance

gh{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : M[i,q,c] * (((Alm * (T[i,q,c] - Tref) + Blm / 2 * (T[i,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[i,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[i,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[i,q,c]- Tref) + Bln / 2 * (T[i,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[i,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[i,q,c] ^ 4 - Tref ^ 4))) * xdot[i,q,c] + (x[i,q,c] * (Alm + Blm * T[i,q,c] + Clm * T[i,q,c] ^ 2 + Dlm * T[i,q,c] ^ 3) + (1 - x[i,q,c]) * (Aln + Bln * T[i,q,c] + Cln * T[i,q,c] ^ 2 + Dln * T[i,q,c] ^ 3)) * Tdot[i,q,c]) = V[i-1,q,c]*(hv[i-1,q,c]-hl[i,q,c])+L[i+1,q,c]*(hl[i+1,q,c]-hl[i,q,c])-V[i,q,c]*(hv[i,q,c]-hl[i,q,c]);

ghfeed{q in fe, c in cp} : M[feedtray,q,c] * (((Alm * (T[feedtray,q,c] - Tref) + Blm / 2 * (T[feedtray,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[feedtray,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[feedtray,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[feedtray,q,c]- Tref) + Bln / 2 * (T[feedtray,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[feedtray,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[feedtray,q,c] ^ 4 - Tref ^ 4))) * xdot[feedtray,q,c] + (x[feedtray,q,c] * (Alm + Blm * T[feedtray,q,c] + Clm * T[feedtray,q,c] ^ 2 + Dlm * T[feedtray,q,c] ^ 3) + (1 - x[feedtray,q,c]) * (Aln + Bln * T[feedtray,q,c] + Cln * T[feedtray,q,c] ^ 2 + Dln * T[feedtray,q,c] ^ 3)) * Tdot[feedtray,q,c]) = V[feedtray-1,q,c] * (hv[feedtray-1,q,c]-hl[feedtray,q,c]) + L[feedtray+1,q,c] * (hl[feedtray+1,q,c]-hl[feedtray,q,c]) - V[feedtray,q,c] * (hv[feedtray,q,c]-hl[feedtray,q,c]) + feed[q,c]*(hf-hl[feedtray,q,c]);


ghb{q in fe, c in cp}   :   M[1,q,c] * (((Alm * (T[1,q,c] - Tref) + Blm / 2 * (T[1,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[1,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[1,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[1,q,c]- Tref) + Bln / 2 * (T[1,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[1,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[1,q,c] ^ 4 - Tref ^ 4))) * xdot[1,q,c] + (x[1,q,c] * (Alm + Blm * T[1,q,c] + Clm * T[1,q,c] ^ 2 + Dlm * T[1,q,c] ^ 3) + (1 - x[1,q,c]) * (Aln + Bln * T[1,q,c] + Cln * T[1,q,c] ^ 2 + Dln * T[1,q,c] ^ 3)) * Tdot[1,q,c])= L[2,q,c]*(hl[2,q,c]-hl[1,q,c]) - V[1,q,c]*(hv[1,q,c]-hl[1,q,c])+Wrebl[q,c]*(hlout[q,c]-hl[1,q,c])+ Wrebv[q,c]*(hvout[q,c]-hl[1,q,c]);


ghc{q in fe, c in cp}   :  M[Ntray,q,c] * (((Alm * (T[Ntray,q,c] - Tref) + Blm / 2 * (T[Ntray,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[Ntray,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[Ntray,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[Ntray,q,c]- Tref) + Bln / 2 * (T[Ntray,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[Ntray,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[Ntray,q,c] ^ 4 - Tref ^ 4))) * xdot[Ntray,q,c] + (x[Ntray,q,c] * (Alm + Blm * T[Ntray,q,c] + Clm * T[Ntray,q,c] ^ 2 + Dlm * T[Ntray,q,c] ^ 3) + (1 - x[Ntray,q,c]) * (Aln + Bln * T[Ntray,q,c] + Cln * T[Ntray,q,c] ^ 2 + Dln * T[Ntray,q,c] ^ 3)) * Tdot[Ntray,q,c])= V[Ntray-1,q,c]*(hv[Ntray-1,q,c]-hl[Ntray,q,c])-Qc[q,c]*1e6;



# define the equilibrium coefficients

#gy{i in tray, q in fe ,c in cp} : y[i,q,c] = ((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* Eff + x[i,q,c];

#Minimizatin function 
var sumeps_trays=sum{i in tray}(sum{q in fe} (h[q]*sum{c in cp}(Meps[i,q,c]+xeps[i,q,c]+Teps[i,q,c]+Veps[i,q,c]+Leps[i,q,c]+yeps[i,q,c])*omega[c,ncp]));
var sumeps_others=sum{q in fe}(h[q]*sum{c in cp}(Receps[q,c]+Qreps[q,c]+Qceps[q,c]+Tseps[q,c]+Fseps[q,c]+Pchesteps[q,c]+Wrebeps[q,c]+Wrebleps[q,c]+Wrebveps[q,c]+rorebeps[q,c]+rorebmeps[q,c]+rowrebeps[q,c]+Mvrebeps[q,c]+Mlrebeps[q,c]+Vvrebeps[q,c]+Vlrebeps[q,c])+Deps[q]);
#minimize profit: sum{i in tray}(sum{q in fe} (h[q]*sum{c in cp} (0.001*((D[q] +eD[q])-Dset[q])^2 + (100000*x[Ntray,q,c]- 100000*xntrayref[q])^2+(1000*x[1,q,c]- 1000*x1ref[q])^2+1e5*Meps[i,q,c]+1e5*xeps[i,q,c])*omega[c,ncp]))+sum{f in fe diff{nfe}}(10000*(B[f+1]-B[f])^2+10000*(Psteam[f+1]-Psteam[f])^2);
minimize profit: sum{i in tray}(sum{q in fe} (0*((D[q] +eD[q])-Dset[q])^2 +h[q]*sum{c in cp} ( (1000000*x[Ntray,q,c]- 1000000*xntrayref[q])^2+(1000*x[1,q,c]- 1000*x1ref[q])^2)*omega[c,ncp]))+sum{f in fe diff{nfe}}(10000*(B[f+1]-B[f])^2+10000*(Psteam[f+1]-Psteam[f])^2)+1e6*sumeps_trays+1e6*sumeps_others;

#+1E8*(0.001*M[1,q,c]- 0.001*300)^2+1E8*(0.001*M[Ntray,q,c]- 0.001*300)^2

#Hydrodynamics

var rol{i in tray, q in fe, c in cp} = x[i,q,c]*1000/MWm*(Adm*(Bdm^-(1-T[i,q,c]/Tkm)^Ndm))+(1-x[i,q,c])*1000/MWn*(Adn*(Bdn^-(1-T[i,q,c]/Tkn)^Ndn));

var Vl{i in tray diff{1,Ntray},q in fe,c in cp} = M[i,q,c]/rol[i,q,c];

Hyd{i in tray diff{1,Ntray},q in fe, c in cp}: L[i,q,c]/rol[i,q,c] = (70+eL[q,c])*(Vl[i,q,c]-(0.8+eWeir[q,c]))^1.5;

#Reboiler

cFs{q in fe, c in cp}:  Fs[q,c] = 14.531*(Psteam[q]-Pchest[q,c])^0.5- 213.94;  #flow of vapor
var Hwv{q in fe, c in cp} = ((2.7376-2.7247) * (Psteam[q]-300))/(400-300)+ 2.7247; #eq in GJ/ton
var Hwl{q in fe, c in cp} = 3.7982E-3*(Ts[q,c]-273.15); 

cchest{q in fe, c in cp}: Ts[q,c]=((341.85-318.55)*(Pchest[q,c]-10))/(30-10)+318.55;

#reboiler volume 4.835784 m^3
#var Areb {q in fe, c in cp} = Atotal*(1-(Vreb[q,c]/4.835784));
param Areb {q in fe, c in cp} = Atotal*(1-(3/4.835784));

cQr0{q in fe, c in cp}: kQr*Fs[q,c]*(Hwv[q,c]- Hwl[q,c])= (Uheat+eU[q,c])*Areb[q,c]*(Ts[q,c]- ((Tout[q,c]+ T[1,q,c])/2));

cQr1{q in fe, c in cp}: Qr[q,c]= (Uheat+eU[q,c])*Areb[q,c]*(Ts[q,c]- ((Tout[q,c]+ T[1,q,c])/2));  #flux of heat

rYreboiler1{q in fe, c in cp}: yout[q,c] = ((xout[q,c]*1.1)/(1+(1.1-1)*xout[q,c]));

var rolm{i in tray, q in fe, c in cp} = x[i,q,c]*1000*(Adm*(Bdm^-(1-T[i,q,c]/Tkm)^Ndm))+(1-x[i,q,c])*1000*(Adn*(Bdn^-(1-T[i,q,c]/Tkn)^Ndn));
#rolm=rol*Mwn?

var rorebl{q in fe, c in cp} = pout*100/((xout[q,c]*0.056+(1-xout[q,c])*0.059)*1*r*Tout[q,c]);

var rorebv{q in fe, c in cp} = pout*100/((yout[q,c]*0.6+(1-yout[q,c])*0.63)*1*r*Tout[q,c]);

rro0{q in fe, c in cp} : rorebm[q,c]=roreb[q,c]*(x[1,q,c]*MWm+(1- x[1,q,c])*MWn);

rro{q in fe, c in cp} : Wreb[q,c]/roreb[q,c]= Wrebv[q,c]/rorebv[q,c]+ Wrebl[q,c]/rorebl[q,c];

#pipe area 0.2918 m^2 tube 24�
var Vel{q in fe, c in cp} =(Wreb[q,c]/3600)/(roreb[q,c]*0.2918);

rmbg5{q in fe, c in cp}: ((p[1]*100000+rolm[1,q,c]*9.81*3))/ rorebm[q,c]-(p[1]*100000+rorebm[q,c]*9.81*3)/rorebm[q,c]=0.5*Vel[q,c]^2;

rebreboiler0{q in fe, c in cp}: 0 = Wreb[q,c]- Wrebl[q,c]- Wrebv[q,c];

rebreboiler1{q in fe, c in cp}: 0 = Wreb[q,c]*x[1,q,c]-Wrebl[q,c]*xout[q,c]-Wrebv[q,c]*yout[q,c];

rebreboiler{q in fe, c in cp}: 0=1E6*((Uheat+eU[q,c])*Areb[q,c]*(Ts[q,c]-((Tout[q,c]+ T[1,q,c])/2)))+Wreb[q,c]*hl[1,q,c]- Wrebl[q,c]*hlout[q,c]- Wrebv[q,c]*hvout[q,c];

