
#only add soft constraints to M,x,M_0,x_0, M_nom and x_nom

param nfe;  
param ncp;
param time; # horizon length

set fe := 1..nfe ;  
set cp := 1..ncp ; 

param omega{cp,cp};   param h{fe}; 
param alpha;
param beta;

param Ntray          ;
param feedtray;
set tray := 1..Ntray  ;
#let Ntray := 158;
#let feedtray := 43;
#let nfe       := 18; 
#let ncp       := 3;
let ncp :=2;

param kQr>=0,<=3;     #cte for Efficiency of Reboiler
param Eff;
param kwreb >=1E-6,<=1E6;
param kFlow >=1E-6,<=1E6;
param Uheat >=0,<=10;		#Coefficient heat transfer


# define the hold in each tray
var   M{tray,fe,cp};
#var Meps{tray,fe,cp}>=0;
Mupper{i in tray, j in fe, k in cp}: M[i,j,k]<=500;
Mlower{i in tray, j in fe, k in cp}: M[i,j,k]>=20;

var   Mdot{tray, fe,cp};
var M_0{tray, fe};
#var M0eps{tray,fe}>=0;
M0upper{i in tray, j in fe}: M_0[i,j]<=500;
M0lower{i in tray, j in fe}: M_0[i,j]>=20;
##var M_0{tray, fe}>=20,<=500;


# define the temperature in each tray
##var T{tray, fe, cp};
##var Teps{tray,fe,cp}>=0;
##Tupper{i in tray, j in fe, k in cp}: T[i,j,k]<=335+Teps[i,j,k];
##Tlower{i in tray, j in fe, k in cp}: T[i,j,k]>=280-Teps[i,j,k];
##var Tbottom{fe, cp};
##var Tboeps{fe,cp}>=0;
##Tboupper{j in fe, k in cp}: Tbottom[j,k]<=335+Tboeps[j,k];
##Tbolower{j in fe, k in cp}: Tbottom[j,k]>=280-Tboeps[j,k];
##var Ttop{fe, cp};
##var Ttoeps{fe,cp}>=0;
##Ttoupper{j in fe, k in cp}: Ttop[j,k]<=335+Ttoeps[j,k];
##Ttolower{j in fe, k in cp}: Ttop[j,k]>=280-Ttoeps[j,k];

var T{tray, fe, cp} >=280, <=335;
#var Tbottom{fe, cp} >=280, <=335;
#var Ttop{fe, cp} >=280, <=335;

var Tdot{tray, fe,cp} ;

param Tref;

param Ttop_des; param Tbottom_des;

param feednmpc{fe};

#parameters to MHE
#these two lines are just if I intend to consider error on states
#var extop{fe,cp}>=-0.05,<=0.05; var exbottom{fe,cp}>=-0.1,<=0.1; 
#var eMtop{fe,cp}>=-5,<=5; var eMbottom{fe,cp}>=-5,<=5; 
#var eD{fe,cp}>=-10,<=10; var eReflux{fe,cp}>=-40,<=40; 
#var eB{fe,cp}>=-10,<=10; 
param eD{fe}:=0; param eReflux{fe}:=0; param eB{fe}:=0;
#var eTtop{fe,cp}>=-2,<=2; var eTbottom{fe,cp}>=-2,<=2; 
#var eEff{fe,cp}>=-0.2,<=0.2; var eU{fe,cp}>=-0.1,<=0.1; 
#var eL{fe,cp}>=-10,<=10; var eWeir{fe,cp}>=-0.1,<=0.1; 
param eEff{fe,cp}:=0; param eU{fe,cp}:=0; 
param eL{fe,cp}:=0; param eWeir{fe,cp}:=0; 
#var eTout{fe,cp}>=-3,<=3; 

#parameters to simulation
param xtop{fe,cp}; param xbottom{fe,cp}; param Mtop{fe,cp}; param Mbottom{fe,cp}; param Mset{tray}; param xset{tray}; param Refluxset{fe}; param feedset{fe}; param Psteamset{fe}; param Dset{fe}; param Qcset{fe,cp}; param Bset{fe};
param Ttopset{fe,cp}; param Tbottomset{fe,cp}; 

# define the pressure
param p{tray}>=15, <= 20;
#this max restriction to P was added to avoid negative sqrt
#param p{tray}>=0, <= 3000000;
param prect; param pstrip; 

param CapAm; param CapBm; 
param CapAn; param CapBn; 


# define the liquid and vapour flowrate in each tray
##var V{tray,fe,cp};
##var Veps{tray,fe,cp}>=0;
##Vupper{i in tray, j in fe, k in cp}: V[i,j,k]<=12000+Veps[i,j,k];
##Vlower{i in tray, j in fe, k in cp}: V[i,j,k]>=-Veps[i,j,k];
param Vup{tray}; param Vlow{tray};    #bounds on V
param vlow;
var V{i in tray,fe,cp}>=Vlow[i],<=Vup[i];
#var V{i in tray,fe,cp}>=0,<=12000;
param Lup{tray};#:=12000;
param Llow{tray};#:=0;
var L{i in tray,fe,cp}>=Llow[i],<=Lup[i];
##var L{i in tray,fe,cp};
##var Leps{tray,fe,cp}>=0;
##Lupper{i in tray, j in fe, k in cp}: L[i,j,k]<=Lup[i]+Leps[i,j,k];
##Llower{i in tray, j in fe, k in cp}: L[i,j,k]>=Llow[i]-Leps[i,j,k];
#hmb1{q in fe, c in cp} : 80 <= L[1,q,c]<= 380;

##var Rec{fe,cp};
##var Receps{fe,cp}>=0;
##Recupper{i in fe, j in cp}: Rec[i,j] <=18+Receps[i,j];
##Reclower{i in fe, j in cp}: Rec[i,j] >=8-Receps[i,j];
var Rec{fe,cp}<=18;
#var Rec{fe,cp}>=8,<=80;
var Reflux{fe};


# define the molar fraction in liquid and vapour
param xup{tray}; param xlow{tray};      
var x{tray,fe,cp};
#var xeps{tray,fe,cp}>=0;
xupper{i in tray diff{1,Ntray}, j in fe, k in cp}: x[i,j,k]<=xup[i];
xlower{i in tray diff{1,Ntray}, j in fe, k in cp}: x[i,j,k]>=xlow[i];
var x1eps{fe,cp}>=0;
var xNtrayeps{fe,cp}>=0;
xlower1{j in fe, c in cp}:x[Ntray,j,c]>=0.99-xNtrayeps[j,c];
xupper1{j in fe, c in cp}:x[1,j,c]<=0.05+x1eps[j,c];
xlower2{j in fe, c in cp}:x[1,j,c]>=0;
xupper2{j in fe, c in cp}:x[Ntray,j,c]<=1;


##var x_0{tray, fe} >=0 , <=1;
##var x_nom{tray} >=0 , <=1;
var y{tray,fe,cp}>=0, <=1;

var xdot{tray, fe, cp};
var x_0{tray, fe};
#var x0eps{tray,fe}>=0;
x0upper{i in tray, j in fe}: x_0[i,j]<=xup[i];
x0lower{i in tray, j in fe}: x_0[i,j]>=xlow[i];

##var y{tray,fe,cp};
##var yeps{tray,fe,cp}>=0;
##yupper{i in tray, j in fe, k in cp}: y[i,j,k]<=1+yeps[i,j,k];
##ylower{i in tray, j in fe, k in cp}: y[i,j,k]>=-yeps[i,j,k];

# define the enthalpy 
#var hl{tray, fe,cp}>=5000;
#var hv{tray,fe,cp}>=15000 ;

# define the 
##var Qr{fe,cp};
##var Qreps{fe,cp}>=0;
##Qrupper{j in fe, k in cp}: Qr[j,k]<=140+Qreps[j,k];
##Qrlower{j in fe, k in cp}: Qr[j,k]>=40-Qreps[j,k];

##var Qc{fe,cp};
##var Qceps{fe,cp}>=0;
##Qclower{j in fe, k in cp}: Qc[j,k]>=40-Qceps[j,k];

##var D{fe};
##var Deps{fe}>=0;
##Dupper{j in fe}: D[j]<=500+Deps[j];
##Dlower{j in fe}: D[j]>=200-Deps[j];

##var B{fe};
##var Beps{fe}>=0;
##Bupper{j in fe}: B[j]<=500+Beps[j];
##Blower{j in fe}: B[j]>=50-Beps[j];

var Qr{fe,cp}>=40,<=140;
var Qc{fe,cp}>=40;
#var D{fe}>=200,<=600;
#var B{fe}>=50,<=500;
var D{fe}>=300;#,<=500;
var B{fe}>=50,<=380;

# define the feed
#var feed{fe,cp};
#cfeed0{q in fe,c in cp}: feed[q]=feednmpc[q];
param feed{q in fe}=feednmpc[q];

param xf; # feed molar fraction
param hf; # feed enthalpy

#Reboiler 
##var Psteam{fe};	#Vapor pressure of Steam
##var Peps{fe}>=0;
##Pupper{j in fe}: Psteam[j]<=400+Peps[j];
##Plower{j in fe}: Psteam[j]>=300-Peps[j];

##var Ts{fe,cp} ; 	#Temperature steam
##var Tseps{fe,cp}>=0;
##Tsupper{j in fe, k in cp}: Ts[j,k]<=360+Tseps[j,k];
##Tslower{j in fe, k in cp}: Ts[j,k]>=310-Tseps[j,k];

##var Fs{fe,cp}>=0,<=100;		#Flow vapor reboiler
##var Fseps{fe,cp}>=0;
##Fsupper{j in fe, k in cp}: Fs[j,k]<=100+Fseps[j,k];
##Fslower{j in fe, k in cp}: Fs[j,k]>=-Fseps[j,k];

var Psteam{fe}>=100,<=400;	#Vapor pressure of Steam
var Ts{fe,cp} >=310 , <=360; 	#Temperature steam
var Fs{fe,cp}>=0,<=100;		#Flow vapor reboiler

param pout>=15,<=60;			#Pressure reboiler
##var Pchest{fe,cp};		#Pressure leaving reboiler
##var Pceps{fe,cp}>=0;
##Pcupper{j in fe, k in cp}: Pchest[j,k]<=30+Pceps[j,k];
##Pclower{j in fe, k in cp}: Pchest[j,k]>=10-Pceps[j,k];

##var Tout{fe,cp}; 	#Out Temperature of reboiler
##var Toeps{fe,cp}>=0;
##Toupper{j in fe, k in cp}: Tout[j,k]<=430+Toeps[j,k];
##Tolower{j in fe, k in cp}: Tout[j,k]>=290-Toeps[j,k];

##var xout{fe,cp};		#Pressure reboiler
##var xoeps{fe,cp}>=0;
##xoupper{j in fe, k in cp}: xout[j,k]<=1+xoeps[j,k];
##xolower{j in fe, k in cp}: xout[j,k]>=-xoeps[j,k];

##var yout{fe,cp}>=0,<=1; 		#Out Temperature of reboiler
##var yoeps{fe,cp}>=0;
##youpper{j in fe, k in cp}: yout[j,k]<=1+yoeps[j,k];
##yolower{j in fe, k in cp}: yout[j,k]>=-yoeps[j,k];

##var Wreb{fe,cp};		#Flow reboiler
##var Weps{fe,cp}>=0;
##Wupper{j in fe, k in cp}: Wreb[j,k]<=500000+Weps[j,k];
##Wlower{j in fe, k in cp}: Wreb[j,k]>=-Weps[j,k];

##var Wrebl{fe,cp};	#Flow reboiler liquid
##var Wleps{fe,cp}>=0;
##Wlupper{j in fe, k in cp}: Wrebl[j,k]<=500000+Wleps[j,k];
##Wllower{j in fe, k in cp}: Wrebl[j,k]>=-Wleps[j,k];

##var Wrebv{fe,cp};	#Flow reboiler vapor
##var Wveps{fe,cp}>=0;
##Wvupper{j in fe, k in cp}: Wrebv[j,k]<=500000+Wveps[j,k];
##Wvlower{j in fe, k in cp}: Wrebv[j,k]>=-Wveps[j,k];

##var roreb{fe,cp};		#ro leaving reboiler
##var reps{fe,cp}>=0;
##rupper{j in fe, k in cp}: roreb[j,k]<=15+reps[j,k];
##rlower{j in fe, k in cp}: roreb[j,k]>=0.1-reps[j,k];

##var rorebm{fe,cp};		#ro leaving reboiler
##var rmeps{fe,cp}>=0;
##rmupper{j in fe, k in cp}: rorebm[j,k]<=550+rmeps[j,k];
##rmlower{j in fe, k in cp}: rorebm[j,k]>=-rmeps[j,k];

##var rowreb{fe,cp};	#ro reboiler total
##var rbeps{fe,cp}>=0;
##rbupper{j in fe, k in cp}: rowreb[j,k]<=15+rbeps[j,k];
##rblower{j in fe, k in cp}: rowreb[j,k]>=0.1-rbeps[j,k];

##var Mvreb{fe,cp};		#Hold up vapor
##var Mveps{fe,cp}>=0;
##Mvupper{j in fe, k in cp}: Mvreb[j,k]<=100+Mveps[j,k];
##Mvlower{j in fe, k in cp}: Mvreb[j,k]>=-Mveps[j,k];

##var Mlreb{fe,cp};		#Hold up liquid
##var Mleps{fe,cp}>=0;
##Mlupper{j in fe, k in cp}: Mlreb[j,k]<=100+Mleps[j,k];
##Mllower{j in fe, k in cp}: Mlreb[j,k]>=-Mleps[j,k];

##var Vvreb{fe,cp};		#Volume vapor on reboiler
##var Vveps{fe,cp}>=0;
##Vvupper{j in fe, k in cp}: Vvreb[j,k]<=100+Vveps[j,k];
##Vvlower{j in fe, k in cp}: Vvreb[j,k]>=-Vveps[j,k];

##var Vlreb{fe,cp};		#Volume liquid on reboiler
##var Vleps{fe,cp}>=0;
##Vlupper{j in fe, k in cp}: Vlreb[j,k]<=100+Vleps[j,k];
##Vllower{j in fe, k in cp}: Vlreb[j,k]>=-Vleps[j,k];

var Pchest{fe,cp}>=10,<=30;		#Pressure leaving reboiler
var Tout{fe,cp} >=290 , <=430; 	#Out Temperature of reboiler
var xout{fe,cp}>=0,<=1;		#Pressure reboiler
var yout{fe,cp}>=0,<=1; 		#Out Temperature of reboiler
var Wreb{fe,cp}>=0,<=500000;		#Flow reboiler
var Wrebl{fe,cp}>=0,<=500000;	#Flow reboiler liquid
var Wrebv{fe,cp}>=0,<=500000;	#Flow reboiler vapor

var roreb{fe,cp}>=0.1, <=15;		#ro leaving reboiler
var rorebm{fe,cp}>=0, <=550;		#ro leaving reboiler
var rowreb{fe,cp}>=0.1, <=15;	#ro reboiler total


#var Mvreb{fe,cp}>=0, <=100;		#Hold up vapor
#var Mlreb{fe,cp}>=0, <=100;		#Hold up liquid
#var Vvreb{fe,cp}>=0, <=100;		#Volume vapor on reboiler
#var Vlreb{fe,cp}>=0, <=100;		#Volume liquid on reboiler

param Atotal;				#Total area reboiler

param x1ref{fe};
param xntrayref{fe};
param M1ref{fe};
param Mntrayref{fe};
param xref{tray,fe};
param yref{tray,fe};
param Mref{tray,fe};
param Dref{fe};
param Bref{fe};
param Psteamref{fe};
param Refluxref{fe};
param Trefe{tray,fe};
param Lref{tray,fe};
param Vref{tray diff{Ntray},fe};
param Recref{fe};
param Qcref{fe};
param Tsref{fe};
param Pchestref{fe};
param Toutref{fe};
param Fsref{fe};
param xoutref{fe};
param youtref{fe};
param Wrebref{fe};
param rorebref{fe};
param Wrebvref{fe};
param Wreblref{fe};
param rorebmref{fe};
param Qrref{fe};

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

#terminal constraints
#conx1l{k in cp}: x[1,nfe,k]=x1ref[nfe];
#conx1u{k in cp}: x[1,nfe,k]<=0.21;
#conxNtrayl{k in cp}: x[Ntray,nfe,k]=xntrayref[nfe];
#conxNtrayu{k in cp}: x[Ntray,nfe,k]<=0.996;

# pressure
#partial pressure p, see Diehl's P116

var pm{i in tray,q in fe,c in cp}=	 exp(CapAm-CapBm/(T[i,q,c]));

var pn{i in tray,q in fe,c in cp}= exp(CapAn-CapBn/(T[i,q,c]));   

dp{i in tray, q in fe, c in cp} : p[i] = pm[i,q,c]*x[i,q,c]+(1-x[i,q,c])*pn[i,q,c];


# define the equilibrium coefficients

yr{i in tray, q in fe ,c in cp}: y[i,q,c] = ((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* (Eff+eEff[q,c]) + x[i,q,c] ;

#var y{i in tray,q in fe, c in cp}=((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* (Eff+eEff[q,c]) + x[i,q,c] ;
#yr{i in tray, q in fe, c in cp} : 0 <= y[i,q,c] <=1 ;


#restriction monotomic Temperature

#mT{i in tray diff{1,Ntray}, q in fe, c in cp} : T[i+1,q,c]<= T[i,q,c];

# calculate the hkl and hkv

var hl{i in tray, q in fe ,c in cp}= x[i,q,c]*(Alm*(T[i,q,c]-Tref)  + Blm/2*(T[i,q,c]^2-Tref^2) + Clm/3*(T[i,q,c]^3-Tref^3) + Dlm/4*(T[i,q,c]^4-Tref^4) ) + (1-x[i,q,c])*(Aln*(T[i,q,c]-Tref) + Bln/2*(T[i,q,c]^2-Tref^2) + Cln/3*(T[i,q,c]^3-Tref^3) + Dln/4*(T[i,q,c]^4-Tref^4) );

var hlout{q in fe ,c in cp}= xout[q,c]*(Alm*(Tout[q,c]-Tref)  + Blm/2*(Tout[q,c]^2-Tref^2) + Clm/3*(Tout[q,c]^3-Tref^3) + Dlm/4*(Tout[q,c]^4-Tref^4) ) + (1-xout[q,c])*(Aln*(Tout[q,c]-Tref) + Bln/2*(Tout[q,c]^2-Tref^2) + Cln/3*(Tout[q,c]^3-Tref^3) + Dln/4*(Tout[q,c]^4-Tref^4) );


var hv{i in tray, q in fe ,c in cp}= y[i,q,c]*((h0m+Avm*(T[i,q,c]-Tref)  + Bvm/2*(T[i,q,c]^2-Tref^2) + Cvm/3*(T[i,q,c]^3-Tref^3) + Dvm/4*(T[i,q,c]^4-Tref^4) + Evm/5*(T[i,q,c]^5-Tref^5))+(4.18*(r*Tkm/Pkm*(0.083-0.422/(T[i,q,c]/Tkm)^(8/5)+gm*(0.139-0.172/(T[i,q,c]/Tkm)^(21/5)))-T[i,q,c]*(0.672/(T[i,q,c]/Tkm)^(13/5)/Tkm+0.714*gm/(T[i,q,c]/Tkm)^(26/5)/Tkm))*p[i]))+(1-y[i,q,c])*((h0n+Avn*(T[i,q,c]-Tref)  + Bvn/2*(T[i,q,c]^2-Tref^2) + Cvn/3*(T[i,q,c]^3-Tref^3) + Dvn/4*(T[i,q,c]^4-Tref^4) + Evn/5*(T[i,q,c]^5-Tref^5))+(4.18*(r*Tkn/Pkn*(0.083-0.422/(T[i,q,c]/Tkn)^(8/5)+gn*(0.139-0.172/(T[i,q,c]/Tkn)^(21/5)))-T[i,q,c]*(0.672/(T[i,q,c]/Tkn)^(13/5)/Tkn+0.714*gn/(T[i,q,c]/Tkn)^(26/5)/Tkn))*p[i]));

var hvout{q in fe ,c in cp}= yout[q,c]*((h0m+Avm*(Tout[q,c]-Tref)  + Bvm/2*(Tout[q,c]^2-Tref^2) + Cvm/3*(Tout[q,c]^3-Tref^3) + Dvm/4*(Tout[q,c]^4-Tref^4) + Evm/5*(Tout[q,c]^5-Tref^5))+(4.18*(r*Tkm/Pkm*(0.083-0.422/(Tout[q,c]/Tkm)^(8/5)+gm*(0.139-0.172/(Tout[q,c]/Tkm)^(21/5)))-Tout[q,c]*(0.672/(Tout[q,c]/Tkm)^(13/5)/Tkm+0.714*gm/(Tout[q,c]/Tkm)^(26/5)/Tkm))*p[1]))+(1-yout[q,c])*((h0n+Avn*(Tout[q,c]-Tref)  + Bvn/2*(Tout[q,c]^2-Tref^2) + Cvn/3*(Tout[q,c]^3-Tref^3) + Dvn/4*(Tout[q,c]^4-Tref^4) + Evn/5*(Tout[q,c]^5-Tref^5))+(4.18*(r*Tkn/Pkn*(0.083-0.422/(Tout[q,c]/Tkn)^(8/5)+gn*(0.139-0.172/(Tout[q,c]/Tkn)^(21/5)))-Tout[q,c]*(0.672/(Tout[q,c]/Tkn)^(13/5)/Tkn+0.714*gn/(Tout[q,c]/Tkn)^(26/5)/Tkn))*pout));


# mass balance
#var M_0{i in tray, q in fe}=(if q==1 then (Mset[i])
#                             else (M[i,q-1,ncp]));
cm1{i in tray}: M_0[i,1] = Mset[i];

cm{i in tray, q in fe: q > 1}: M_0[i,q] = M[i,q-1,ncp];

hmm{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : Mdot[i,q,c] = V[i-1,q,c] - V[i,q,c] + L[i+1,q,c]-L[i,q,c];

hmb{q in fe, c in cp}: Mdot[1,q,c] = L[2,q,c] - L[1,q,c] - V[1,q,c]-Wreb[q,c]+ Wrebl[q,c] + Wrebv[q,c];

hmfeed{q in fe, c in cp} : Mdot[feedtray,q,c] = V[feedtray-1,q,c] - V[feedtray,q,c] + L[feedtray+1,q,c]-L[feedtray,q,c]+feed[q];

#var Mdot{i in tray,q in fe, c in cp}= (if i==1 then (L[2,q,c] - L[1,q,c] - V[1,q,c]-Wreb[q,c]+ Wrebl[q,c] + Wrebv[q,c])
#                                        else if i==feedtray then (V[feedtray-1,q,c] - V[feedtray,q,c] + L[feedtray+1,q,c]-L[feedtray,q,c]+feed[q])
#                                        else if i==Ntray then (V[Ntray-1,q,c] - L[Ntray,q,c]-( D[q] + eD[q] ))
#                                        else (V[i-1,q,c] - V[i,q,c] + L[i+1,q,c]-L[i,q,c]));


hmc{q in fe, c in cp}: Mdot[Ntray,q,c] = V[Ntray-1,q,c] - L[Ntray,q,c]-( D[q] + eD[q] ) ;
#hmc{q in fe, c in cp}: 0= V[Ntray-1,q,c] - L[Ntray,q,c] - D[q,c] ;

fm{i in tray, q in fe, c in cp}: M[i,q,c] = M_0[i,q] + h[q]* sum{k in cp} omega[k,c]*Mdot[i,q,k];

hrc{q in fe, c in cp}: 0 = (D[q]+eD[q])*Rec[q,c] - L[Ntray,q,c];

hrc1{q in fe, c in cp}: 0 = (Reflux[q]+eReflux[q]) - L[Ntray,q,c];




hrc2{q in fe, c in cp}: 0 = (B[q]+eB[q]) - L[1,q,c];

# component balance


#inicial constraints

#var x_0{i in tray, q in fe}=(if q==1 then (xset[i])
#                             else (x[i,q-1,ncp]));
cx1{i in tray}: x_0[i,1] = xset[i];

cx{i in tray, q in fe: q > 1} : x_0[i,q] = x[i,q-1,ncp];

#var xdot{i in tray, q in fe, c in cp}=(if i==1 then ((L[2,q,c]*(x[2,q,c]-x[1,q,c])-V[1,q,c]*(y[1,q,c]-x[1,q,c])+ Wrebl[q,c]*(xout[q,c]-x[1,q,c]) + Wrebv[q,c]*(yout[q,c]-x[1,q,c]))/M[1,q,c])
#                                       else if i==feedtray then ((V[feedtray-1,q,c]*(y[feedtray-1,q,c]-x[feedtray,q,c])+L[feedtray+1,q,c]*(x[feedtray+1,q,c]-x[feedtray,q,c])-V[feedtray,q,c]*(y[feedtray,q,c]-x[feedtray,q,c]) + feed[q]*(xf-x[feedtray,q,c]))/M[feedtray,q,c])
#                                       else if i==Ntray then ((V[Ntray-1,q,c]*(y[Ntray-1,q,c]-x[Ntray,q,c]))/M[Ntray,q,c])
#                                       else ((V[i-1,q,c]*(y[i-1,q,c]-x[i,q,c])+L[i+1,q,c]*(x[i+1,q,c]-x[i,q,c])-V[i,q,c]*(y[i,q,c]-x[i,q,c]))/M[i,q,c]));

gx{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : M[i,q,c]*xdot[i,q,c] = V[i-1,q,c]*(y[i-1,q,c]-x[i,q,c])+L[i+1,q,c]*(x[i+1,q,c]-x[i,q,c])-V[i,q,c]*(y[i,q,c]-x[i,q,c]);



gxfeed{q in fe, c in cp} : M[feedtray,q,c]*xdot[feedtray,q,c] = V[feedtray-1,q,c]*(y[feedtray-1,q,c]-x[feedtray,q,c])+L[feedtray+1,q,c]*(x[feedtray+1,q,c]-x[feedtray,q,c])-V[feedtray,q,c]*(y[feedtray,q,c]-x[feedtray,q,c]) + feed[q]*(xf-x[feedtray,q,c]);

gxb{q in fe, c in cp} : M[1,q,c]*xdot[1,q,c] = L[2,q,c]*(x[2,q,c]-x[1,q,c])-V[1,q,c]*(y[1,q,c]-x[1,q,c])+ Wrebl[q,c]*(xout[q,c]-x[1,q,c]) + Wrebv[q,c]*(yout[q,c]-x[1,q,c]);

gxc{q in fe, c in cp} : M[Ntray,q,c]*xdot[Ntray,q,c] = V[Ntray-1,q,c]*(y[Ntray-1,q,c]-x[Ntray,q,c]);

fx{i in tray, q in fe, c in cp} : x[i,q,c] = x_0[i,q] + h[q]*sum{k in cp} omega[k,c]*xdot[i,q,k];

#Constant level
Clb{q in fe, c in cp} : M[1,q,c] = 300;
Clt{q in fe, c in cp} : M[Ntray,q,c] = 300;
#Clb{q in fe, c in cp} : 290 <= M[1,q,c]<= 310;
#Clt{q in fe, c in cp} : 290 <= M[Ntray,q,c]<= 310;

#define Tdot
lTdot{i in tray, q in fe, c in cp} : Tdot[i,q,c] = -(pm[i,q,c] - pn[i,q,c])*xdot[i,q,c]/(x[i,q,c]*exp(CapAm-CapBm/(T[i,q,c]))*CapBm/(T[i,q,c])^2+(1-x[i,q,c])*exp(CapAn-CapBn/(T[i,q,c]))*CapBn/(T[i,q,c])^2);

#Clb{q in fe, c in cp} : M[1,q,c] = 300;
#Clt{q in fe, c in cp} : M[Ntray,q,c] = 300;
# energy balance

gh{i in tray diff{1,feedtray,Ntray}, q in fe, c in cp} : M[i,q,c] * (((Alm * (T[i,q,c] - Tref) + Blm / 2 * (T[i,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[i,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[i,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[i,q,c]- Tref) + Bln / 2 * (T[i,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[i,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[i,q,c] ^ 4 - Tref ^ 4))) * xdot[i,q,c] + (x[i,q,c] * (Alm + Blm * T[i,q,c] + Clm * T[i,q,c] ^ 2 + Dlm * T[i,q,c] ^ 3) + (1 - x[i,q,c]) * (Aln + Bln * T[i,q,c] + Cln * T[i,q,c] ^ 2 + Dln * T[i,q,c] ^ 3)) * Tdot[i,q,c]) = V[i-1,q,c]*(hv[i-1,q,c]-hl[i,q,c])+L[i+1,q,c]*(hl[i+1,q,c]-hl[i,q,c])-V[i,q,c]*(hv[i,q,c]-hl[i,q,c]);

ghfeed{q in fe, c in cp} : M[feedtray,q,c] * (((Alm * (T[feedtray,q,c] - Tref) + Blm / 2 * (T[feedtray,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[feedtray,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[feedtray,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[feedtray,q,c]- Tref) + Bln / 2 * (T[feedtray,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[feedtray,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[feedtray,q,c] ^ 4 - Tref ^ 4))) * xdot[feedtray,q,c] + (x[feedtray,q,c] * (Alm + Blm * T[feedtray,q,c] + Clm * T[feedtray,q,c] ^ 2 + Dlm * T[feedtray,q,c] ^ 3) + (1 - x[feedtray,q,c]) * (Aln + Bln * T[feedtray,q,c] + Cln * T[feedtray,q,c] ^ 2 + Dln * T[feedtray,q,c] ^ 3)) * Tdot[feedtray,q,c]) = V[feedtray-1,q,c] * (hv[feedtray-1,q,c]-hl[feedtray,q,c]) + L[feedtray+1,q,c] * (hl[feedtray+1,q,c]-hl[feedtray,q,c]) - V[feedtray,q,c] * (hv[feedtray,q,c]-hl[feedtray,q,c]) + feed[q]*(hf-hl[feedtray,q,c]);


ghb{q in fe, c in cp}   :   M[1,q,c] * (((Alm * (T[1,q,c] - Tref) + Blm / 2 * (T[1,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[1,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[1,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[1,q,c]- Tref) + Bln / 2 * (T[1,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[1,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[1,q,c] ^ 4 - Tref ^ 4))) * xdot[1,q,c] + (x[1,q,c] * (Alm + Blm * T[1,q,c] + Clm * T[1,q,c] ^ 2 + Dlm * T[1,q,c] ^ 3) + (1 - x[1,q,c]) * (Aln + Bln * T[1,q,c] + Cln * T[1,q,c] ^ 2 + Dln * T[1,q,c] ^ 3)) * Tdot[1,q,c])= L[2,q,c]*(hl[2,q,c]-hl[1,q,c]) - V[1,q,c]*(hv[1,q,c]-hl[1,q,c])+Wrebl[q,c]*(hlout[q,c]-hl[1,q,c])+ Wrebv[q,c]*(hvout[q,c]-hl[1,q,c]);


ghc{q in fe, c in cp}   :  M[Ntray,q,c] * (((Alm * (T[Ntray,q,c] - Tref) + Blm / 2 * (T[Ntray,q,c] ^ 2 - Tref ^ 2) + Clm / 3 * (T[Ntray,q,c] ^ 3 - Tref ^ 3) + Dlm / 4 * (T[Ntray,q,c] ^ 4 - Tref ^ 4)) - (Aln * (T[Ntray,q,c]- Tref) + Bln / 2 * (T[Ntray,q,c] ^ 2 - Tref ^ 2) + Cln / 3 * (T[Ntray,q,c] ^ 3 - Tref ^ 3) + Dln / 4 * (T[Ntray,q,c] ^ 4 - Tref ^ 4))) * xdot[Ntray,q,c] + (x[Ntray,q,c] * (Alm + Blm * T[Ntray,q,c] + Clm * T[Ntray,q,c] ^ 2 + Dlm * T[Ntray,q,c] ^ 3) + (1 - x[Ntray,q,c]) * (Aln + Bln * T[Ntray,q,c] + Cln * T[Ntray,q,c] ^ 2 + Dln * T[Ntray,q,c] ^ 3)) * Tdot[Ntray,q,c])= V[Ntray-1,q,c]*(hv[Ntray-1,q,c]-hl[Ntray,q,c])-Qc[q,c]*1e6;

# define Tdot
 
#lTdot{i in tray, q in fe, c in cp} : Tdot[i,q,c] = -(pm[i,q,c] - pn[i,q,c])*xdot[i,q,c]/(x[i,q,c]*exp(CapAm-CapBm/(T[i,q,c]))*CapBm/(T[i,q,c])^2+(1-x[i,q,c])*exp(CapAn-CapBn/(T[i,q,c]))*CapBn/(T[i,q,c])^2);

# define the equilibrium coefficients

#gy{i in tray, q in fe ,c in cp} : y[i,q,c] = ((x[i,q,c]*pm[i,q,c]/p[i])-x[i,q,c])* Eff + x[i,q,c];

####end point constraint: the last Nsamp elements should all be identical as set point

#Minimizatin function 

param dd{1..965};
#minimize profit: sum{q in fe} (h[q]*sum{c in cp} (alpha*Qr[q,c]+beta*sum{i in tray diff{1,Ntray}}(dd[i-1]*(M[i,q,c]-Mref[i,q])^2)+beta*sum{j in tray}(dd[j+156]*(T[j,q,c]-Trefe[j,q])^2) + beta*sum{k in tray diff{Ntray}}(dd[k+314]*(V[k,q,c]-Vref[k,q])^2)+beta*sum{l in tray diff{1}}(dd[l+470]*(L[l,q,c]-Lref[l,q])^2) + beta*sum{m in tray}(dd[m+629]*(x[m,q,c]-xref[m,q])^2)+ beta*sum{n in tray}(dd[n+787]*(y[n,q,c]-yref[n,q])^2)+beta*(dd[629]*(Rec[q,c]-Recref[q])^2+dd[947]*(Ts[q,c]-Tsref[q])^2+dd[948]*(Fs[q,c]-Fsref[q])^2+dd[949]*(Pchest[q,c]-Pchestref[q])^2+dd[950]*(Tout[q,c]-Toutref[q])^2+dd[951]*(xout[q,c]-xoutref[q])^2+dd[952]*(yout[q,c]-youtref[q])^2+dd[953]*(Wreb[q,c]-Wrebref[q])^2+dd[954]*(Wrebl[q,c]-Wreblref[q])^2+dd[955]*(Wrebv[q,c]-Wrebvref[q])^2+dd[956]*(roreb[q,c]-rorebref[q])^2+dd[957]*(rorebm[q,c]-rorebmref[q])^2+dd[958]*(M[1,q,c]-Mref[1,q])^2+dd[959]*(M[Ntray,q,c]-Mref[Ntray,q])^2+dd[960]*(L[1,q,c]-Lref[1,q])^2+dd[962]*(y[Ntray,q,c]-yref[Ntray,q])^2+dd[963]*(Qr[q,c]-Qrref[q])^2+dd[964]*(Qc[q,c]-Qcref[q])^2))*omega[c,ncp]+beta*(dd[945]*(D[q]-Dref[q])^2+dd[946]*(Psteam[q]-Psteamref[q])^2+dd[961]*(Reflux[q]-Refluxref[q])^2+dd[965]*(B[q]-Bref[q])^2));
var cost{q in fe}=h[q]*sum{c in cp}((alpha*Qr[q,c]+1e8*x1eps[q,c]+1e8*xNtrayeps[q,c])*omega[c,ncp]);
#var deltaB=sum{q in fe diff{1}}(1e3*(B[q]-B[q-1]));
minimize profit: sum{q in fe} (cost[q]+h[q]*sum{c in cp} (beta*sum{i in tray diff{1,Ntray}}(dd[i-1]*(M[i,q,c]-Mref[i,q])^2)+beta*sum{j in tray}(dd[j+156]*(T[j,q,c]-Trefe[j,q])^2) + beta*sum{k in tray diff{Ntray}}(dd[k+314]*(V[k,q,c]-Vref[k,q])^2)+beta*sum{l in tray diff{1}}(dd[l+470]*(L[l,q,c]-Lref[l,q])^2) + beta*sum{m in tray}(dd[m+629]*(x[m,q,c]-xref[m,q])^2)+ beta*sum{n in tray}(dd[n+787]*(y[n,q,c]-yref[n,q])^2)+beta*dd[629]*(Rec[q,c]-Recref[q])^2+beta*dd[947]*(Ts[q,c]-Tsref[q])^2+beta*dd[948]*(Fs[q,c]-Fsref[q])^2+beta*dd[949]*(Pchest[q,c]-Pchestref[q])^2+beta*dd[950]*(Tout[q,c]-Toutref[q])^2+beta*dd[951]*(xout[q,c]-xoutref[q])^2+beta*dd[952]*(yout[q,c]-youtref[q])^2+beta*dd[953]*(Wreb[q,c]-Wrebref[q])^2+beta*dd[954]*(Wrebl[q,c]-Wreblref[q])^2+beta*dd[955]*(Wrebv[q,c]-Wrebvref[q])^2+beta*dd[956]*(roreb[q,c]-rorebref[q])^2+beta*dd[957]*(rorebm[q,c]-rorebmref[q])^2+beta*dd[958]*(M[1,q,c]-Mref[1,q])^2+beta*dd[959]*(M[Ntray,q,c]-Mref[Ntray,q])^2+beta*dd[960]*(L[1,q,c]-Lref[1,q])^2+beta*dd[962]*(y[Ntray,q,c]-yref[Ntray,q])^2+beta*dd[963]*(Qr[q,c]-Qrref[q])^2+beta*dd[964]*(Qc[q,c]-Qcref[q])^2)*omega[c,ncp]+beta*dd[945]*(D[q]-Dref[q])^2+beta*dd[946]*(Psteam[q]-Psteamref[q])^2+beta*dd[961]*(Reflux[q]-Refluxref[q])^2+beta*dd[965]*(B[q]-Bref[q])^2);
#minimize profit: sum{q in fe} (cost[q]+h[q]*sum{c in cp} (beta*sum{i in tray diff{1,Ntray}}(dd[i-1]*(M[i,q,c]-Mref[i,q])^2)+beta*sum{j in tray}(dd[j+156]*(T[j,q,c]-Trefe[j,q])^2) + beta*sum{k in tray diff{Ntray}}(dd[k+314]*(V[k,q,c]-Vref[k,q])^2)+beta*sum{l in tray diff{1}}(dd[l+470]*(L[l,q,c]-Lref[l,q])^2) + beta*sum{m in tray}(dd[m+629]*(x[m,q,c]-xref[m,q])^2)+ beta*sum{n in tray}(dd[n+787]*(y[n,q,c]-yref[n,q])^2)+beta*dd[629]*(Rec[q,c]-Recref[q])^2+beta*dd[947]*(Ts[q,c]-Tsref[q])^2+beta*dd[948]*(Fs[q,c]-Fsref[q])^2+beta*dd[949]*(Pchest[q,c]-Pchestref[q])^2+beta*dd[950]*(Tout[q,c]-Toutref[q])^2+beta*dd[951]*(xout[q,c]-xoutref[q])^2+beta*dd[952]*(yout[q,c]-youtref[q])^2+beta*dd[953]*(Wreb[q,c]-Wrebref[q])^2+beta*dd[954]*(Wrebl[q,c]-Wreblref[q])^2+beta*dd[955]*(Wrebv[q,c]-Wrebvref[q])^2+beta*dd[956]*(roreb[q,c]-rorebref[q])^2+beta*dd[957]*(rorebm[q,c]-rorebmref[q])^2+beta*dd[958]*(M[1,q,c]-Mref[1,q])^2+beta*dd[959]*(M[Ntray,q,c]-Mref[Ntray,q])^2+beta*dd[960]*(L[1,q,c]-Lref[1,q])^2+beta*dd[962]*(y[Ntray,q,c]-yref[Ntray,q])^2+beta*dd[963]*(Qr[q,c]-Qrref[q])^2+beta*dd[964]*(Qc[q,c]-Qcref[q])^2)*omega[c,ncp]+beta*dd[945]*(D[q]-Dref[q])^2+beta*dd[946]*(Psteam[q]-Psteamref[q])^2+beta*dd[961]*(Reflux[q]-Refluxref[q])^2+beta*dd[965]*(B[q]-Bref[q])^2);

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

