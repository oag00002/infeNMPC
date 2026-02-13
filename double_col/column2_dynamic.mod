##########################
# The distillation model #
##########################

### Two columns in series
### Three components: A (light), B and C (heavy).
### Based on Column A. Steady state model.
### made by Roald Brck Leer
### Revised by Xue Yang into dynamic model

#--------------------------------------------
# Setting parameters as variables
#--------------------------------------------
param ncp:=3;
param nfe;
set fe:=1..nfe;
set cp:=1..ncp;

param nominal_F;
param nominal_pV;
param nominal_qF;


#-----------------
# Parameters
#-----------------

param NT := 41 ; #stages
param NC := 3  ; #number of components
param NF := 21 ; #feed enter at stage 21
param h{fe};
param time;
param omega{cp,cp};

param Fref;
param qFref;
param y1ref{1..NT-1,1..2};
param y2ref{1..NT-1,1..2};
param x1ref{1..NT,1..2};
param x2ref{1..NT,1..2};
param M1ref{1..NT};
param M2ref{1..NT};
param V1ref{1..NT-1};
param V2ref{1..NT-1};
param L1ref{2..NT};
param L2ref{2..NT};
param D1ref;
param B1ref;
param D2ref;
param B2ref;
param y_1_1ref{1..NT-1,1..2};
param y_1_2ref{1..NT-1,1..2};
param y_2_1ref{1..NT-1,1..2};
param y_2_2ref{1..NT-1,1..2};
param pVref;
param VB1ref;
param VB2ref;
param LT1ref;
param LT2ref;
param TC1ref{1..NT};
param TC2ref{1..NT};
param a;
param b;

param F_w;
param qF_w;
param y1_w{1..NT-1,1..2};
param y2_w{1..NT-1,1..2};
param x1_w{1..NT,1..2};
param x2_w{1..NT,1..2};
param M1_w{1..NT};
param M2_w{1..NT};
param V1_w{1..NT-1};
param V2_w{1..NT-1};
param L1_w{2..NT};
param L2_w{2..NT};
param D1_w;
param B1_w;
param D2_w;
param B2_w;
param y_1_1_w{1..NT-1,1..2};
param y_1_2_w{1..NT-1,1..2};
param y_2_1_w{1..NT-1,1..2};
param y_2_2_w{1..NT-1,1..2};
param pV_w;
param VB1_w;
param VB2_w;
param LT1_w;
param LT2_w;
param TC1_w{1..NT};
param TC2_w{1..NT}; 

# Disturbances
var F{fe} >= 0.1 ; #feed
var qF{fe} ; #liquid fraction in feed
var pV{fe} ; #energy price


### Relative volatility, alpha(A/B) alpha(B/C)
#relative volatility different from paper?
param alpha {1..2} ;
    let alpha[1] := 2.00 ;
    let alpha[2] := 1.50 ;

### Diagonal matrix of relative volatilities
param m_alpha {1..2,1..2} ;
for {i in 1..2, j in 1..2}
    let m_alpha[i,j] := 0 ;
for {i in 1..2}
    let m_alpha[i,i] := alpha[i] ;

### Parameters for Franci's Weir Formula L(i) = K**Mow(i)^1.5
param Kuf := 21.65032 ; # constant above feed
param Kbf := 29.65032 ; # constant below feed
param Muw := 0.25     ; # liquid holdup under weir (kmol)

### P-Controllers for control of reboiler and condenser hold up
param KcB, := 10  ; # controller gain reboiler
param KcD, := 10  ; # controller gain condenser
param MDs, := 0.5 ; # nominal holdup condenser
param MBs, := 0.5 ; # nominal holdup reboiler
param Ds , := 0.5 ; # nominal flow condenser
param Bs , := 0.5 ; # nominal flow reboiler
param M0 , := 0.5 ; # nominal liquid holdups

### Define prices
param pF, := 1       ; # feed price ($)
param pA, := 1       ; # light comp price ($)
param pB, := 2       ; # medium comp price ($)
param pC, := 1       ; # heavy comp price ($)
param zF {1..1,1..3} ; # feed composition
let zF[1,1] := 0.4 ;#
let zF[1,2] := 0.2 ;#
let zF[1,3] := 0.4 ;

#------------
# Variables
#------------

### FIRST COLUMN ----------------------------------------------------------
    var LT1{fe} >= 0 <= 10 := 3.43656 ;
    var VB1{fe} >= 0 <= 10 := 4.008 ;

### comp frac liquid and vap
    var y1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} >= 0 := 0.3 ; # vapor comp
    var x1 {i in 1..NT, j in 1..2, f in fe, c in cp} >= 0 :=0.4    ; # liquid comp
    var M1 {i in 1..NT, f in fe, c in cp} >= 0 := 0.5              ; # holdup
    var V1 {1..NT-1,fe,cp} >= 0 := 1                   ; # vapor flow
    var L1 {2..NT, fe,cp} >= 0 := 1                     ; # liquid flow
    var D1{fe,cp} >= 0 := 0.57                          ; # distillate flow
    var B1{fe,cp} >= 0 := 0.83                          ; # bottoms flow

# VLE equation split
    var y_1_1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} >= 0 := 0.2 ;
    var y_1_2 {i in 1..NT-1, j in 1..2, f in fe, c in cp} >= 0 := 0.7 ;
   #TC1, SET1 and SET2 are not the original model
    var TC1 {i in 1..NT, f in fe, c in cp} >= 350 <= 412 := 380;
    #var SET1 :=-150;
    #var SET2 :=-100;
    #var SET3 :=-270;

### SECOND COLUMN----------------------------------------------------------
    var LT2{fe} >= 0 <= 10 := 2.13827 ;
    var VB2{fe} >= 0 <= 10 := 2.40367 ;
    
### comp frac liquid and gas
    var y2 {i in 1..NT-1, j in 1..2, f in fe, c in cp} >= 0 := 0.5 ; # vapor comp
    var x2 {i in 1..NT, j in 1..2, f in fe, c in cp} >= 0 := 0.4 ;   # liquid comp
    var M2 {i in 1..NT, f in fe, c in cp} >= 0 := 0.5 ;              # holdup
    var V2 {1..NT-1, fe,cp} >= 0 :=1;                     # vapor flow
    var L2 {2..NT,fe,cp} >= 0 :=1;                       # liquid flow
    var D2{fe,cp} >= 0 := 0.26 ;                          # distillate flow
    var B2{fe,cp} >= 0 := 0.56 ;                          # bottoms flow
    
# VLE equation split
    var y_2_1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} := 0.3 ;
    var y_2_2 {i in 1..NT-1, j in 1..2, f in fe, c in cp} := 0.8 ;
    #TC2 is not the in original model
    var TC2 {i in 1..NT,f in fe, c in cp} >= 350 <= 412 := 380;

#----------------------------
# Model/Constraints 1
#----------------------------

### VLE equation split

const1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} :
    y_1_1[i,j,f,c] = x1[i,j,f,c]*m_alpha[j,j] ;

const2 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:
    y_1_2[i,j,f,c] = ((x1[i,1,f,c]*(alpha[1]-1)+x1[i,2,f,c]*(alpha[2]-1))+1) ;
    
### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const3 {i in 1..NT-1, j in 1..2,f in fe, c in cp}:
    y1[i,j,f,c] = y_1_1[i,j,f,c]/y_1_2[i,j,f,c] ;
    
### Vapor flows assuming constant molar flows##############################

const4 {i in 1..NF-1,f in fe, c in cp}:
    V1[i,f,c] = VB1[f] ; # vapor flow below feed

const5 {i in NF..NT-1,f in fe, c in cp}:
    V1[i,f,c] = VB1[f] + (1-qF[f])*F[f] ; # vapor flow above feed
    
### Liquid flows are given by Franci's Weir Formula L(i)=K*Mow(i)^1.5######
# Liquid flow L(i) dependent only on the holdup over the weir Mow(i) ######
#M(i)= Mow(i) + Muw(i) (Total holdup = holdup over weir +            ######
#holdup below weir)                                                  ######

const6 {i in 2..NF, f in fe, c in cp}:
    L1[i,f,c] = Kbf*(M1[i,f,c] - Muw)^1.5 ; # Liquid flow below feed

const7 {i in NF+1..NT-1,f in fe, c in cp}:
    L1[i,f,c] = Kuf*(M1[i,f,c] - Muw)^1.5 ;# Liquid flows above feed

const8{f in fe, c in cp}:
    L1[NT,f,c] = LT1[f]; # Condenser's liquid flow

### Distillate and bottom

const9{f in fe, c in cp}:
    B1[f,c] = Bs + (M1[1,f,c] - MBs)*KcB ;

const10{f in fe, c in cp}:
    D1[f,c] = Ds + (M1[NT,f,c] - MDs)*KcD ;

### Material balances for total holdup and component holdup################
var M1dot{1..NT,fe,cp};
const11 {i in 2..NF-1,f in fe, c in cp}:
    M1dot[i,f,c]=L1[i+1,f,c] - L1[i,f,c] + V1[i-1,f,c] - V1[i,f,c]; # dM below feed

const12 {i in NF+1..NT-1,f in fe, c in cp}:
    M1dot[i,f,c]=L1[i+1,f,c] - L1[i,f,c] + V1[i-1,f,c] - V1[i,f,c]; # dM above feed

var Mx1dot{1..NT,1..2,fe,cp};
const13 {i in 2..NF-1, j in 1..2,f in fe, c in cp}:
    Mx1dot[i,j,f,c]=L1[i+1,f,c]*x1[i+1,j,f,c] - L1[i,f,c]*x1[i,j,f,c] + V1[i-1,f,c]*y1[i-1,j,f,c] - V1[i,f,c]*y1[i,j,f,c];

const14 {i in NF+1..NT-1, j in 1..2,f in fe, c in cp}:
Mx1dot[i,j,f,c]=L1[i+1,f,c]*x1[i+1,j,f,c] - L1[i,f,c]*x1[i,j,f,c] + V1[i-1,f,c]*y1[i-1,j,f,c] - V1[i,f,c]*y1[i,j,f,c];

### Correction for feed at the feed stage: The feed is assumed to
#be mixed into the feed stage

const15{f in fe, c in cp}:
    M1dot[NF,f,c]=L1[NF+1,f,c] - L1[NF,f,c] + V1[NF-1,f,c] - V1[NF,f,c] + F[f];

const16{f in fe, c in cp}:
    Mx1dot[NF,1,f,c]=L1[NF+1,f,c]*x1[NF+1,1,f,c] - L1[NF,f,c]*x1[NF,1,f,c] + V1[NF-1,f,c]*y1[NF-1,1,f,c] - V1[NF,f,c]*y1[NF,1,f,c] + F[f]*zF[1,1];

const17{f in fe, c in cp}:
    Mx1dot[NF,2,f,c]=L1[NF+1,f,c]*x1[NF+1,2,f,c] - L1[NF,f,c]*x1[NF,2,f,c] + V1[NF-1,f,c]*y1[NF-1,2,f,c] - V1[NF,f,c]*y1[NF,2,f,c] + F[f]*zF[1,2];

### Reboiler (assumed to be an equilibrium stage)

const18{f in fe, c in cp}:
    M1dot[1,f,c]=L1[2,f,c] - V1[1,f,c] - B1[f,c];
    
const19 {j in 1..2,f in fe, c in cp}:
    Mx1dot[1,j,f,c]=L1[2,f,c]*x1[2,j,f,c] - V1[1,f,c]*y1[1,j,f,c] - B1[f,c]*x1[1,j,f,c];

### Total condenser (no equilibrium stage)

const20{f in fe, c in cp}:
    M1dot[NT,f,c]=V1[NT-1,f,c] - LT1[f] - D1[f,c];

const21 {j in 1..2,f in fe, c in cp}:
    Mx1dot[NT,j,f,c]=V1[NT-1,f,c]*y1[NT-1,j,f,c] - L1[NT,f,c]*x1[NT,j,f,c] - D1[f,c]*x1[NT,j,f,c];

var x1dot{i in 1..NT, j in 1..2,f in fe, c in cp}=(Mx1dot[i,j,f,c]-x1[i,j,f,c]*M1dot[i,f,c])/M1[i,f,c];

param x1set{1..NT,1..2};
param M1set{1..NT};
var x1_nom{1..NT,1..2};
var M1_nom{1..NT};
var x1_0{1..NT,1..2,fe};
var M1_0{1..NT,fe};

x1_init_sens{i in 1..NT,j in 1..2}: x1_nom[i,j]=x1set[i,j];
x1_init_constr{i in 1..NT, j in 1..2}:x1_0[i,j,1]=x1_nom[i,j];
x1_constr{i in 1..NT, j in 1..2,f in fe diff{1}}:x1_0[i,j,f]=x1[i,j,f-1,ncp];
x1_Lag{i in 1..NT,j in 1..2,q in fe, c in cp}: x1[i,j,q,c] = x1_0[i,j,q] + h[q]* sum{k in cp} omega[k,c]*x1dot[i,j,q,k];


M1_init_sens{i in 1..NT}: M1_nom[i]=M1set[i];
M1_init_constr{i in 1..NT}:M1_0[i,1]=M1_nom[i];
M1_constr{i in 1..NT,f in fe diff{1}}:M1_0[i,f]=M1[i,f-1,ncp];
M1_Lag{i in 1..NT,q in fe, c in cp}: M1[i,q,c] = M1_0[i,q] + h[q]* sum{k in cp} omega[k,c]*M1dot[i,q,k];


#----------------------------
# Model/Constraints 2
#----------------------------

### VLE equation split

const22 {i in 1..NT-1, j in 1..2, f in fe, c in cp} :
    y_2_1[i,j,f,c] = x2[i,j,f,c]*m_alpha[j,j] ;
    
const23 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:
    y_2_2[i,j,f,c] = ((x2[i,1,f,c]*(alpha[1]-1)+x2[i,2,f,c]*(alpha[2]-1))+1) ;

### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const24 {i in 1..NT-1, j in 1..2,f in fe, c in cp}:
    y2[i,j,f,c] = y_2_1[i,j,f,c]/y_2_2[i,j,f,c] ;
    
### Vapor flows assuming constant molar flows

const25 {i in 1..NF-1,f in fe, c in cp}:
    V2[i,f,c] = VB2[f] ; # vapor flow below feed

const26 {i in NF..NT-1,f in fe, c in cp}:
    V2[i,f,c] = VB2[f] ; # vapor flow above

const27 {i in 2..NF,f in fe, c in cp}:
    L2[i,f,c] = Kbf*(M2[i,f,c] - Muw)^1.5 ; # Liquid flow below feed

const28 {i in NF+1..NT-1,f in fe, c in cp}:
    L2[i,f,c] = Kuf*(M2[i,f,c] - Muw)^1.5 ; # Liquid flow above feed

const29{f in fe, c in cp}:
    L2[NT,f,c] = LT2[f]; # Condenser's liquid flow

const30{f in fe, c in cp}:
    B2[f,c] = Bs + (M2[1,f,c] - MBs)*KcB ;

const31{f in fe, c in cp}:
    D2[f,c] = Ds + (M2[NT,f,c] - MDs)*KcD ;

### Material balances for total holdup and component holdup
var M2dot{1..NT,fe,cp};
const32 {i in 2..NF-1,f in fe, c in cp}:
    M2dot[i,f,c]=L2[i+1,f,c] - L2[i,f,c] + V2[i-1,f,c] - V2[i,f,c]; # dM below feed

const33 {i in NF+1..NT-1,f in fe, c in cp}:
    M2dot[i,f,c]=L2[i+1,f,c] - L2[i,f,c] + V2[i-1,f,c] - V2[i,f,c]; # dM above feed

var Mx2dot{1..NT,1..2,fe,cp};
const34 {i in 2..NF-1, j in 1..2,f in fe, c in cp}: # dMxdt below feed
    Mx2dot[i,j,f,c]=L2[i+1,f,c]*x2[i+1,j,f,c] - L2[i,f,c]*x2[i,j,f,c] + V2[i-1,f,c]*y2[i-1,j,f,c] - V2[i,f,c]*y2[i,j,f,c];

const35 {i in NF+1..NT-1, j in 1..2,f in fe, c in cp}: # dMxdt above feed
    Mx2dot[i,j,f,c]=L2[i+1,f,c]*x2[i+1,j,f,c] - L2[i,f,c]*x2[i,j,f,c] + V2[i-1,f,c]*y2[i-1,j,f,c] - V2[i,f,c]*y2[i,j,f,c];

### Correction for feed at the feed stage: The feed is assumed to be
# mixed into the feed stage

const36{f in fe, c in cp}:
    M2dot[NF,f,c]=L2[NF+1,f,c] - L2[NF,f,c] + V2[NF-1,f,c] - V2[NF,f,c] + B1[f,c];

const37{f in fe, c in cp}:
    Mx2dot[NF,1,f,c]=L2[NF+1,f,c]*x2[NF+1,1,f,c] - L2[NF,f,c]*x2[NF,1,f,c] + V2[NF-1,f,c]*y2[NF-1,1,f,c] - V2[NF,f,c]*y2[NF,1,f,c] + B1[f,c]*x1[1,1,f,c];

const38{f in fe, c in cp}:
    Mx2dot[NF,2,f,c]=L2[NF+1,f,c]*x2[NF+1,2,f,c] - L2[NF,f,c]*x2[NF,2,f,c] + V2[NF-1,f,c]*y2[NF-1,2,f,c] - V2[NF,f,c]*y2[NF,2,f,c] + B1[f,c]*x1[1,2,f,c];

### Reboiler (assumed to be an equilibrium stage)

const39{f in fe, c in cp}:
M2dot[1,f,c]=L2[2,f,c] - V2[1,f,c] - B2[f,c];

const40 {i in 1..2,f in fe, c in cp}:
Mx2dot[1,i,f,c]=L2[2,f,c]*x2[2,i,f,c] - V2[1,f,c]*y2[1,i,f,c] - B2[f,c]*x2[1,i,f,c];

### Total condenser (no equilibrium stage)

const41{f in fe, c in cp}:
    M2dot[NT,f,c]=V2[NT-1,f,c] - LT2[f] - D2[f,c];

const42 {i in 1..2,f in fe, c in cp}:
    Mx2dot[NT,i,f,c]=V2[NT-1,f,c]*y2[NT-1,i,f,c] - L2[NT,f,c]*x2[NT,i,f,c] - D2[f,c]*x2[NT,i,f,c];

var x2dot{i in 1..NT, j in 1..2,f in fe, c in cp}=(Mx2dot[i,j,f,c]-x2[i,j,f,c]*M2dot[i,f,c])/M2[i,f,c];

param x2set{1..NT,1..2};
param M2set{1..NT};
var x2_nom{1..NT,1..2};
var M2_nom{1..NT};
var x2_0{1..NT,1..2,fe};
var M2_0{1..NT,fe};

x2_init_sens{i in 1..NT,j in 1..2}: x2_nom[i,j]=x2set[i,j];
x2_init_constr{i in 1..NT, j in 1..2}:x2_0[i,j,1]=x2_nom[i,j];
x2_constr{i in 1..NT, j in 1..2,f in fe diff{1}}:x2_0[i,j,f]=x2[i,j,f-1,ncp];
x2_Lag{i in 1..NT,j in 1..2,q in fe, c in cp}: x2[i,j,q,c] = x2_0[i,j,q] + h[q]* sum{k in cp} omega[k,c]*x2dot[i,j,q,k];


M2_init_sens{i in 1..NT}: M2_nom[i]=M2set[i];
M2_init_constr{i in 1..NT}:M2_0[i,1]=M2_nom[i];
M2_constr{i in 1..NT,f in fe diff{1}}:M2_0[i,f]=M2[i,f-1,ncp];
M2_Lag{i in 1..NT,q in fe, c in cp}: M2[i,q,c] = M2_0[i,q] + h[q]* sum{k in cp} omega[k,c]*M2dot[i,q,k];

#const43 and const44 are not in original model
const43 {i in 1..NT,f in fe, c in cp}:
     TC1[i,f,c] = x1[i,1,f,c]*353.3 + x1[i,2,f,c]*383.8 + (1-x1[i,1,f,c]-x1[i,2,f,c])*411.5 ;

const44 {f in fe, c in cp,i in 1..NT}:
     TC2[i,f,c] = x2[i,1,f,c]*353.3 + x2[i,2,f,c]*383.8 + (1-x2[i,1,f,c]-x2[i,2,f,c])*411.5 ; #x1[14,1]*353.3+x1[14,2]*383.8+(1-x1[14,1]-x1[14,2])*411.5 = 382.9450;

### Inequality constraints

const45{f in fe, c in cp}:
    x1[NT,1,f,c] >= 0.95 ; # xA

const46{f in fe, c in cp}:
    x2[NT,2,f,c] >= 0.95 ; # xB (= 0.95)

const47{f in fe, c in cp}:
    1 - (x2[1,1,f,c]+x2[1,2,f,c]) >= 0.95 ; #

const48{f in fe}:
    VB1[f] <= 4.008 ; # max boilup column 1

const49{f in fe}:
    VB2[f] <= 2.405 ; # max boilup column 2

### Initial constraints for parameters

constF{f in fe}: F[f] = nominal_F ;
constpV{f in fe}: pV[f] = nominal_pV ;
# constzF: zF = nominal_zF ;
constqF{f in fe}: qF[f] = nominal_qF ;

#--------------------
# OBJECTIVE FUNCTION
#--------------------
#var deltau=sum{i in fe diff{fe}}((LT1[i+1]-LT1[i])^2+(VB1[i+1]-VB1[i])^2+(LT2[i+1]-LT2[i])^2+(VB2[i+1]-VB2[i])^2); 
var cost{f in fe}=(pF*F[f] + pV[f]*(VB1[f] + VB2[f])+h[f]*sum{c in cp}(-pA*D1[f,c] - pB*D2[f,c] - pC*B2[f,c])*omega[c,ncp]);

var regu_1=sum{f in fe} ((F[f]-Fref)^2 +(qF[f]-qFref)^2+(pV[f]-pVref)^2+(VB1[f]-VB1ref)^2 +(VB2[f]-VB2ref)^2+(LT1[f]-LT1ref)^2 +(LT2[f]-LT2ref)^2);
var regu_2=sum{f in fe}(h[f]*sum{c in cp}((D1[f,c]-D1ref)^2 +(B1[f,c]-B1ref)^2+(D2[f,c]-D2ref)^2 +(B2[f,c]-B2ref)^2)*omega[c,ncp]);
var regu_M_TC=sum{i in 1..NT}(sum{f in fe}(h[f]*sum{c in cp}((M1[i,f,c]-M1ref[i])^2 +(M2[i,f,c]-M2ref[i])^2+(TC1[i,f,c]-TC1ref[i])^2 +(TC2[i,f,c]-TC2ref[i])^2)*omega[c,ncp]));
var regu_x_y=sum{i in 1..NT}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}((x1[i,j,f,c]-x1ref[i,j])^2 +(x2[i,j,f,c]-x2ref[i,j])^2)*omega[c,ncp])))+sum{i in 1..NT-1}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}((y1[i,j,f,c]-y1ref[i,j])^2 +(y2[i,j,f,c]-y2ref[i,j])^2)*omega[c,ncp])));
var regu_V_L=sum{i in 1..NT-1}(sum{f in fe}(h[f]*sum{c in cp}((V1[i,f,c]-V1ref[i])^2 +(V2[i,f,c]-V2ref[i])^2)*omega[c,ncp]))+sum{i in 2..NT}(sum{f in fe}(h[f]*sum{c in cp}((L1[i,f,c]-L1ref[i])^2 +(L2[i,f,c]-L2ref[i])^2)*omega[c,ncp]));
var regu_y=sum{i in 1..NT-1}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}((y_1_1[i,j,f,c]-y_1_1ref[i,j])^2 +(y_1_2[i,j,f,c]-y_1_2ref[i,j])^2+(y_2_1[i,j,f,c]-y_2_1ref[i,j])^2 +(y_2_2[i,j,f,c]-y_2_2ref[i,j])^2)*omega[c,ncp])));
minimize obj: a*sum{f in 1..nfe}cost[f]+b*(regu_1+regu_2+regu_M_TC+regu_x_y+regu_V_L+regu_y);

#var regu_1=sum{f in fe} (F_w*(F[f]-Fref)^2 +qF_w*(qF[f]-qFref)^2+pV_w*(pV[f]-pVref)^2+VB1_w*(VB1[f]-VB1ref)^2 +VB2_w*(VB2[f]-VB2ref)^2+LT1_w*(LT1[f]-LT1ref)^2 +LT2_w*(LT2[f]-LT2ref)^2);
#var regu_2=sum{f in fe}(h[f]*sum{c in cp}(D1_w*(D1[f,c]-D1ref)^2 +B1_w*(B1[f,c]-B1ref)^2+D2_w*(D2[f,c]-D2ref)^2 +B2_w*(B2[f,c]-B2ref)^2)*omega[c,ncp]);
#var regu_M_TC=sum{i in 1..NT}(sum{f in fe}(h[f]*sum{c in cp}(M1_w[i]*(M1[i,f,c]-M1ref[i])^2 +M2_w[i]*(M2[i,f,c]-M2ref[i])^2+TC1_w[i]*(TC1[i,f,c]-TC1ref[i])^2 +TC2_w[i]*(TC2[i,f,c]-TC2ref[i])^2)*omega[c,ncp]));
#var regu_x_y=sum{i in 1..NT}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}(x1_w[i,j]*(x1[i,j,f,c]-x1ref[i,j])^2 +x2_w[i,j]*(x2[i,j,f,c]-x2ref[i,j])^2)*omega[c,ncp])))+sum{i in 1..NT-1}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}(y1_w[i,j]*(y1[i,j,f,c]-y1ref[i,j])^2 +y2_w[i,j]*(y2[i,j,f,c]-y2ref[i,j])^2)*omega[c,ncp])));
#var regu_V_L=sum{i in 1..NT-1}(sum{f in fe}(h[f]*sum{c in cp}(V1_w[i]*(V1[i,f,c]-V1ref[i])^2 +V2_w[i]*(V2[i,f,c]-V2ref[i])^2)*omega[c,ncp]))+sum{i in 2..NT}(sum{f in fe}(h[f]*sum{c in cp}(L1_w[i]*(L1[i,f,c]-L1ref[i])^2 +L2_w[i]*(L2[i,f,c]-L2ref[i])^2)*omega[c,ncp]));
#var regu_y=sum{i in 1..NT-1}(sum{j in 1..2}(sum{f in fe}(h[f]*sum{c in cp}(y_1_1_w[i,j]*(y_1_1[i,j,f,c]-y_1_1ref[i,j])^2 +y_1_2_w[i,j]*(y_1_2[i,j,f,c]-y_1_2ref[i,j])^2+y_2_1_w[i,j]*(y_2_1[i,j,f,c]-y_2_1ref[i,j])^2 +y_2_2_w[i,j]*(y_2_2[i,j,f,c]-y_2_2ref[i,j])^2)*omega[c,ncp])));
#minimize obj: sum{f in 1..nfe}cost[f]+(regu_1+regu_2+regu_M_TC+regu_x_y+regu_V_L+regu_y);
