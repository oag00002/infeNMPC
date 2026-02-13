##########################
# The distillation model #
##########################

### Two columns in series
### Three components: A (light), B and C (heavy).
### Based on Column A. Steady state model.
### made by Roald Brck Leer
### Revised by Xue Yang into dynamic model
## Edited by Devin Griffith


#--------------------------------------------
# Run options
#--------------------------------------------
param sigma := .01;
param K :=25;
param econ_switch := 0;
param track_switch := 1;
param state_noise_switch := 0 ; 
param pen_switch:=1;
param noise_switch := 0;
param regu_frac:= 0;
param rho := 10^4;
param smt := .01;
param F_dev:=.141;
param zF_dev:=.03;
param state_dev := 0.01 ;
param state_w := 10 ;
param cont_w := 1 ;
param rho_L := 1;
param rho_psi:=100000;
param termnorm ; 

#--------------------------------------------
# Setting parameters as variables
#--------------------------------------------
param ncp:=3;
param nfe;
param Nmax:=25;
set fe:=0..nfe-1;
set cp:=1..ncp;
set weights:=1..979 ;

param nominal_F;
param nominal_pV;
param nominal_qF;
param nominal_zF{1..3};
param weight {weights} ;



#-----------------
# Parameters
#-----------------

param time;
param NT := 41 ; #stages
param NC := 3  ; #number of components
param NF := 21 ; #feed enter at stage 21
param h;
param omega{cp,cp};
param numstat := 6*NT;
param state_act{1..numstat,0..K} ; #actual states in order of x11,x12,m1,x21,x22,m2
param xak{0..K};
param xbk{0..K};
param xck{0..K};
param state_act_ss{1..numstat};
param sol_stat{0..K-1};
param sens_sol_stat{0..K-1};
param step_cost{0..K-1};
param tot_cost;
param avg_cost;
param normx{0..K}; #norm of states
set s := {1..numstat} ; 
param P{s,s} ; 

param Nk{0..K} ;
param Ns := 5;
param normz{0..K-1,0..Nmax};
param normz_sens{0..K-1,0..Nmax};
param jj;
param termnorm_b2:=.41;
param termnorm_a2:=.38;
param termnorm_c2:=.36;
param term_step{0..K-1} ;
param sens_term_step{0..K-1};
param lecss := -.223;
var stat_vec{1..6*NT} ;
param zFk{0..K-1,1..3} ; 
param Fk{0..K-1} ;

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
param step_time{0..K-1} ;
param uk{1..8,0..K-1};

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
param F ; #feed
param  qF ; #liquid fraction in feed
param  pV ; #energy price
param zF {1..3} ; # feed composition


param F_r ;
param zF_r {1..3} ; # feed composition
param F_rand {0..K-1} ;
param zF_rand {0..K-1,1..3};



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
param M0  := 0.5 ; # nominal liquid holdups

### Define prices
param pF, := 1       ; # feed price ($)
param pA, := 1       ; # light comp price ($)
param pB, := 2       ; # medium comp price ($)
param pC, := 1       ; # heavy comp price ($)

#------------
# Variables
#------------

### vtr
var vk ;
var lk;

### FIRST COLUMN ----------------------------------------------------------
    var LT1{fe} >= 0 <= 10 := 3.43656 ;
    var VB1{fe} >= 0 <= 10 := 4.008 ;

### comp frac liquid and vap
    var y1 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:= 0.3 ; # vapor comp
    var x1 {i in 1..NT, j in 1..2, f in fe, c in cp}:=0.4    ; # liquid comp
    var x1_0{1..NT,1..2,0..nfe};
    var M1_0{1..NT,0..nfe};
    var M1 {i in 1..NT, f in fe, c in cp}:= 0.5              ; # holdup
    var V1 {1..NT-1,f in fe}:= 1                   ; # vapor flow
    var L1 {2..NT, f in fe, c in cp}:= 1                     ; # liquid flow
    var L1_0{2..NT,f in 0..nfe} ;
    var D1{fe}:= 0.57                          ; # distillate flow
    var B1{fe}:= 0.83                          ; # bottoms flow

    
#slack variables
    var x1eps {i in 1..NT, j in 1..2, 0..nfe} >= 0; 
    var M1eps {i in 1..NT, f in 0..nfe} >= 0           ; 
    var V1eps {1..NT-1,fe} >= 0    ; 
    var L1eps {2..NT, fe} >= 0                 ; 
    var D1eps{fe} >= 0                  ; 
    var B1eps{fe} >= 0                 ;
   	var TC1eps {i in 1..NT, f in 0..nfe} >= 0;

#soft constraints
x1lower {i in 1..NT, j in 1..2, f in 0..nfe}:x1_0[i,j,f]>=0-x1eps[i,j,f];
M1lower {i in 1..NT, f in 0..nfe}:M1_0[i,f]>=.25-M1eps[i,f];
M1upper {i in 1..NT, f in 0..nfe, c in cp}:M1_0[i,f]<=.75+M1eps[i,f];
V1lower {i in 1..NT-1,f in fe}:V1[i,f]>=0-V1eps[i,f];
L1lower {i in 2..NT, f in fe}: L1_0[i,f] >= 0-L1eps[i,f];
D1lower {f in fe}: D1[f]>=0-D1eps[f];
B1lower {f in fe}: B1[f]>=0-B1eps[f];
D1upper {f in fe}: D1[f] <= 2+D1eps[f] ; 
B1upper {f in fe}: B1[f] <= 2+B1eps[f] ; 

# VLE equation split
    var y_1_1 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:= 0.2 ;
    var y_1_2 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:= 0.7 ;
        var y_1_1_0 {i in 1..NT-1, j in 1..2, f in fe}=x1_0[i,j,f]*m_alpha[j,j];
    var y_1_2_0 {i in 1..NT-1, j in 1..2, f in fe}= ((x1_0[i,1,f]*(alpha[1]-1)+x1_0[i,2,f]*(alpha[2]-1))+1) ;    
      var y1_0 {i in 1..NT-1, j in 1..2, f in fe}= y_1_1_0[i,j,f] / y_1_2_0[i,j,f]  ; # vapor comp
   #TC1, SET1 and SET2 are not the original model
    var TC1 {i in 1..NT, f in 0..nfe}:= 380;
    #var SET1 :=-150;
    #var SET2 :=-100;
    #var SET3 :=-270;

#soft constraints
TC1lower {i in 1..NT, f in 0..nfe}: TC1[i,f]>=350-TC1eps[i,f];
TC1upper {i in 1..NT, f in 0..nfe}: TC1[i,f]<=412+TC1eps[i,f];

#real variables with uncertainty

    var y1_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.3 ; # vapor comp
    var x1_r {i in 1..NT, j in 1..2, c in cp}:=0.4    ; # liquid comp
    var M1_r {i in 1..NT, c in cp}:= 0.5              ; # holdup
    var L1_r {2..NT, c in cp}:= 1                     ; # liquid flow
    var y_1_1_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.2 ;
    var y_1_2_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.7 ;

### SECOND COLUMN----------------------------------------------------------
    var LT2{fe} >= 0 <= 10 := 2.13827 ;
    var VB2{fe} >= 0 <= 10 := 2.40367 ;
    
### comp frac liquid and gas
    var y2 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:= 0.5 ; # vapor comp

    var x2 {i in 1..NT, j in 1..2, f in fe, c in cp}:= 0.4 ;   # liquid comp
    var M2 {i in 1..NT, f in fe, c in cp}:= 0.5 ;              # holdup
    var x2_0{1..NT,1..2,0..nfe};
    var M2_0{1..NT,0..nfe};
    var V2 {1..NT-1, fe}:=1;                     # vapor flow
    var L2 {i in 2..NT,f in fe,c in cp}:=1;                       # liquid flow
    var L2_0 {i in 2..NT, f in 0..nfe} ; 
    var D2{fe}:= 0.26 ;                          # distillate flow
    var B2{fe}:= 0.56 ;                          # bottoms flow

#slack variables
    var x2eps {i in 1..NT, j in 1..2, f in 0..nfe} >= 0; 
    var x2ceps{0..nfe}>= 0; 
    var M2eps {i in 1..NT, f in 0..nfe} >= 0       ; 
    var V2eps {1..NT-1,fe} >= 0      ; 
    var L2eps {2..NT, fe} >= 0                 ; 
    var D2eps{fe} >= 0                       ; 
    var B2eps{fe} >= 0;       
    var TC2eps {i in 1..NT, f in 0..nfe} >= 0;      
    var termeps >= 0  ;      
    var desceps >= 0 ;          ;

#soft constraints
x2lower {i in 1..NT, j in 1..2, f in 0..nfe}:x2_0[i,j,f]>=0-x2eps[i,j,f];
M2lower {i in 1..NT, f in 0..nfe}:M2_0[i,f]>=.25-M2eps[i,f];
M2upper {i in 1..NT, f in 0..nfe}:M2_0[i,f]<=.75+M2eps[i,f];
V2lower {i in 1..NT-1,f in fe}:V2[i,f]>=0-V2eps[i,f];
L2lower {i in 2..NT, f in fe}: L2_0[i,f] >= 0-L2eps[i,f];
D2lower {f in fe}: D2[f]>=0-D2eps[f];
B2lower {f in fe}: B2[f]>=0-B2eps[f];
D2upper {f in fe}: D1[f] <= 2+D2eps[f] ; 
B2upper {f in fe}: B1[f] <= 2+B2eps[f] ; 
const45{f in 0..nfe}:  x1_0[NT,1,f] >= 0.95-x1eps[NT,1,f] ; # xA
const46{f in 0..nfe}:   x2_0[NT,2,f] >= 0.95-x2eps[NT,2,f] ; # xB (= 0.95)
const47{f in 0..nfe}:  1 - (x2_0[1,1,f]+x2_0[1,2,f]) >= 0.95-x2ceps[f] ; #
    
# VLE equation split
    var y_2_1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} := 0.3 ;
    var y_2_2 {i in 1..NT-1, j in 1..2, f in fe, c in cp} := 0.8 ;
        var y_2_1_0 {i in 1..NT-1, j in 1..2, f in fe} = x2_0[i,j,f]*m_alpha[j,j] ;
    var y_2_2_0 {i in 1..NT-1, j in 1..2, f in fe}  = ((x2_0[i,1,f]*(alpha[1]-1)+x2_0[i,2,f]*(alpha[2]-1))+1) ;
          var y2_0 {i in 1..NT-1, j in 1..2, f in fe}= y_2_1_0[i,j,f]/y_2_2_0[i,j,f]; # vapor comp
    #TC2 is not the in original model
    var TC2 {i in 1..NT,f in 0..nfe}:= 380;


#soft constraints
TC2lower {i in 1..NT, f in 0..nfe}: TC2[i,f]>=350-TC2eps[i,f];
TC2upper {i in 1..NT, f in 0..nfe}: TC2[i,f]<=412+TC2eps[i,f];

#real variables with uncertainty

    var y2_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.3 ; # vapor comp
    var x2_r {i in 1..NT, j in 1..2, c in cp}:=0.4    ; # liquid comp
    var M2_r {i in 1..NT, c in cp}:= 0.5              ; # holdup
    var L2_r {2..NT, c in cp}:= 1                     ; # liquid flow
    var y_2_1_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.2 ;
    var y_2_2_r {i in 1..NT-1, j in 1..2, c in cp}:= 0.7 ;


#----------------------------
# Model/Constraints 1
#----------------------------

### VLE equation split

const1 {i in 1..NT-1, j in 1..2, f in fe, c in cp} :
    y_1_1[i,j,f,c] = x1[i,j,f,c]*m_alpha[j,j] ;
    
    const1_r {i in 1..NT-1, j in 1..2, c in cp} :
    y_1_1_r[i,j,c] = x1_r[i,j,c]*m_alpha[j,j] ;
    

const2 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:
    y_1_2[i,j,f,c] = ((x1[i,1,f,c]*(alpha[1]-1)+x1[i,2,f,c]*(alpha[2]-1))+1) ;
    
    const2_r {i in 1..NT-1, j in 1..2, c in cp}:
    y_1_2_r[i,j,c] = ((x1_r[i,1,c]*(alpha[1]-1)+x1_r[i,2,c]*(alpha[2]-1))+1) ;
   
    
    
### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const3 {i in 1..NT-1, j in 1..2,f in fe, c in cp}:
    y1[i,j,f,c]*y_1_2[i,j,f,c] = y_1_1[i,j,f,c] ;
    
    const3_r {i in 1..NT-1, j in 1..2, c in cp}:
    y1_r[i,j,c]*y_1_2_r[i,j,c] = y_1_1_r[i,j,c] ;
    
    
    
### Vapor flows assuming constant molar flows##############################

const4 {i in 1..NF-1,f in fe}:
    V1[i,f] = VB1[f] ; # vapor flow below feed

const5 {i in NF..NT-1,f in fe}:
    V1[i,f] = VB1[f] + (1-qF)*F ; # vapor flow above feed
    
### Liquid flows are given by Franci's Weir Formula L(i)=K*Mow(i)^1.5######
# Liquid flow L(i) dependent only on the holdup over the weir Mow(i) ######
#M(i)= Mow(i) + Muw(i) (Total holdup = holdup over weir +            ######
#holdup below weir)                                                  ######

const6 {i in 2..NF, f in fe, c in cp}:
    L1[i,f,c] = Kbf*( (M1[i,f,c] - Muw+sqrt((M1[i,f,c] - Muw)^2+10^(-8)))/2    )^1.5 ; # Liquid flow below feed
    
    const6_r {i in 2..NF, c in cp}:
    L1_r[i,c] = Kbf*( (M1_r[i,c] - Muw+sqrt((M1_r[i,c] - Muw)^2+10^(-8)))/2    )^1.5 ; # Liquid flow below feed

const7 {i in NF+1..NT-1,f in fe, c in cp}:
    L1[i,f,c] = Kuf*( (M1[i,f,c] - Muw+sqrt((M1[i,f,c] - Muw)^2+10^(-8)))/2 )^1.5 ;# Liquid flows above feed
    
    const7_r {i in NF+1..NT-1,f in fe, c in cp}:
    L1_r[i,c] = Kuf*( (M1_r[i,c] - Muw+sqrt((M1_r[i,c] - Muw)^2+10^(-8)))/2 )^1.5 ;# Liquid flows above feed

const8{f in fe, c in cp}:
    L1[NT,f,c] = LT1[f]; # Condenser's liquid flow
    
    const8_r{ c in cp}:
    L1_r[NT,c] = LT1[0]; # Condenser's liquid flow


### Material balances for total holdup and component holdup################
var M1dot{1..NT,fe,cp};
const11 {i in 2..NF-1,f in fe, c in cp}:
    M1dot[i,f,c]=L1[i+1,f,c] - L1[i,f,c] + V1[i-1,f] - V1[i,f]; # dM below feed
    
    var M1dot_r{1..NT,cp};
const11_r {i in 2..NF-1, c in cp}:
    M1dot_r[i,c]=L1_r[i+1,c] - L1_r[i,c] + V1[i-1,0] - V1[i,0]; # dM below feed

const12 {i in NF+1..NT-1,f in fe, c in cp}:
    M1dot[i,f,c]=L1[i+1,f,c] - L1[i,f,c] + V1[i-1,f] - V1[i,f]; # dM above feed
    
    const12_r {i in NF+1..NT-1, c in cp}:
    M1dot_r[i,c]=L1_r[i+1,c] - L1_r[i,c] + V1[i-1,0] - V1[i,0]; # dM above feed

var Mx1dot{1..NT,1..2,fe,cp};
const13 {i in 2..NF-1, j in 1..2,f in fe, c in cp}:
Mx1dot[i,j,f,c]=L1[i+1,f,c]*x1[i+1,j,f,c] - L1[i,f,c]*x1[i,j,f,c] + V1[i-1,f]*y1[i-1,j,f,c] - V1[i,f]*y1[i,j,f,c];

var Mx1dot_r{1..NT,1..2,cp};
const13_r {i in 2..NF-1, j in 1..2, c in cp}:
Mx1dot_r[i,j,c]=L1_r[i+1,c]*x1_r[i+1,j,c] - L1_r[i,c]*x1_r[i,j,c] + V1[i-1,0]*y1_r[i-1,j,c] - V1[i,0]*y1_r[i,j,c];

const14 {i in NF+1..NT-1, j in 1..2,f in fe, c in cp}:
Mx1dot[i,j,f,c]=L1[i+1,f,c]*x1[i+1,j,f,c] - L1[i,f,c]*x1[i,j,f,c] + V1[i-1,f]*y1[i-1,j,f,c] - V1[i,f]*y1[i,j,f,c];

const14_r {i in NF+1..NT-1, j in 1..2, c in cp}:
Mx1dot_r[i,j,c]=L1_r[i+1,c]*x1_r[i+1,j,c] - L1_r[i,c]*x1_r[i,j,c] + V1[i-1,0]*y1_r[i-1,j,c] - V1[i,0]*y1_r[i,j,c];

### Correction for feed at the feed stage: The feed is assumed to
#be mixed into the feed stage

const15{f in fe, c in cp}:
    M1dot[NF,f,c]=L1[NF+1,f,c] - L1[NF,f,c] + V1[NF-1,f] - V1[NF,f] + F;
    
    const15_r{ c in cp}:
    M1dot_r[NF,c]=L1_r[NF+1,c] - L1_r[NF,c] + V1[NF-1,0] - V1[NF,0] + F_r;

const16{f in fe, c in cp}:
    Mx1dot[NF,1,f,c]=L1[NF+1,f,c]*x1[NF+1,1,f,c] - L1[NF,f,c]*x1[NF,1,f,c] + V1[NF-1,f]*y1[NF-1,1,f,c] - V1[NF,f]*y1[NF,1,f,c] + F*zF[1];
    
    const16_r{ c in cp}:
    Mx1dot_r[NF,1,c]=L1_r[NF+1,c]*x1_r[NF+1,1,c] - L1_r[NF,c]*x1_r[NF,1,c] + V1[NF-1,0]*y1_r[NF-1,1,c] - V1[NF,0]*y1_r[NF,1,c] + F_r*zF_r[1];

const17{f in fe, c in cp}:
    Mx1dot[NF,2,f,c]=L1[NF+1,f,c]*x1[NF+1,2,f,c] - L1[NF,f,c]*x1[NF,2,f,c] + V1[NF-1,f]*y1[NF-1,2,f,c] - V1[NF,f]*y1[NF,2,f,c] + F*zF[2];
    
    const17_r{ c in cp}:
    Mx1dot_r[NF,2,c]=L1_r[NF+1,c]*x1_r[NF+1,2,c] - L1_r[NF,c]*x1_r[NF,2,c] + V1[NF-1,0]*y1_r[NF-1,2,c] - V1[NF,0]*y1_r[NF,2,c] + F_r*zF_r[2];

### Reboiler (assumed to be an equilibrium stage)

const18{f in fe, c in cp}:
    M1dot[1,f,c]=L1[2,f,c] - V1[1,f] - B1[f];
    
const19 {j in 1..2,f in fe, c in cp}:
    Mx1dot[1,j,f,c]=L1[2,f,c]*x1[2,j,f,c] - V1[1,f]*y1[1,j,f,c] - B1[f]*x1[1,j,f,c];
    
    const18_r{c in cp}:
    M1dot_r[1,c]=L1_r[2,c] - V1[1,0] - B1[0];
    
const19_r {j in 1..2, c in cp}:
    Mx1dot_r[1,j,c]=L1_r[2,c]*x1_r[2,j,c] - V1[1,0]*y1_r[1,j,c] - B1[0]*x1_r[1,j,c];

### Total condenser (no equilibrium stage)

const20{f in fe, c in cp}:
    M1dot[NT,f,c]=V1[NT-1,f] - LT1[f] - D1[f];
    
    const20_r{c in cp}:
    M1dot_r[NT,c]=V1[NT-1,0] - LT1[0] - D1[0];

const21 {j in 1..2,f in fe, c in cp}:
    Mx1dot[NT,j,f,c]=V1[NT-1,f]*y1[NT-1,j,f,c] - L1[NT,f,c]*x1[NT,j,f,c] - D1[f]*x1[NT,j,f,c];
    
    const21_r {j in 1..2, c in cp}:
    Mx1dot_r[NT,j,c]=V1[NT-1,0]*y1_r[NT-1,j,c] - L1_r[NT,c]*x1_r[NT,j,c] - D1[0]*x1_r[NT,j,c];

var x1dot{1..NT,1..2,fe,cp};
var x1dot_r{1..NT,1..2,cp};

const21_2{i in 1..NT, j in 1..2,f in fe, c in cp}:
x1dot[i,j,f,c]*M1[i,f,c]=(Mx1dot[i,j,f,c]-x1[i,j,f,c]*M1dot[i,f,c]);

const21_2_r{i in 1..NT, j in 1..2, c in cp}:
x1dot_r[i,j,c]*M1_r[i,c]=(Mx1dot_r[i,j,c]-x1_r[i,j,c]*M1dot_r[i,c]);

param x1set{1..NT,1..2};
param M1set{1..NT};

x1_init_constr{i in 1..NT, j in 1..2}:x1_0[i,j,0]=x1set[i,j];
x1_constr{i in 1..NT, j in 1..2,f in 1..nfe}:x1_0[i,j,f]=x1[i,j,f-1,ncp];
x1_Lag{i in 1..NT,j in 1..2,q in fe , c in cp}: x1[i,j,q,c] = x1_0[i,j,q] + h* sum{k in cp} omega[k,c]*x1dot[i,j,q,k];
x1_r_Lag{i in 1..NT,j in 1..2,  c in cp}: x1_r[i,j,c] = x1_0[i,j,0] + h* sum{k in cp} omega[k,c]*x1dot_r[i,j,k];
L1_const{i in 2..NT, f in 1..nfe}:L1_0[i,f] = L1[i,f-1,ncp];
L1_const2 {i in 2..NF}:L1_0[i,0] = Kbf*( (M1_0[i,0] - Muw+sqrt((M1_0[i,0] - Muw)^2+10^(-8)))/2    )^1.5 ; # Liquid flow below feed
L1_const3 {i in NF+1..NT-1}: L1_0[i,0] = Kuf*( (M1_0[i,0] - Muw+sqrt((M1_0[i,0] - Muw)^2+10^(-8)))/2 )^1.5 ;# Liquid flows above feed
M1_init_constr{i in 1..NT}:M1_0[i,0]=M1set[i];
M1_constr{i in 1..NT,f in 1..nfe}:M1_0[i,f]=M1[i,f-1,ncp];
M1_Lag{i in 1..NT,q in fe, c in cp}: M1[i,q,c] = M1_0[i,q] + h* sum{k in cp} omega[k,c]*M1dot[i,q,k];
M1_r_Lag{i in 1..NT, c in cp}: M1_r[i,c] = M1_0[i,0] + h* sum{k in cp} omega[k,c]*M1dot_r[i,k];


#----------------------------
# Model/Constraints 2
#----------------------------

### VLE equation split

const22 {i in 1..NT-1, j in 1..2, f in fe, c in cp} :
    y_2_1[i,j,f,c] = x2[i,j,f,c]*m_alpha[j,j] ;
    
    const22_r {i in 1..NT-1, j in 1..2, c in cp} :
    y_2_1_r[i,j,c] = x2_r[i,j,c]*m_alpha[j,j] ;
   
    
const23 {i in 1..NT-1, j in 1..2, f in fe, c in cp}:
    y_2_2[i,j,f,c] = ((x2[i,1,f,c]*(alpha[1]-1)+x2[i,2,f,c]*(alpha[2]-1))+1) ;
    
    const23_r {i in 1..NT-1, j in 1..2, c in cp}:
    y_2_2_r[i,j,c] = ((x2_r[i,1,c]*(alpha[1]-1)+x2_r[i,2,c]*(alpha[2]-1))+1) ;
    

### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const24 {i in 1..NT-1, j in 1..2,f in fe, c in cp}:
    y2[i,j,f,c]*y_2_2[i,j,f,c] = y_2_1[i,j,f,c] ;
    
    const24_r {i in 1..NT-1, j in 1..2, c in cp}:
    y2_r[i,j,c]*y_2_2_r[i,j,c] = y_2_1_r[i,j,c] ;
    
    
### Vapor flows assuming constant molar flows

const25 {i in 1..NF-1,f in fe}:
    V2[i,f] = VB2[f] ; # vapor flow below feed

const26 {i in NF..NT-1,f in fe}:
    V2[i,f] = VB2[f] ; # vapor flow above

const27 {i in 2..NF,f in fe, c in cp}:
    L2[i,f,c] = Kbf*( (M2[i,f,c] - Muw+sqrt((M2[i,f,c] - Muw)^2+10^(-8)))/2 )^1.5 ; # Liquid flow below feed
    
    const27_r {i in 2..NF, c in cp}:
    L2_r[i,c] = Kbf*( (M2_r[i,c] - Muw+sqrt((M2_r[i,c] - Muw)^2+10^(-8)))/2 )^1.5 ; # Liquid flow below feed

const28 {i in NF+1..NT-1,f in fe, c in cp}:
    L2[i,f,c] = Kuf*( (M2[i,f,c] - Muw+sqrt((M2[i,f,c] - Muw)^2+10^(-8)))/2 )^1.5 ; # Liquid flow above feed
    
    const28_r {i in NF+1..NT-1, c in cp}:
    L2_r[i,c] = Kuf*( (M2_r[i,c] - Muw+sqrt((M2_r[i,c] - Muw)^2+10^(-8)))/2 )^1.5 ; # Liquid flow above feed

const29{f in fe,c in cp}:
    L2[NT,f,c] = LT2[f]; # Condenser's liquid flow
    
    const29_r{c in cp}:
    L2_r[NT,c] = LT2[0]; # Condenser's liquid flow

### Material balances for total holdup and component holdup
var M2dot{1..NT,fe,cp};
var M2dot_r{1..NT,cp};

const32 {i in 2..NF-1,f in fe, c in cp}:
    M2dot[i,f,c]=L2[i+1,f,c] - L2[i,f,c] + V2[i-1,f] - V2[i,f]; # dM below feed
    
    const32_r {i in 2..NF-1, c in cp}:
    M2dot_r [i,c]=L2_r[i+1,c] - L2_r[i,c] + V2[i-1,0] - V2[i,0]; # dM below feed

const33 {i in NF+1..NT-1,f in fe, c in cp}:
    M2dot[i,f,c]=L2[i+1,f,c] - L2[i,f,c] + V2[i-1,f] - V2[i,f]; # dM above feed
    
    const33_r {i in NF+1..NT-1, c in cp}:
    M2dot_r[i,c]=L2_r[i+1,c] - L2_r[i,c] + V2[i-1,0] - V2[i,0]; # dM above feed

var Mx2dot{1..NT,1..2,fe,cp};
var Mx2dot_r{1..NT,1..2,cp};
const34 {i in 2..NF-1, j in 1..2,f in fe, c in cp}: # dMxdt below feed
    Mx2dot[i,j,f,c]=L2[i+1,f,c]*x2[i+1,j,f,c] - L2[i,f,c]*x2[i,j,f,c] + V2[i-1,f]*y2[i-1,j,f,c] - V2[i,f]*y2[i,j,f,c];
    
    const34_r {i in 2..NF-1, j in 1..2, c in cp}: # dMxdt below feed
    Mx2dot_r[i,j,c]=L2_r[i+1,c]*x2_r[i+1,j,c] - L2_r[i,c]*x2_r[i,j,c] + V2[i-1,0]*y2_r[i-1,j,c] - V2[i,0]*y2_r[i,j,c];
    

const35 {i in NF+1..NT-1, j in 1..2,f in fe, c in cp}: # dMxdt above feed
    Mx2dot[i,j,f,c]=L2[i+1,f,c]*x2[i+1,j,f,c] - L2[i,f,c]*x2[i,j,f,c] + V2[i-1,f]*y2[i-1,j,f,c] - V2[i,f]*y2[i,j,f,c];
    
    const35_r {i in NF+1..NT-1, j in 1..2, c in cp}: # dMxdt above feed
    Mx2dot_r[i,j,c]=L2_r[i+1,c]*x2_r[i+1,j,c] - L2_r[i,c]*x2_r[i,j,c] + V2[i-1,0]*y2_r[i-1,j,c] - V2[i,0]*y2_r[i,j,c];

### Correction for feed at the feed stage: The feed is assumed to be
# mixed into the feed stage

const36{f in fe, c in cp}:
    M2dot[NF,f,c]=L2[NF+1,f,c] - L2[NF,f,c] + V2[NF-1,f] - V2[NF,f] + B1[f];
    
    const36_r{c in cp}:
    M2dot_r[NF,c]=L2_r[NF+1,c] - L2_r[NF,c] + V2[NF-1,0] - V2[NF,0] + B1[0];

const37{f in fe,j in 1..2, c in cp}:
    Mx2dot[NF,j,f,c]=L2[NF+1,f,c]*x2[NF+1,j,f,c] - L2[NF,f,c]*x2[NF,j,f,c] + V2[NF-1,f]*y2[NF-1,j,f,c] - V2[NF,f]*y2[NF,j,f,c] + B1[f]*x1[1,j,f,c];
    
    const37_r{j in 1..2, c in cp}:
    Mx2dot_r[NF,j,c]=L2_r[NF+1,c]*x2_r[NF+1,j,c] - L2_r[NF,c]*x2_r[NF,j,c] + V2[NF-1,0]*y2_r[NF-1,j,c] - V2[NF,0]*y2_r[NF,j,c] + B1[0]*x1_r[1,j,c];

### Reboiler (assumed to be an equilibrium stage)

const39{f in fe, c in cp}:
M2dot[1,f,c]=L2[2,f,c] - V2[1,f] - B2[f];

const39_r{c in cp}:
M2dot_r[1,c]=L2_r[2,c] - V2[1,0] - B2[0];

const40 {i in 1..2,f in fe, c in cp}:
Mx2dot[1,i,f,c]=L2[2,f,c]*x2[2,i,f,c] - V2[1,f]*y2[1,i,f,c] - B2[f]*x2[1,i,f,c];

const40_r {i in 1..2, c in cp}:
Mx2dot_r[1,i,c]=L2_r[2,c]*x2_r[2,i,c] - V2[1,0]*y2_r[1,i,c] - B2[0]*x2_r[1,i,c];

### Total condenser (no equilibrium stage)

const41{f in fe, c in cp}:
    M2dot[NT,f,c]=V2[NT-1,f] - LT2[f] - D2[f];
    
    const41_r{ c in cp}:
    M2dot_r[NT,c]=V2[NT-1,0] - LT2[0] - D2[0];

const42 {i in 1..2,f in fe, c in cp}:
    Mx2dot[NT,i,f,c]=V2[NT-1,f]*y2[NT-1,i,f,c] - L2[NT,f,c]*x2[NT,i,f,c] - D2[f]*x2[NT,i,f,c];
    
    const42_r {i in 1..2, c in cp}:
    Mx2dot_r[NT,i,c]=V2[NT-1,0]*y2_r[NT-1,i,c] - L2_r[NT,c]*x2_r[NT,i,c] - D2[0]*x2_r[NT,i,c];

var x2dot{1..NT,1..2,fe,cp};
var x2dot_r{1..NT,1..2,cp};
const42_2{i in 1..NT, j in 1..2,f in fe, c in cp}:
x2dot[i,j,f,c]*M2[i,f,c]=(Mx2dot[i,j,f,c]-x2[i,j,f,c]*M2dot[i,f,c]);

const42_2_r{i in 1..NT, j in 1..2, c in cp}:
x2dot_r[i,j,c]*M2_r[i,c]=(Mx2dot_r[i,j,c]-x2_r[i,j,c]*M2dot_r[i,c]);

#var x2dot{i in 1..NT, j in 1..2,f in fe, c in cp}=(Mx2dot[i,j,f,c]-x2[i,j,f,c]*M2dot[i,f,c])/M2[i,f,c];

param x2set{1..NT,1..2};
param M2set{1..NT};

x2_init_constr{i in 1..NT, j in 1..2}:x2_0[i,j,0]=x2set[i,j];
x2_constr{i in 1..NT, j in 1..2,f in 1..nfe}:x2_0[i,j,f]=x2[i,j,f-1,ncp];
x2_Lag{i in 1..NT,j in 1..2,q in fe, c in cp}: x2[i,j,q,c] = x2_0[i,j,q] + h* sum{k in cp} omega[k,c]*x2dot[i,j,q,k];
x2_r_Lag{i in 1..NT,j in 1..2, c in cp}: x2_r[i,j,c] = x2_0[i,j,0] + h* sum{k in cp} omega[k,c]*x2dot_r[i,j,k];


M2_init_constr{i in 1..NT}:M2_0[i,0]=M2set[i];
M2_constr{i in 1..NT,f in 1..nfe}:M2_0[i,f]=M2[i,f-1,ncp];
L2_const{i in 2..NT, f in 1..nfe}:L2_0[i,f] = L2[i,f-1,ncp];
L2_const2 {i in 2..NF}:L2_0[i,0] = Kbf*( (M2_0[i,0] - Muw+sqrt((M2_0[i,0] - Muw)^2+10^(-8)))/2    )^1.5 ; # Liquid flow below feed
L2_const3 {i in NF+1..NT-1}: L2_0[i,0] = Kuf*( (M2_0[i,0] - Muw+sqrt((M2_0[i,0] - Muw)^2+10^(-8)))/2 )^1.5 ;# Liquid flows above feed
M2_Lag{i in 1..NT,q in fe, c in cp}: M2[i,q,c] = M2_0[i,q] + h* sum{k in cp} omega[k,c]*M2dot[i,q,k];
M2_r_Lag{i in 1..NT, c in cp}: M2_r[i,c] = M2_0[i,0] + h* sum{k in cp} omega[k,c]*M2dot_r[i,k];

#const43 and const44 are not in original model
const43 {i in 1..NT,f in 0..nfe}:
     TC1[i,f] = x1_0[i,1,f]*353.3 + x1_0[i,2,f]*383.8 + (1-x1_0[i,1,f]-x1_0[i,2,f])*411.5 ;

const44 {i in 1..NT, f in 0..nfe}:
     TC2[i,f] = x2_0[i,1,f]*353.3 + x2_0[i,2,f]*383.8 + (1-x2_0[i,1,f]-x2_0[i,2,f])*411.5 ; #x1[14,1]*353.3+x1[14,2]*383.8+(1-x1[14,1]-x1[14,2])*411.5 = 382.9450;

#const48{f in fe}:
 #   VB1[f] <= 4.008 ; # max boilup column 1

#const49{f in fe}:
#    VB2[f] <= 2.405 ; # max boilup column 2

### Initial constraints for parameters


#--------------------
# Terminal Constraint
#--------------------

#termcon_x1 {i in 1..NT,j in 1..2} : x1_0[i,j,nfe]=x1ref[i,j] ; 
#termcon_x2 {i in 1..NT,j in 1..2} : x2_0[i,j,nfe]=x2ref[i,j] ; 
#termcon_m1 {i in 1..NT} : M1_0[i,nfe]=M1ref[i] ; 
#termcon_m2 {i in 1..NT} : M2_0[i,nfe]=M2ref[i] ; 
termcon : sum {i in 1..NT}(sum{j in 1..2}((x1_0[i,j,nfe]-x1ref[i,j])^2 +(x2_0[i,j,nfe]-x2ref[i,j])^2 )) + sum {i in 1..NT}( (M1_0[i,nfe]-M1ref[i])^2 + (M2_0[i,nfe]-M2ref[i])^2) <=    termnorm^2 +termeps ; 

#--------------------
# Descent Constraint
#--------------------

# 2 norm
#descent_con1: lk =(D1[0]-D1ref)^2+(D2[0]-D2ref)^2+(B1[0]-B1ref)^2+(B2[0]-B1ref)^2+(VB1[0]-VB1ref)^2+(VB2[0]-VB2ref)^2+(LT1[0]-LT1ref)^2+(LT2[0]-LT2ref)^2 + sum{i in 1..NT}( (M1_0[i,0]-M1ref[i])^2 +(M2_0[i,0]-M2ref[i])^2  + sum{j in 1..2}(   (x1_0[i,j,0]-x1ref[i,j])^2 +(x2_0[i,j,0]-x2ref[i,j])^2)    )  ;

#descent_con2: vk= sum{f in fe}((D1[f]-D1ref)^2+(D2[f]-D2ref)^2+(B1[f]-B1ref)^2+(B2[f]-B2ref)^2+(VB1[f]-VB1ref)^2+(VB2[f]-VB2ref)^2+(LT1[f]-LT1ref)^2+(LT2[f]-LT2ref)^2+  sum{i in 1..NT}( (M1_0[i,f]-M1ref[i])^2 +(M2_0[i,f]-M2ref[i])^2  + sum{j in 1..2}(   (x1_0[i,j,f]-x1ref[i,j])^2 +(x2_0[i,j,f]-x2ref[i,j])^2)    ) )   ; 
 

# 1 norm

#descent_con1: lk = sqrt((D1[0]-D1ref)^2+smt^2)-smt+sqrt((D2[0]-D2ref)^2+smt^2)-smt+sqrt((B1[0]-B1ref)^2+smt^2)-smt+sqrt((B2[0]-B1ref)^2+smt^2)-smt+sqrt((VB1[0]-VB1ref)^2+smt^2)-smt+sqrt((VB2[0]-VB2ref)^2+smt^2)-smt+sqrt((LT1[0]-LT1ref)^2+smt^2)-smt+sqrt((LT2[0]-LT2ref)^2+smt^2)-smt + sum{i in 1..NT}( sqrt((M1_0[i,0]-M1ref[i])^2+smt^2)-smt +sqrt((M2_0[i,0]-M2ref[i])^2+smt^2)-smt  + sum{j in 1..2}(   sqrt((x1_0[i,j,0]-x1ref[i,j])^2+smt^2)-smt +sqrt((x2_0[i,j,0]-x2ref[i,j])^2+smt^2)-smt   )    )  ;

#descent_con2: vk= sum{f in fe}(sqrt((D1[f]-D1ref)^2+smt^2)-smt+sqrt((D2[f]-D2ref)^2+smt^2)-smt+sqrt((B1[f]-B1ref)^2+smt^2)-smt+sqrt((B2[f]-B1ref)^2+smt^2)-smt+sqrt((VB1[f]-VB1ref)^2+smt^2)-smt+sqrt((VB2[f]-VB2ref)^2+smt^2)-smt+sqrt((LT1[f]-LT1ref)^2+smt^2)-smt+sqrt((LT2[f]-LT2ref)^2+smt^2)-smt + sum{i in 1..NT}( sqrt((M1_0[i,f]-M1ref[i])^2+smt^2)-smt +sqrt((M2_0[i,f]-M2ref[i])^2+smt^2)-smt  + sum{j in 1..2}(   sqrt((x1_0[i,j,f]-x1ref[i,j])^2+smt^2)-smt +sqrt((x2_0[i,j,f]-x2ref[i,j])^2+smt^2)-smt   )    )   )   ; 


# descent constraint 

#descent_con3: vk <= vkmin1 - sigma*lkmin1 + desceps ; 

#--------------------
# OBJECTIVE FUNCTION
#--------------------

tcon1{j in 1..NT}:  stat_vec[j] = x1_0[j,1,nfe]-x1ref[j,1];
tcon2{j in 1..NT}: stat_vec[NT+j] = x1_0[j,2,nfe]-x1ref[j,2];
tcon3{j in 1..NT}:  stat_vec[2*NT+j]  = x2_0[j,1,nfe]-x2ref[j,1];
tcon4{j in 1..NT}: stat_vec[3*NT+j]  = x2_0[j,2,nfe]-x2ref[j,2];
tcon5{j in 1..NT}: stat_vec[4*NT+j]  = M1_0[j,nfe]-M1ref[j] ; 
tcon6{j in 1..NT}: stat_vec[5*NT+j]  = M2_0[j,nfe]-M2ref[j] ; 


#var cost{f in fe}=pF*F + pV*(VB1[f] + VB2[f])-pA*D1[f] - pB*D2[f] - pC*B2[f];

var termcost = sum{i in s, j in s} stat_vec[i]*P[i,j]*stat_vec[j] ; 

var penalty_2{f in fe}=D1eps[f] +B1eps[f]+D2eps[f] +B2eps[f];
var penalty_M_TC{f in 0..nfe}=sum{i in 1..NT}(M1eps[i,f] +M2eps[i,f]+TC1eps[i,f]+TC2eps[i,f]);
var penalty_x_y{f in 0..nfe}=x2ceps[f]+sum{i in 1..NT}sum{j in 1..2}(x1eps[i,j,f]+x2eps[i,j,f]);
var penalty_V_L{f in fe}=sum{i in 1..NT-1}(V1eps[i,f]+V2eps[i,f])+sum{i in 2..NT}(L1eps[i,f]+L2eps[i,f]);

var tracking {f in fe} = cont_w*((VB1[f]-VB1ref)^2 + (VB2[f]-VB2ref)^2+ (LT1[f]-LT1ref)^2 + (LT2[f]-LT2ref)^2+(D1[f]-D1ref)^2 +(B1[f]-B1ref)^2+(D2[f]-D2ref)^2 + (B2[f]-B2ref)^2) + state_w*(sum {i in 1..NT}((M1_0[i,f]-M1ref[i])^2 +(M2_0[i,f]-M2ref[i])^2+sum{j in 1..2}( (x1_0[i,j,f]-x1ref[i,j])^2 +(x2_0[i,j,f]-x2ref[i,j])^2  ) )  );

#var regu_1{f in fe}=F_w*(F-Fref)^2 +qF_w*(qF-qFref)^2+pV_w*(pV-pVref)^2+VB1_w*(VB1[f]-VB1ref)^2 +VB2_w*(VB2[f]-VB2ref)^2+LT1_w*(LT1[f]-LT1ref)^2 +LT2_w*(LT2[f]-LT2ref)^2;
#var regu_2{f in fe}=(D1_w*(D1[f]-D1ref)^2 +B1_w*(B1[f]-B1ref)^2+D2_w*(D2[f]-D2ref)^2 +B2_w*(B2[f]-B2ref)^2);
#var regu_M_TC{f in fe}=sum{i in 1..NT}(M1_w[i]*(M1_0[i,f]-M1ref[i])^2 +M2_w[i]*(M2_0[i,f]-M2ref[i])^2+TC1_w[i]*(TC1[i,f]-TC1ref[i])^2 +TC2_w[i]*(TC2[i,f]-TC2ref[i])^2);
#var regu_x_y{f in fe}=sum{i in 1..NT}(sum{j in 1..2}((x1_w[i,j]*(x1_0[i,j,f]-x1ref[i,j])^2 +x2_w[i,j]*(x2_0[i,j,f]-x2ref[i,j])^2)))+sum{i in 1..NT-1}(sum{j in 1..2}((y1_w[i,j]*(y1_0[i,j,f]-y1ref[i,j])^2 +y2_w[i,j]*(y2_0[i,j,f]-y2ref[i,j])^2)));
#var regu_V_L{f in fe}=sum{i in 1..NT-1}(V1_w[i]*(V1[i,f]-V1ref[i])^2 +V2_w[i]*(V2[i,f]-V2ref[i])^2)+sum{i in 2..NT}(L1_w[i]*(L1_0[i,f]-L1ref[i])^2 +L2_w[i]*(L2_0[i,f]-L2ref[i])^2);
#var regu_y{f in fe}=sum{i in 1..NT-1}(sum{j in 1..2}(y_1_1_w[i,j]*(y_1_1_0[i,j,f]-y_1_1ref[i,j])^2 +y_1_2_w[i,j]*(y_1_2_0[i,j,f]-y_1_2ref[i,j])^2+y_2_1_w[i,j]*(y_2_1_0[i,j,f]-y_2_1ref[i,j])^2 +y_2_2_w[i,j]*(y_2_2_0[i,j,f]-y_2_2ref[i,j])^2));
#var regu{f in fe}=regu_1[f]+regu_2[f]+regu_M_TC[f]+regu_x_y[f]+regu_V_L[f]+regu_y[f];

var penalty1{f in fe}=(penalty_2[f]+penalty_M_TC[f]+penalty_x_y[f]+penalty_V_L[f]);

var penalty2{f in 0..nfe}=(penalty_M_TC[f]+penalty_x_y[f]);

minimize obj: rho_psi*termcost+rho*pen_switch*termeps+sum{f in 0..nfe}rho*pen_switch*penalty2[f]+sum{f in fe}(rho*pen_switch*penalty1[f]+track_switch*tracking[f]/rho_L);
