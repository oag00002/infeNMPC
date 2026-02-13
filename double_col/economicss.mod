##########################
# The distillation model #
##########################

### Two columns in series
### Three components: A (light), B and C (heavy).
### Based on Column A. Steady state model.
### made by Roald Brck Leer

#comment by Xue: This is a steady state problem

#--------------------------------------------
# Setting parameters as variables
#--------------------------------------------

# Disturbances
var F >= 0.1 ; #feed
var qF ; #liquid fraction in feed
var pV ; #energy price

#-----------------
# Parameters
#-----------------

param NT := 41 ; #stages
param NC := 3  ; #number of components
param NF := 21 ; #feed enter at stage 21

### Relative volatility, alpha(A/B) alpha(B/C)
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
#let zF[1,1] := 0.4 ;#
#let zF[1,2] := 0.2 ;#
#let zF[1,3] := 0.4 ;

#------------
# Variables
#------------

### FIRST COLUMN ----------------------------------------------------------
    var LT1 >= 0 <= 10 := 3.43656 ;
    var VB1 >= 0 <= 10 := 4.008 ;

### comp frac liquid and vap
    var y1 {i in 1..NT-1, j in 1..2} >= 0 := 0.3 ; # vapor comp
    var x1 {i in 1..NT, j in 1..2} >= 0 :=0.4    ; # liquid comp
    var M1 {i in 1..NT} >= 0 := 0.5              ; # holdup
    var V1 {1..NT-1} >= 0 := 1                   ; # vapor flow
    var L1 {2..NT} >= 0 := 1                     ; # liquid flow
    var D1 >= 0 := 0.57                          ; # distillate flow
    var B1 >= 0 := 0.83                          ; # bottoms flow

# VLE equation split
    var y_1_1 {i in 1..NT-1, j in 1..2} >= 0 := 0.2 ;
    var y_1_2 {i in 1..NT-1, j in 1..2} >= 0 := 0.7 ;
   #TC1, SET1 and SET2 are not the original model
    var TC1 {i in 1..NT} >= 350 <= 412 := 380;
    var SET1 :=-150;
    var SET2 :=-100;
    #var SET3 :=-270;

### SECOND COLUMN----------------------------------------------------------
    var LT2 >= 0 <= 10 := 2.13827 ;
    var VB2 >= 0 <= 10 := 2.40367 ;
    
### comp frac liquid and gas
    var y2 {i in 1..NT-1, j in 1..2} >= 0 := 0.5 ; # vapor comp
    var x2 {i in 1..NT, j in 1..2} >= 0 := 0.4 ;   # liquid comp
    var M2 {i in 1..NT} >= 0 := 0.5 ;              # holdup
    var V2 {1..NT-1} >= 0 :=1;                     # vapor flow
    var L2 {2..NT} >= 0 :=1;                       # liquid flow
    var D2 >= 0 := 0.26 ;                          # distillate flow
    var B2 >= 0 := 0.56 ;                          # bottoms flow
    
# VLE equation split
    var y_2_1 {i in 1..NT-1, j in 1..2} := 0.3 ;
    var y_2_2 {i in 1..NT-1, j in 1..2} := 0.8 ;
    #TC2 is not the in original model
    var TC2 {i in 1..NT} >= 350 <= 412 := 380;

#----------------------------
# Model/Constraints 1
#----------------------------

### VLE equation split

const1 {i in 1..NT-1, j in 1..2} :
    y_1_1[i,j] = x1[i,j]*m_alpha[j,j] ;

const2 {i in 1..NT-1, j in 1..2}:
    y_1_2[i,j] = ((x1[i,1]*(alpha[1]-1)+x1[i,2]*(alpha[2]-1))+1) ;
    
### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const3 {i in 1..NT-1, j in 1..2}:
    y1[i,j] = y_1_1[i,j]/y_1_2[i,j] ;
    
### Vapor flows assuming constant molar flows##############################

const4 {i in 1..NF-1}:
    V1[i] = VB1 ; # vapor flow below feed

const5 {i in NF..NT-1}:
    V1[i] = VB1 + (1-qF)*F ; # vapor flow above feed
    
### Liquid flows are given by Franci's Weir Formula L(i)=K*Mow(i)^1.5######
# Liquid flow L(i) dependent only on the holdup over the weir Mow(i) ######
#M(i)= Mow(i) + Muw(i) (Total holdup = holdup over weir +            ######
#holdup below weir)                                                  ######

const6 {i in 2..NF}:
    L1[i] = Kbf*(M1[i] - Muw)^1.5 ; # Liquid flow below feed

const7 {i in NF+1..NT-1}:
    L1[i] = Kuf*(M1[i] - Muw)^1.5 ;# Liquid flows above feed

const8:
    L1[NT] = LT1; # Condenser's liquid flow

### Distillate and bottom

const9:
    B1 = Bs + (M1[1] - MBs)*KcB ;

const10:
    D1 = Ds + (M1[NT] - MDs)*KcD ;

### Material balances for total holdup and component holdup################

const11 {i in 2..NF-1}:
    L1[i+1] - L1[i] + V1[i-1] - V1[i] = 0 ; # dM below feed

const12 {i in NF+1..NT-1}:
    L1[i+1] - L1[i] + V1[i-1] - V1[i] = 0 ; # dM above feed

const13 {i in 2..NF-1, j in 1..2}:
    L1[i+1]*x1[i+1,j] - L1[i]*x1[i,j] + V1[i-1]*y1[i-1,j] - V1[i]*y1[i,j] = 0 ;

const14 {i in NF+1..NT-1, j in 1..2}:
L1[i+1]*x1[i+1,j] - L1[i]*x1[i,j] + V1[i-1]*y1[i-1,j] - V1[i]*y1[i,j] = 0 ;

### Correction for feed at the feed stage: The feed is assumed to
#be mixed into the feed stage

const15:
    L1[NF+1] - L1[NF] + V1[NF-1] - V1[NF] + F = 0 ;

const16:
    L1[NF+1]*x1[NF+1,1] - L1[NF]*x1[NF,1] + V1[NF-1]*y1[NF-1,1] - V1[NF]*y1[NF,1] + F*zF[1,1] = 0 ;

const17:
    L1[NF+1]*x1[NF+1,2] - L1[NF]*x1[NF,2] + V1[NF-1]*y1[NF-1,2] - V1[NF]*y1[NF,2] + F*zF[1,2] = 0 ;

### Reboiler (assumed to be an equilibrium stage)

const18:
    L1[2] - V1[1] - B1 = 0 ;
    
const19 {j in 1..2}:
    L1[2]*x1[2,j] - V1[1]*y1[1,j] - B1*x1[1,j] = 0 ;

### Total condenser (no equilibrium stage)

const20:
    V1[NT-1] - LT1 - D1 = 0 ;

const21 {j in 1..2}:
    V1[NT-1]*y1[NT-1,j] - L1[NT]*x1[NT,j] - D1*x1[NT,j] = 0 ;

#----------------------------
# Model/Constraints 2
#----------------------------

### VLE equation split

const22 {i in 1..NT-1, j in 1..2} :
    y_2_1[i,j] = x2[i,j]*m_alpha[j,j] ;
    
const23 {i in 1..NT-1, j in 1..2}:
    y_2_2[i,j] = ((x2[i,1]*(alpha[1]-1)+x2[i,2]*(alpha[2]-1))+1) ;

### Vapour-liquid equilibria (multicomponent ideal VLE,
# Stichlmair-Fair, 'Distillation', p. 36, 1998)

const24 {i in 1..NT-1, j in 1..2}:
    y2[i,j] = y_2_1[i,j]/y_2_2[i,j] ;
    
### Vapor flows assuming constant molar flows

const25 {i in 1..NF-1}:
    V2[i] = VB2 ; # vapor flow below feed

const26 {i in NF..NT-1}:
    V2[i] = VB2 ; # vapor flow above

const27 {i in 2..NF}:
    L2[i] = Kbf*(M2[i] - Muw)^1.5 ; # Liquid flow below feed

const28 {i in NF+1..NT-1}:
    L2[i] = Kuf*(M2[i] - Muw)^1.5 ; # Liquid flow above feed

const29:
    L2[NT] = LT2; # Condenser's liquid flow

const30:
    B2 = Bs + (M2[1] - MBs)*KcB ;

const31:
    D2 = Ds + (M2[NT] - MDs)*KcD ;

### Material balances for total holdup and component holdup

const32 {i in 2..NF-1}:
    L2[i+1] - L2[i] + V2[i-1] - V2[i] = 0 ; # dM below feed

const33 {i in NF+1..NT-1}:
    L2[i+1] - L2[i] + V2[i-1] - V2[i] = 0 ; # dM above feed

const34 {i in 2..NF-1, j in 1..2}: # dMxdt below feed
    L2[i+1]*x2[i+1,j] - L2[i]*x2[i,j] + V2[i-1]*y2[i-1,j] - V2[i]*y2[i,j] = 0 ;

const35 {i in NF+1..NT-1, j in 1..2}: # dMxdt above feed
    L2[i+1]*x2[i+1,j] - L2[i]*x2[i,j] + V2[i-1]*y2[i-1,j] - V2[i]*y2[i,j] = 0 ;

### Correction for feed at the feed stage: The feed is assumed to be
# mixed into the feed stage

const36:
    L2[NF+1] - L2[NF] + V2[NF-1] - V2[NF] + B1 = 0 ;

const37:
    L2[NF+1]*x2[NF+1,1] - L2[NF]*x2[NF,1] + V2[NF-1]*y2[NF-1,1] - V2[NF]*y2[NF,1] + B1*x1[1,1] = 0 ;

const38:
    L2[NF+1]*x2[NF+1,2] - L2[NF]*x2[NF,2] + V2[NF-1]*y2[NF-1,2] - V2[NF]*y2[NF,2] + B1*x1[1,2] = 0 ;

### Reboiler (assumed to be an equilibrium stage)

const39:
L2[2] - V2[1] - B2 = 0 ;

const40 {i in 1..2}:
L2[2]*x2[2,i] - V2[1]*y2[1,i] - B2*x2[1,i] = 0 ;

### Total condenser (no equilibrium stage)

const41:
    V2[NT-1] - LT2 - D2 = 0 ;

const42 {i in 1..2}:
    V2[NT-1]*y2[NT-1,i] - L2[NT]*x2[NT,i] - D2*x2[NT,i] = 0 ;

#const43 and const44 are not in original model
const43 {i in 1..NT}:
     TC1[i] = x1[i,1]*353.3 + x1[i,2]*383.8 + (1-x1[i,1]-x1[i,2])*411.5 ;

const44 {i in 1..NT}:
     TC2[i] = x2[i,1]*353.3 + x2[i,2]*383.8 + (1-x2[i,1]-x2[i,2])*411.5 ; #x1[14,1]*353.3+x1[14,2]*383.8+(1-x1[14,1]-x1[14,2])*411.5 = 382.9450;

### Inequality constraints

const45:
    x1[NT,1] >= 0.95 ; # xA

const46:
    x2[NT,2] >= 0.95 ; # xB (= 0.95)

const47:
    1 - (x2[1,1]+x2[1,2]) >= 0.95 ; #

#const48:
#    VB1 <= 4.008 ; # max boilup column 1

#const49:
#    VB2 <= 2.405 ; # max boilup column 2

### Initial constraints for parameters
param nominal_F;
param nominal_pV;
param nominal_qF;

constF: F = nominal_F ;
constpV: pV = nominal_pV ;
# constzF: zF = nominal_zF ;
constqF: qF = nominal_qF ;

#--------------------
# OBJECTIVE FUNCTION
#--------------------

minimize cost: pF*F + pV*(VB1 + VB2) - pA*D1 - pB*D2 - pC*B2 ;