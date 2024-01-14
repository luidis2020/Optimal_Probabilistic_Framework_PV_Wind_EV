set N;
set L in N cross N;
set F := 1..3;
#set T
set H; # conjunto de horas

param factor_gen {H};
param factor_load;

param SE;
param vb; 
param pi := 3.141592653589793;
param r3 := 1.732050807568877;



param Rp      {L,1..18};		# Tabela de valores de R e X
param cp      {L,1..9};			# Tabela de valores C
param B       {L,F};
param hours {H};
param fatorPV {H};
param OM_PV; #O&M cost (US$/kWh)


param Vmin    {N};
param Vmax    {N};
param Pd      {N,F};
param Qd      {N,F};

param R       {L,F,F};
param X       {L,F,F};
param Bshunt  {L,F,F};

param Imax    {L};
param at      {L};				


set G;
param pdgmax  {G};
var Pg        {G} >= 0;
var Qg        {G};

var Vre       {N,F};
var Vim       {N,F};
var Idre      {N,F}; # Real part of load current
var Idim      {N,F};
var Igre      {N,F};
var Igim      {N,F};
var Isubre      {N,F};
var Isubim      {N,F};
var Ire       {L,F}; # Real part of line current
var Iim       {L,F};
param fem_pv;
param fem_se;
param COST_emission;
param CENER;#cost of substation energy


#-sum{(m,n) in L} (Bshunt[m,n,f,f]*Vim[m,f])
#subject to Ire_linhas {n in N, f in F}:
#  - sum{(n,m) in L} Ire[n,m,f] + sum{(m,n) in L} Ire[m,n,f] = Idre[n,f] - Igre[n,f];
  


# For some reason, the + simbol (Bshunt) gives nice resutls *** read theory documention
#

subject to Ire_linhas {n in N, f in F}:
  - sum{(n,m) in L} Ire[n,m,f] + sum{(m,n) in L} Ire[m,n,f] 
   +sum{h in F}((sum{(n,m) in L} Bshunt[n,m,f,h] + sum{(m,n) in L} Bshunt[m,n,f,h])* Vim[n,h]/2)
   =  Idre[n,f] - Igre[n,f]-Isubre[n,f];  
   
#subject to Ire_linhas {n in N, f in F}:
#  - sum{(n,m) in L} Ire[n,m,f] + sum{(m,n) in L} Ire[m,n,f] 
#   -sum{h in F}((sum{(n,m) in L} Bshunt[n,m,f,h])* Vim[n,h])
#   =  Idre[n,f] - Igre[n,f];     
   

#+sum{(m,n) in L}(Bshunt[m,n,f,f]*Vre[m,f])
#subject to Iim_linhas {n in N, f in F}:
#  - sum{(n,m) in L} Iim[n,m,f] + sum{(m,n) in L} Iim[m,n,f] = Idim[n,f] - Igim[n,f];
  
# For some reason, the - simbol (Bshunt) gives nice resutls *** read theory documention
#  
subject to Iim_linhas {n in N, f in F}:
  - sum{(n,m) in L} Iim[n,m,f] + sum{(m,n) in L} Iim[m,n,f] 
  -sum{h in F}((sum{(n,m) in L} Bshunt[n,m,f,h] + sum{(m,n) in L} Bshunt[m,n,f,h])* Vre[n,h]/2)
  =  Idim[n,f] - Igim[n,f]-Isubim[n,f];  
  
# subject to Iim_linhas {n in N, f in F}:
#  - sum{(n,m) in L} Iim[n,m,f] + sum{(m,n) in L} Iim[m,n,f] 
#  +sum{h in F}((sum{(n,m) in L} Bshunt[n,m,f,h])* Vre[n,h])
#  =  Idim[n,f] - Igim[n,f];  
   
  
subject to Vre_linha {(m,n) in L, f in F}:
  Vre[m,f] - Vre[n,f] = sum{r in F} ( R[m,n,f,r] * Ire[m,n,r] - X[m,n,f,r] * Iim[m,n,r] );
  
#subject to Vre_linha {(m,n) in L, f in F}:
#  Vre[m,f] - Vre[n,f] = sum{r in F} ( R[m,n,f,r] * (Ire[m,n,r]-(Bshunt[m,n,f,r])* Vim[n,r]) - #X[m,n,f,r] * (Iim[m,n,r]+ (Bshunt[m,n,f,r])* Vre[n,r]));

  

subject to Vim_linha {(m,n) in L, f in F}:
  Vim[m,f] - Vim[n,f] = sum{r in F} ( X[m,n,f,r] * Ire[m,n,r] + R[m,n,f,r] * Iim[m,n,r] );

#subject to Vim_linha {(m,n) in L, f in F}:
#  Vim[m,f] - Vim[n,f] = sum{r in F} ( X[m,n,f,r] * (Ire[m,n,r]-(Bshunt[m,n,f,r])* Vim[n,r]) + #R[m,n,f,r] * (Iim[m,n,r]+ (Bshunt[m,n,f,r])* Vre[n,r]) );



subject to Ire_CARGAS_nolinear {n in N, f in F}:
  Idre[n,f] = factor_load*( Pd[n,f] * Vre[n,f] + Qd[n,f] * Vim[n,f] ) / ( Vre[n,f]^2 + Vim[n,f]^2 );

subject to Iim_CARGAS_nolinear {n in N, f in F}:
  Idim[n,f] = factor_load*( -Qd[n,f] * Vre[n,f] + Pd[n,f] * Vim[n,f] ) / ( Vre[n,f]^2 + Vim[n,f]^2 );


subject to Pg_gerador_limit {n in G,h in H}:
  Pg[n] <= pdgmax[n]*factor_gen[h];
# Pg[n] = 0;
subject to Qg_gerador_limit {n in G}:
# Qg[n] <= Pg[n] * tan ( acos (0.80) );
Qg[n] = 0;

#subject to Pg_gerador_limitpenetration :
#  sum {n in G} Pg[n] <= sum {n in N,f in F} 0.5*Pd[n,f];

subject to Ire_gerador {n in G, f in F}:
  Igre[n,f] = (  Pg[n] * Vre[n,f] + Qg[n] * Vim[n,f] ) / ( Vre[n,f]^2 + Vim[n,f]^2 ) / 3;

subject to Iim_gerador {n in G, f in F}:
  Igim[n,f] = ( -Qg[n] * Vre[n,f] + Pg[n] * Vim[n,f] ) / ( Vre[n,f]^2 + Vim[n,f]^2 ) / 3;
  
  
subject to Vlimit_min {n in N, f in F}:
(Vmin[n]*vb)^2 <= Vre[n,f]^2+Vim[n,f]^2;

subject to Vlimit_max {n in N, f in F}:
(Vmax[n]*vb)^2 >= Vre[n,f]^2+Vim[n,f]^2;


#Power reversion
#subject to Power_reversion {f in F}:
#(Vre[SE,f]*Isubre[SE,f]+ Vim[SE,f]*Isubim[SE,f])>=0;

var OMPV = sum{n in G,h in H} Pg[n]*OM_PV*hours[h];
var PSE=  sum{f in F}  (Vre[SE,f]*Isubre[SE,f]+ Vim[SE,f]*Isubim[SE,f]);
var energySE=  sum{f in F,h in H}  (Vre[SE,f]*Isubre[SE,f]+ Vim[SE,f]*Isubim[SE,f])*hours[h];
var cOSTPSE=  (abs (energySE))*CENER;
var EMdg_pv = sum {n in G,h in H}(hours[h]*fem_pv*Pg[n]/1000)*COST_emission;
var EM_se = sum {h in H}((hours[h]*fem_se*(abs(PSE)))/1000)*COST_emission;



#minimize COST : OMPV + EMdg_pv + cOSTPSE + EM_se;
#minimize  COST : OMPV+EMdg_pv
maximize COST: sum{n in G,h in H} Pg[n];









  
  
