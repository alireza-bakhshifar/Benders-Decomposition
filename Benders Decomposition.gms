********************************************************
************** BD  *************************************
********************************************************

Set
i /i1*i50/
j /j1*j100/
m /m1*m40/
n /n1*n20/
;

*********** Data
Parameters
c(i)
f(j)

A(m,i)
D(m,j)
b(m)

K(n,j)
e(n)
;


c(i) = uniform(3,10);
f(j) = uniform(40,60);

A(m,i)= uniform(2,5);
D(m,j)= uniform(50,100);
b(m)  = uniform(1000,10000);

K(n,j)= uniform(1,5);
e(n)  = uniform(50,100);

Display
c
f
A
D
b
K
e
;


********************************************************************************
* Step 0: Initialization of BD**************************************************

Set iter /1*20/

Parameters
MaxRE   /1E-3/
LB
UB
;
LB=0;
UB=inf;
********************************************************************************

********************************************************************************
* Step 1: Subproblem (SP) of BD*************************************************

Free variables
Z_SP    ;

Positive variable
x(i) ;

Parameter
y_fix(j) ;

Equations
obj_SP       'obective function of SP'
Cons1_SP      'Constraints of SP'
;

obj_SP..          Z_SP =e= sum(i,c(i)*x(i)) + sum(j,f(j)*y_fix(j));

Cons1_SP(m)..     sum(i,A(m,i)*x(i)) =g= b(m) - sum(j,D(m,j)*y_fix(j)) ;

Model SP
/
obj_SP
Cons1_SP
/
;
********************************************************************************
********************************************************************************


********************************************************************************
* Step 2: Dual Subproblem (DSP) of BD*******************************************

Free variables
Z_DSP    ;

Positive variable
u(m) ;


Equations
obj_DSP       'obective function of DSP'
Cons_DSP      'Constraints of DSP'
;

obj_DSP..         Z_DSP =e=  sum{m,( b(m) - sum(j,D(m,j)*y_fix(j)))*u(m)} + sum(j,f(j)*y_fix(j));

Cons_DSP(i)..     sum(m,A(m,i)*u(m)) =l= c(i) ;

Model DSP
/
obj_DSP
Cons_DSP
/
;
********************************************************************************
********************************************************************************


********************************************************************************
* Step 3: Homogeneous/Modified Dual Subproblem (MDSP) of BD*********************

Free variables
Z_MDSP    ;

Equations
obj_MDSP       'obective function of MDSP'
Cons_MDSP      'Constraints of MDSP'
Bounded_Cons
;

obj_MDSP..        Z_MDSP =e= sum{m,( b(m) - sum(j,D(m,j)*y_fix(j)))*u(m)};
Cons_MDSP(i)..    sum(m,A(m,i)*u(m)) =l= 0 ;
Bounded_Cons..    sum(m,u(m)) =l= card(m);

Model MDSP
/
obj_MDSP
Cons_MDSP
Bounded_Cons
/
;
********************************************************************************
********************************************************************************



********************************************************************************
* Step 4: Relaxed Master Problem (RMP) of BD************************************

Free variables
Z_RMP    ;

Binary variable
y(j)
;

Positive Variable Say;

Sets
OC(iter)
FC(iter)
;

OC(iter)= NO;
FC(iter)= NO;

Parameter
u_fix(m,iter)
;



Equations
obj_RMP       'obective function of MP'
Cons2_RMP

OptimalityCut(iter)
FeasibilityCut(iter)
;

obj_RMP..                Z_RMP =e= sum(j,f(j)*y(j)) + Say;

Cons2_RMP(n)..               sum(j,K(n,j)*y(j)) =l= e(n) ;

OptimalityCut(OC)..      Say  =g=  sum{m,(b(m) - sum(j,D(m,j)*y(j)))*u_fix(m,oc)} ;

FeasibilityCut(FC)..     sum{m,(b(m) - sum(j,D(m,j)*y(j)))*u_fix(m,fc)} =l= 0;

Model RMP
/
obj_RMP
Cons2_RMP
OptimalityCut
FeasibilityCut
/
;
********************************************************************************
********************************************************************************



********************************************************************************
* Step 5: Main Loop for implentation of BD**************************************


Parameters
Result(iter,*)
Converged /NO/
Iteration
Gap
Y_Feasibility
;

y_fix(j)=0;

Options
LP  =  CPLEX
MIP = CPLEX
OPTCR = 0
RESLIM = 300
;

Loop(iter$(NOT(Converged)),

***** Solve DSP or MDSP to find u and update UB
Solve DSP us LP MAX Z_DSP;

*Infeasibility of Main OP
Abort$(DSP.ModelStat = 2) "Your OP Model is not fasible"

* Bounded Situtation
if( DSP.ModelStat <> 3,

Y_Feasibility = YES;
OC(iter)=YES;
UB= Z_DSP.l;
Result(iter,'UB')=UB;

ELSE
* Unbounded Situtation
Y_Feasibility = NO;
FC(iter)=YES;
Solve MDSP us LP MAX Z_MDSP;

)
* end of If
;

Result(iter,'Feasible')=Y_Feasibility;
u_fix(m,iter)=u.l(m);
************************************************


*************** Solve RMP to  fine new y and y and update LB

Solve RMP us MIP MIN Z_RMP ;

Abort$(RMP.ModelStat = 2) "Your OP Model is not fasible" ;

LB=Z_RMP.l;
Result(iter,'LB')=LB;

y_fix(j)=y.l(j);
******************************************


* Stop Criteria

Gap = abs((UB - LB)/ UB );

Result(iter,'Gap')=Gap;

IF(Gap <= MaxRE,
Converged = YES;
)
;

Iteration=ord(iter);
Display
"-----------Iteration-----------"
Iteration
Y_Feasibility
LB
UB
Gap
"-----------Varable-------------"
"Y"
y_fix
"X"
Cons_DSP.m

)
;
*End of BD Main Loop


Display
Result


