var pi          ${\pi}$                 (long_name='inflation')
    y_gap       ${\tilde y}$            (long_name='output gap')
    y_nat       ${y^{nat}}$             (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y           ${y}$                   (long_name='output')
    r_nat       ${r^{nat}}$             (long_name='natural interest rate')   
    i           ${i}$                   (long_name='nominal interrst rate')
    n           ${n}$                   (long_name='hours worked')
    nu          ${\nu}$                 (long_name='AR(1) monetary policy shock process')    
    a           ${a}$                   (long_name='AR(1) technology shock process')
    p           ${p}$                   (long_name='price level')
    w           ${w}$                   (long_name='nominal wage')
    c           ${c}$                   (long_name='consumption')
    mu          ${\mu}$                 (long_name='markup')
;     

varexo  eps_a       ${\varepsilon_a}$       (long_name='technology shock')
        eps_nu  ${\varepsilon_\nu}$     (long_name='monetary policy shock')
       ;

parameters alppha       ${\alpha}$     (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
    kappa
    psi_n_ya
    ;
%----------------------------------------------------------------
% Parametrization, p. 67-75
%----------------------------------------------------------------
siggma = 1;        
varphi=5;          
phi_pi = 1.5;       
phi_y  = 0.125;     
theta=3/4;          
rho_nu =0.5;    	  
rho_a  = 0.9;       
betta  = 0.99;       
alppha=1/4;     	
epsilon=9;      
//Composite parameters
Omega=(1-alppha)/(1-alppha+alppha*epsilon);        %defined on page 60
psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha);   %defined on page 62
lambda=(1-theta)*(1-betta*theta)/theta*Omega;      %defined on page 61
kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));     %defined on page 63    

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
[name='New Keynesian Phillips Curve eq. (22)']
pi=betta*pi(+1)+kappa*y_gap;

[name='Dynamic IS Curve eq. (23)']
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1);

[name='Interest Rate Rule eq. (26)']
i=phi_pi*pi+phi_y*y_gap+nu;

[name='Definition natural rate of interest eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a;

[name='Definition natural output, eq. (20)']
y_nat=psi_n_ya*a;

[name='Definition output gap']
y_gap=y-y_nat;

[name='Monetary policy shock']
nu=rho_nu*nu(-1)+eps_nu;

[name='TFP shock']
a=rho_a*a(-1)+eps_a;

[name='Production function (eq. 14)']
y=a+(1-alppha)*n;

[name='Definition price level']
pi=p-p(-1);

[name='resource constraint, eq. (12)']
y=c;

[name='FOC labor, eq. (2)']
w-p=siggma*c+varphi*n;

[name='average price markup, eq. (18)']
mu=-(siggma+(varphi+alppha)/(1-alppha))*y+(1+varphi)/(1-alppha)*a + log(1-alppha);
end;


%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------


shocks;
        var eps_nu = 0.25^2; 
	  var eps_a  = 1^2;
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid;
steady;
check;

stoch_simul(order = 1,irf=15,irf_plot_threshold=0) y_gap pi y n w p i r_nat;



