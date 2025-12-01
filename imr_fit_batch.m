%                   Simulation of Laser-Induced Cavitation
%               
% Zhiren Zhu (zhiren@umich.edu)
% Dec. 2024
% =========================================================================
% Usage:
%
%   This is intended to be a clean update to the IMR simulation code
%   written by Jon Estrada leading up to his 2018 paper. (Ref [1])
%   The code structure is partially inspired by the code written by Jin
%   Yang during his post-doc at Wisconsin.
%
%   The present script searches for an optimized set of viscoelastic
%   parameters to match input data. An entire batch of experiments can be
%   treated.
%
%   Only do NHKV model for now.
%
% =========================================================================
% References:
%   [1] Estrada, Barajas, Henann, Johnsen & Franck (2018) JMPS 
%           (https://doi.org/10.1016/j.jmps.2017.12.006)
%   [2] Prosperetti, Crum & Commander (1988) JASA
%           (https://doi.org/10.1121/1.396145)
%   [3] Preston (2004) Ph.D. Thesis @ Caltech
%           (https://www.proquest.com/dissertations-theses/modeling-heat-mass-transfer-bubbly-cavitating/docview/305200526/se-2)
%   [4] Ascher & Greif, A First Course in Numerical Methods (SIAM)
%
% =========================================================================

clc; close all; clearvars;
addpath('../graphics');

%% Import Experiment

input_file = 'data/PA_10%_0.06%_WithRefined_T_c';
load(input_file);

% Unless specified otherwise, treat all experiments in the batch
n_samples = length(expts);

% Set up storage structure:
imr_batch_fit(n_samples) = struct;

for pick = 80:80 % n_samples

    % Announce which experiment, so that we can keep track:
    disp("Processing Experiment # " + pick + " ...")

    % Extract info:
    t_all = expts(pick).t;      % Time recorded (s)
    R_all = expts(pick).Roft;   % History of bubble radius (m)
    
    Rmax = expts(pick).R0; % Max. radius of bubble (m)
    Req = expts(pick).Req;  % Equilibrium radius (m)
    Lmax = Rmax/Req; % % Max. stretch of bubble (=Rmax/Req)


%% Define parameters

    % Note: do here, so that we can adjust the params structure later
    
    % (A) Far-field/atmosphere
    T_inf = 298.15;             % Far-field temperature (K)
    p_inf = 101325;             % Far-field pressure (Pa) 
    
    % (B) Surrounding material
    rho = 1020; % 998.2;                % Mass density (kg/m^3)
    c_wave = 1430;              % Longitudinal wave speed (m/s)
    gam = 7.0E-2;               % Surface tension (N/m)
    
    % (C) Bubble content
    kap = 1.4;                  % Heat capacity ratio
    D_bin = 2.42E-5;            % Binary diffusion coefficient (m^2/s)
    Ru = 8.3144598;             % Universal gas constant (J/mol*K)
    mwv = 18.01528E-3;          % Molecular weight of vapor (kg/mol)
    mwg = 28.966E-3;            % Molecular weight of non-vapor gas (kg/mol) [~ Air]
    Rv = Ru/mwv;                % Gas constant for vapor
    Rg = Ru/mwg;                % Gas constant for non-vapor gas
    
    % (C-i) Empirical constant for thermal conductivity 
    %   Per Ref [2], K = A*T + B works well for air at 200 to 3000 Kelvin
    A_thm = 5.28E-5;            % (W/m*K^2)
    B_thm = 1.17E-2;            % (W/m*K)
    
    % (C-ii) Empirical constant for saturated temperature of vapor
    %   Per Appendix A3 of Ref [3], pvsat(T) = p_ref*exp(-T_ref/T) 
    p_ref = 1.17E11;            % Reference pressure (Pa)
    T_ref = 5.2E3;              % Reference temperature (K)
    
%% Intermediate & non-dimensionalization parameters
    Req = Rmax/Lmax;                        % Equilibrium radius (m) [Not used for calculation]
    vc = sqrt(p_inf/rho);                   % Characteristic velocity (m/s)
    tc = Rmax/vc;                           % Characteristic time scale (s)
    K_inf = A_thm*T_inf + B_thm;            % Thermal conductivity (W/m*K)
    pv_inf = p_ref*exp(-T_ref/T_inf);       % Vapor pressure at equilibrium/far-field (Pa)
    
    c_star = c_wave/vc;                     % Non-dim of wave speed
    We = p_inf*Rmax/(2*gam);                % Weber number
    Fom = D_bin/(vc*Rmax);                  % Mass Fourier number
    chi = K_inf*T_inf/(p_inf*vc*Rmax);      % Lockhart-Martinelli number
    A_star = A_thm*T_inf/K_inf;             % Non-dim of A (= 1 - B-star)
    B_star = B_thm/K_inf;                   % Non-dim of B
    pv_star = pv_inf/p_inf;                 % Non-dim of far-field vapor pressure
    
    
%% Pick time span to simulate & package experiment data
    % (Define down here so that we can use characteristic scale flexibly)
    
    % Find optima:
    [max_data,min_data] = find_optima(t_all,R_all);
    tcut = max_data(3,1)/tc;  % Non-dimensional cut-off time
    tspan_star = 1.2*tcut; % Simulate a little more for fitting purpose 
    
    % Scale experiment data:
    nmin = find(t_all > 0, 1);
    nmax = find(t_all/tc > tcut, 1) - 1;
    ncount = nmax - nmin + 1;
    targ = [reshape(t_all(nmin:nmax),[ncount,1])/tc, reshape(R_all(nmin:nmax),[ncount,1])/Rmax];

%% Package parameters
    params = struct;
    
    params.We = We;
    params.Fom = Fom;
    params.chi = chi;
    params.kap = kap;
    params.A_star = A_star;
    params.B_star = B_star;
    params.c_star = c_star;
    params.pv_star = pv_star;
    params.Rv = Rv;
    params.Rg = Rg;
    params.Lmax = Lmax;
    params.tspan_star = tspan_star;
    
    % Extra dimensional info for this solver:
    params.Rmax = Rmax;
    params.rho = rho;
    params.p_inf = p_inf;

%% Run optimization

    % Define bounds:
    % Elastic shear modulus, G (Pa)
    G_start = 1.0E4;
    G_min = 1.0E2;
    G_max = 1.0E6;
    
    % Viscous shear modulus, mu (Pa*s)
    mu_start = 0.01;
    mu_min = 0.0;
    mu_max = 2.0;
    
    cost_fn = @(X) get_cost(X(1),X(2),params,targ);
    % search_options = optimset('OutputFcn', @outfun); % For debug - Add as an extra argument below
    
    tic
    search_sol = fminsearchbnd(cost_fn,[G_start,mu_start],[G_min,mu_min],[G_max,mu_max]);
    toc
    
    G_opt = search_sol(1);
    mu_opt = search_sol(2);
    
    % Display:
    disp("KV Best Fit: G = " + G_opt + " Pa, mu = " + mu_opt + " Pa*s.")
    
    % Run simulation again for quick check later:
    Ca_opt = p_inf/G_opt;                     
    Re_opt = sqrt(p_inf*rho)*Rmax/(mu_opt);    
    params.Ca = Ca_opt;
    params.Re = Re_opt;
    [tS,RS] = run_imr(params);
    tS_dim = tS*tc;
    RS_dim = RS*Rmax;
    
%% Output result

    imr_batch_fit(pick).G = G_opt;
    imr_batch_fit(pick).mu = mu_opt;
    imr_batch_fit(pick).exp_data = [t_all',R_all'];
    imr_batch_fit(pick).sim_data = [tS_dim,RS_dim];
    imr_batch_fit(pick).tcut = tcut*tc;

end

%% Save data
savename = ['batch-fit-test_PA10_',datestr(now,'yyyymmdd')];
save([savename,'.mat'],'imr_batch_fit');

%%
% =========================================================================
%                              SUBROUTINES
% =========================================================================

%% (I) Evaluate cost function

function psi = get_cost(G_elast, mu_visc, params, exp_data)

% Convert parameters to non-dim
p_inf = params.p_inf;
rho = params.rho;
Rmax = params.Rmax;

Ca = p_inf/G_elast;                     % Cauchy number
Re = sqrt(p_inf*rho)*Rmax/(mu_visc);    % Reynolds number

% Update parameters
params.Ca = Ca;
params.Re = Re;

% Run simulation:
[t_sim,R_sim] = run_imr(params);

% Refine simulation
div = 1E3;
dt = 1/div;
tmin = 0.0;
tmax = max(t_sim);
tmesh = (tmin:dt:tmax)';
RS_mesh = interp1(t_sim,R_sim,tmesh);
nmesh = length(tmesh);

% Loop through experiment points and identify closest point:
tX = exp_data(:,1);
RX = exp_data(:,2);
nX = length(tX);

psi = 0.0; % Total error

for ii = 1:nX

    delta = sum( ([tmesh,RS_mesh] - [tX(ii)*ones(nmesh,1),RX(ii)*ones(nmesh,1)]).^2, 2);
    psi = psi + min(delta);

end

psi = psi/nX; % Normalize by # of experiment points

end

%% (II) Simulation Driver

function [t_sol,R_sol] = run_imr(params)

%% Update parameters
% Parameters for computation
NY = 500;                   % # of mesh points for bubble content calculation
solver_RelTolX = 1E-7;      % Relative tolerance for ODE solver

params.NY = NY;

%% Extract parameters
pv_star = params.pv_star;
Lmax = params.Lmax;
We = params.We;
Rv = params.Rv;
Rg = params.Rg;
Ca = params.Ca;
tspan_star = params.tspan_star;

%% Initial conditions
Rb0 = 1;                                            % Bubble radius
Vb0 = 0;                                            % Bubble wall velocity
Theta0 = zeros(1,NY);                               % Thermal field - See Eqn.(52) of Ref [1] 
pb0 = pv_star + (1 + Lmax/We - pv_star)/(Lmax^3);   % Isobaric pressure
kv0_mag = 1/(1 + (Rv/Rg)*(pb0/pv_star - 1));        % Magnitude of vapor mass fraction - See Eqn.(56) of Ref [1] 
kv0 = kv0_mag*ones(1,NY);                           % Vapor mass fraction field
SI0 = -(5 - Lmax^(-4) - 4/Lmax)/(2*Ca);             % Stress integral

X0 = [Rb0,Vb0,pb0,SI0,Theta0,kv0]';

%% March Keller-Miksis forward

options = odeset('RelTol', solver_RelTolX);

[t_sol,X_sol] = ode23tb(@(t,X) KM_bubble(t,X,params), [0,tspan_star], X0, options);

%% Post-process

% Extract (dimensionless) solution:
R_sol = X_sol(:,1);

% Note: Output non-dim data directly

end

%% (III) Keller-Miksis steps
% Call this function to march governing equation forward

function dxdt = KM_bubble(t,X,params)

% READ INFO ================================================
% Unpack params:
NY = params.NY;
Ca = params.Ca;
Re = params.Re;
We = params.We;
Fom = params.Fom;
chi = params.chi;
kap = params.kap;
A_star = params.A_star;
B_star = params.B_star;
c_star = params.c_star;
pv_star = params.pv_star;
Rv = params.Rv;
Rg = params.Rg;
Lmax = params.Lmax;

% Current values:
Rb = X(1);      % Bubble radius
Vb = X(2);      % Bubble wall velocity
pb = X(3);      % Isobaric pressure
% SI = X(4);    % Stress integral. Evaluate directly from kinematics.

indx_T = 4;
indx_k = indx_T + NY;
Theta = X((indx_T + 1):(indx_T + NY));  % Thermal field ('theta')
kv = X((indx_k + 1):(indx_k + NY));     % Vapor mass fraction ('k')

% BUBBLE CONTENT ================================================

% Mesh for bubble content
dy = 1/(NY - 1);                % Uniform spacing
Ymesh = (dy*((1:NY)-1))';       % Dimensionless mesh from center to wall of bubble

% Dirichlet BC for vapor mass fraction
kv(end) = 1/(1 + (Rv/Rg)*(pb/pv_star - 1));

% Calculate mixture field
Tb = (-B_star + sqrt(1 + (2*A_star)*Theta) )./A_star;        % Temperature field (normalized by T_inf) - Integrate Eqn. (52) of Ref [1] to get quadratic equation
K_star = A_star.*Tb + B_star;      % Mixture thermal conductivity
Rmix = kv.*Rv + (1 - kv).*Rg;      % Mixture gas-constant field

% Find spatial derivatives
%       Neumann BC at center + Central difference for interior points + Backward difference at bubble wall
%       See: Ref [4] - Section 16.4 for backward difference formula. The formula here is BDF2.
%       For Laplacian in axisymmetric problem, nabla^2(f) = (1/r^2)*d(r^2*df/dr)/dr.
%           Expand: nabla^2(f) = d^2(f)/dr^2 + (2/r)*df/dr. Hence, the DD shown below.
%       For center of bubble, use L'Hopital to obtain Laplacian: 
%           nabla^2(f[0]) = 3 d^2(f[0])/dr^2 ~= (df[h/2]/dr - 0)/(h/2)  
%       Laplacian at bubble wall is not actually used, since we impose BC.
%           Equation shown uses (f'(1-h/2) - f'(1-h))/(h/2) for d^2(f[1])/dr^2.

BC1 = 0;
CD1= (Theta(3:end) - Theta(1:(end-2)))/(2*dy);             
BD1 = (3*Theta(end) - 4*Theta(end - 1) + Theta(end-2))/(2*dy);

DTheta = [BC1; CD1; BD1];

BC2 = (6*(Theta(2) - Theta(1)))/(dy^2);
CD2 = diff(diff(Theta))/(dy^2) + (2./Ymesh(2:(end-1))).*DTheta(2:(end-1));
BD2 = 0;        % (2*Theta(end) - 5*Theta(end-1) + 4*Theta(end-2) - Theta(end - 3))/(dy^2) + (2/Ymesh(end))*DTheta(end); % Not used
DDTheta = [BC2; CD2; BD2];

BC3 = 0;
CD3= (kv(3:end) - kv(1:(end-2)))/(2*dy);             
BD3 = (3*kv(end) - 4*kv(end - 1) + kv(end-2))/(2*dy);
Dkv = [BC3; CD3; BD3];

BC4 = (6*(kv(2) - kv(1)))/(dy^2);
CD4 = diff(diff(kv))/(dy^2) + (2./Ymesh(2:(end-1))).*Dkv(2:(end-1));
BD4 = 0;        % (2*kv(end) - 5*kv(end-1) + 4*kv(end-2) - kv(end - 3))/(dy^2) + (2/Ymesh(end))*Dkv(end); % Not used
DDkv = [BC4; CD4; BD4];

% Time rate of isobaric bubble pressure
%       See Eqn. (55) of Ref [1]
pb_part1 = -kap*pb*Vb + (kap-1)*(chi/Rb)*DTheta(end);
pb_part2 = kap*pb*Rv*Fom*Dkv(end)/(Rmix(end)*(1-kv(end))*Rb); 
pbdot = (3/Rb)*(pb_part1 + pb_part2);

% Mixture velocity field
%       See Eqn. (54) of Ref [1]
Vm_part1 = ((kap - 1)*(chi/Rb)*DTheta - Rb*pbdot*Ymesh/3)/(kap*pb);
Vm_part2 = Fom*(Rv - Rg)*Dkv./(Rb*Rmix);
Vmix = Vm_part1 + Vm_part2;

% Evolution of thermal and vapor mass fraction
%       See Eqn. (53) of Ref [1]
rhs_theta = pbdot + chi*DDTheta/(Rb^2) + (kap/(kap-1))*(pb./(K_star.*Tb))*(Fom/(Rb^2)).*((Rv-Rg)./Rmix).*Dkv.*DTheta;
denom_theta = (kap/(kap-1))*pb./(K_star.*Tb);
subtr_theta = (Vmix - Vb*Ymesh).*DTheta/Rb;
Theta_prime = rhs_theta./denom_theta - subtr_theta;
Theta_prime(end) = 0;   % Equilibrium imposed at boundary
   
rhomy = (1./Tb).*DTheta./sqrt(1 + 2*A_star*Theta) + (Rv-Rg)./Rmix.*Dkv; % = -(1/rho_m)*d(rho_m)/dy. Note: from Eqn. (20) of Ref [1], take spatial derivative to get this result.
rhs_kv = Fom/(Rb^2) * (DDkv - Dkv .* rhomy);
subtr_kv = (Vmix - Vb*Ymesh).*(Dkv)/Rb;
kv_prime = rhs_kv - subtr_kv;
kv_prime(end) = -(Rv/Rg)*(pbdot/pv_star)*(kv(end))^2;   % Equilibrium imposed at boundary (!= 0)

% SURROUNDING MATERIAL ================================================
Lb = Lmax*Rb;      % Current stretch of bubble
S_Ca = -(5 - Lb^(-4) - 4/Lb)/(2*Ca);
S_Re = -4*Vb/(Re*Rb);
SI = S_Ca + S_Re;

Sdot_Ca = (2/Ca)*(Vb/Rb)*(Lb^(-4) + 1/Lb); 
Sdot_Re1 = (4/Re)*(Vb/Rb)^2; % Acceleration-independent part
Sdot_Re2 = -(4/Re)/Rb; % Acceleration-dependent part, divided by acceleration
SIdot = Sdot_Ca + Sdot_Re1;

% INCREMENTS ================================================
Rdot = Vb;

KM_RHS = (1 + Vb/c_star) * (pb - 1/(We*Rb) + SI - 1) + (Rb/c_star)*(pbdot + SIdot + Vb/(We*Rb^2));
KM_LHS = (3/2)*(1 - Vb/(3*c_star))*(Vb^2);
KM_denom = (1 - Vb/c_star)*Rb - Sdot_Re2*(Rb/c_star);
Vdot = (KM_RHS - KM_LHS)/KM_denom;

SIdot = SIdot + Sdot_Re2*Vdot;  % Not actually used, but update properly

dxdt = [Rdot; Vdot; pbdot; SIdot; Theta_prime; kv_prime];

end

%% (IV) Plot history of bubble dynamics for comparison
function [] = plot_history(tX,RX,tS,RS,tcut,rgb1,rgb2,rgb3)

hold on; box on;

% Parameters:
to_micro = 1E6;
lw = 1;         % Line width
fs = 12;        % General font size
fs_ax = 16;     % Axis label font size

% Color options:
pb = [8 81 156]/255; % Blue
pr = [165,15,21]/255; % Red
pg = [35,139,69]/255; % Green
blk = [0,0,0]; % Black

% Assign default color
if nargin<6
    rgb1 = blk;
    rgb2 = pb;
    rgb3 = pg;
end

plot(tS*to_micro,RS*to_micro,'-','LineWidth',lw,'Color',rgb2)
plot(tX*to_micro,RX*to_micro,'s','LineWidth',lw,'Color',rgb1)
xline(tcut*to_micro,'--','LineWidth',lw*2,'Color',rgb3)

set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xl = xlabel('$t ({\rm \mu s})$');
yl = ylabel('$R ({\rm \mu m})$');
set(xl,'Interpreter','Latex','FontSize',fs_ax);
set(yl,'Interpreter','Latex','FontSize',fs_ax)

pbaspect([2, 1, 1])

end

%% (V) Output function during Nelder-Mead search
function stop = outfun(X,optimValues,~)
    stop = false;
    disp(X(1));
    disp(X(2));
    disp(optimValues.fval);
end

%% (VI) Find local minima and maxima of bubble dynamics
function [RMAX,RMIN] = find_optima(tX,RX)

nmax = 1;
nmin = 0;

ii = find(RX==max(RX)); % Assume we start at global max
tnow = tX(ii);
Rnow = RX(ii);

nX = length(tX);

big = 1E3;

RMAX = zeros(big,2); % [t, R]
RMIN = zeros(big,2); % [t, R]

% Initialize:
RMAX(1,:) = [tnow,Rnow];

find_max = 0; % Status

while ii < nX && min(nmin,nmax) < big

    if find_max == 0
        % Look for minimum
        slope = (RX((ii+1):end) - Rnow)./(tX((ii+1):end) - tnow);
        jj = find(slope == min(slope));
        jj = jj(1); % In case we have multiple results ...

        ii = ii + jj;

        if ii > nX
            % Don't update
            warning('Reached end while searching for min.');
        else
            tnow = tX(ii);
            Rnow = RX(ii);
            
            nmin = nmin + 1;
            RMIN(nmin,:) = [tnow,Rnow];
    
            dt = tnow - RMAX(nmin,1); 
        end
    else
        % Look for maximum
        scale = 1.5;
        icut = find(tX>(min(tnow + scale*dt,tX(end-1))),1);
        jj = find(RX((ii+1):icut) == max(RX((ii+1):icut)));
        jj = jj(1); % In case we have multiple results ...

        ii = ii + jj;

        if ii > nX
            % Don't update
            warning('Reached end while searching for max.');
        else
            tnow = tX(ii);
            Rnow = RX(ii);
    
            nmax = nmax + 1;
            RMAX(nmax,:) = [tnow,Rnow];
        end
    end

    find_max = 1 - find_max; % Switch

end

% Clean up:
RMIN = RMIN(1:nmin,:);
RMAX = RMAX(1:nmax,:);
end