%%%%%%%%%%%%%%% Transition Dynamics in the Aiyagari Model %%%%%%%%%%%%%%%%%
% Autor: Yanran Guo
% Date:  10/8/2019
% Task:  I first replicate the steady-state of the Aiyagari (1994) model, 
%        and then analyze the model's transition dynamics after an
%        unexpected shock hits the economy. 

% An "MIT shock" is an unexpected shock that hits an economy at its 
% steady state, leading to a transition path back towards the economy's steady state.

% I compute the steady state before/after the shock and the transition path.
% In particular, the shock I'm using as an example here is a change in
% borrowing constraint.


%% Environment Setup
clc;                     % Clear screen
clear;                   % Clear memory
close all;               % Close open windows
addpath(genpath(pwd));   % Include subfolders of the current folder
          
                          
%% Parameterization
cS        = SetParameterValues;
cS.dbg    = 1;            % dbg is debugging parameter
                          % dbg=1, all functions will debug automatically


%% Probability Transition Matrices
[cS.Pi, cS.Z] ...
          = StationaryDis_MarkovProcess(cS.s, cS.P, cS.dbg); 


%% Solve for Stationary Equilibrium
% Compute the initial steady state, i.e. before the shock
SS1       = Aiyagari_Calib(cS);

% An unexpected change in the borrowing constraint
% e.g. debt limit in new steady state: k'>=0
cS.phi    = 0; 
% Have to change the grid vector for k because it starts from phi
cS.kMin   = cS.phi; 
cS.kGridV = linspace(cS.kMin, cS.kMax, cS.nk);  
cS.kGridV = cS.kGridV';
% Compute the new steady state, i.e. after the shock
SS2       = Aiyagari_Calib(cS);


%% Transition Dynamics: MIT Shock to Borrowing Constraint
% Find the sequence of prices, optimal decistions, distributions between
% two stationary equilibria (SE)

%------------------------- Define input -----------------------------------
% Initial stationary equilibrium
Ts.kPolM_initial  = SS1.kPolM;
Ts.Ind_initial    = SS1.IndM;
Ts.valueM_initial = SS1.valueM;
Ts.mu_initial     = SS1.mu;
Ts.r_initial      = SS1.r;
Ts.K_initial      = SS1.wealth;

% End stationary equilibrium
Ts.kPolM_end      = SS2.kPolM;
Ts.Ind_end        = SS2.IndM;
Ts.valueM_end     = SS2.valueM;
Ts.mu_end         = SS2.mu;
Ts.r_end          = SS2.r;
Ts.K_end          = SS2.wealth;

% Others
Ts.phi            = 0;                  % Borrowing constraint in new stationary equilibrium
Ts.kMin           = Ts.phi;             % Due to the borrowing constraint: k >= phi
Ts.kGridV         = linspace(Ts.kMin, cS.kMax, cS.nk); 
Ts.kGridV         = Ts.kGridV';
Ts.T              = 50;                 % Transition periods between two stationary equilibria


%--------------------- Run Transition Dynamics Routine --------------------
Seq               = TransitionDynamics(cS, Ts);


%----------------------- Plot the Transition Results ----------------------
% Plots interest rates during the transition period
t                 = 0 : 1 : Ts.T;
r                 = cell2mat({Seq([1:size(Seq,1)]).r}');
r                 = [Ts.r_initial; r; Ts.r_end];

figure
plot(t, r, 'red', 'LineWidth', 1)
ylabel('Equilibrium Interest Rate', 'FontSize', 11)
xlabel('Period t', 'FontSize', 11)
title('Transition Dynamics of Interest Rate', 'FontSize', 11)
set(gca, 'FontSize', 11, 'LineWidth', 1, 'Box', 'on', 'FontName', 'Times New Roman');
grid on
dim = [0.3 0.05 0.3 0.3];
str = {'t=0 is the initial stationary equilibrium, r1=0.03687', ...
       't=50 is the new stationary equilibrium, r2=0.03684'};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

