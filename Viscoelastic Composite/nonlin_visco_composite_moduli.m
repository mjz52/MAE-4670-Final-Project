% Michael Zakoworotny
% Plot a viscoelastic composite with a "hard" and "soft" nonlinear response

clear all; clc; close all;
% formatMike(); % plot formatting

% Hard PDMS Material constants
mu_hard = 0.6563; Jm_hard = 2.537; % Hard Hyperelastic spring 1
Eb_hard = 2; % Linear spring 2
eta_hard = 0.5; % Dashpot
% Soft PDMS Material constants
mu_soft = 0.042; Jm_soft = 30; % Soft Hyperelastic spring 1
Eb_soft = 0.1;
eta_soft = 0.2;

% Fiber volume fraction
Vf = 0.5;

% Loading conditions
e_max = 0.4; % max strain
e_freq = 2*pi/192; % frequency - period of 200 s
strain_model = @(t)sin_loading(t, e_max, e_freq); % use sinusoidal strain

%%% HARD PDMS STRESS
% Visco model component functions
sigAFun_h = @(e) sigAModel(e+1, mu_hard, Jm_hard);
sigDotAFun_h = @(e,eDot) sigDotAModel(e+1, eDot, mu_hard, Jm_hard);
sigBFun_h = @(e) sigBModel(e, Eb_hard);
sigDotBFun_h = @(eDot)sigDotBModel(eDot, Eb_hard); % may also make a function of eB if generalized
sigCFun_h = @(eDot)sigCModel(eDot, eta_hard);
mat_models_hard = {sigAFun_h, sigDotAFun_h, sigBFun_h, sigDotBFun_h, sigCFun_h};
%%% SOFT PDMS STRESS
% Visco model component functions
sigAFun_s = @(e) sigAModel(e+1, mu_soft, Jm_soft);
sigDotAFun_s = @(e,eDot) sigDotAModel(e+1, eDot, mu_soft, Jm_soft);
sigBFun_s = @(e) sigBModel(e, Eb_soft);
sigDotBFun_s = @(eDot)sigDotBModel(eDot, Eb_soft); % may also make a function of eB if generalized
sigCFun_s = @(eDot)sigCModel(eDot, eta_soft);
mat_models_soft = {sigAFun_s, sigDotAFun_s, sigBFun_s, sigDotBFun_s, sigCFun_s};
% Simulation
t_end = 500;
% Numerically integrate to obtain total stress over time
sDotFun = @(t,s)getSDot(t,s,mat_models_hard,mat_models_soft,strain_model);
tol = 1e-8; opts = odeset('RelTol',tol,'AbsTol',tol,'Events',@(t,s)cycle_event(t,s,strain_model));
[t_array,s_array] = ode45(sDotFun,[0,t_end],[0,0]',opts);

% Get corresponding stretches to plot with
[e_array,~] = strain_model(t_array);
lam_array = e_array + 1;

% Get pure Gent model curve - no viscoelasticity
s_fit_hard = sigAFun_h(sort(e_array));
s_fit_soft = sigAFun_s(sort(e_array));

% Composite stress
s_comp_array = Vf*s_array(:,1) + (1-Vf)*s_array(:,2);

% initial_slope = sDotFun(0,0)/e_rate; % small strain slope - should approach 

% Plot stress vs strain
figure; box on; hold on;
plot(lam_array, s_comp_array,'-k');
plot(sort(e_array)+1, s_fit_hard, '--b','Linewidth',1);
plot(lam_array,s_array(:,1),'-b','Linewidth',1);
plot(sort(e_array)+1, s_fit_soft, '--r','Linewidth',1);
plot(lam_array,s_array(:,2),'-r','LineWidth',1);
plot([min(xlim),max(xlim)],[0,0],'-k','LineWidth',0.5);
ylim([-0.5,max(ylim)]);
xlabel('$\lambda$'); ylabel('Nominal stress (MPa)');
legend({'Visco Comp.','Pure - Hard','Visco - Hard','Pure - Soft','Visco - Soft'},'location','nw');
% Plot applied strain and response
figure;
colororder({['k'],['b']});
yyaxis left
plot(t_array,lam_array);
ylabel('$\lambda$');
yyaxis right
plot(t_array,s_comp_array);
xlim([0,max(t_array)]);
ylabel('$\sigma_{comp}$ (MPa)');
xlabel('t (s)');

% Obtain viscoelastic moduli
% Storage modulus
t_stor = pi/2/e_freq;
s_stor = interp1(t_array,s_comp_array,t_stor);
E_stor = s_stor / e_max;
% Loss modulus
t_loss = 2*pi/e_freq;
s_loss = interp1(t_array,s_comp_array,t_loss);
E_loss = s_loss / e_max;
% delta
delta = atan(E_loss/E_stor);
% Obtain tan delta from plot
[~,stress_max_ind] = findpeaks(abs(s_comp_array));
t_round = (pi*(floor(e_freq*t_array(stress_max_ind(end))/pi))+pi/2)/e_freq;
delta_plot = (t_round - t_array(stress_max_ind(end)))*e_freq;

% Printout moduli
fprintf('Storage modulus: %.4f MPa\n',E_stor);
fprintf('Loss modulus: %.3f MPa\n',E_loss);
fprintf('Tan delta: %.4f\n',tan(delta));
fprintf('Graphical tan delta: %.4f\n',tan(delta_plot));


%%%%%% VISCOELASTIC MODEL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the propagation function for stress - return stress rate in fiber
% and matrix
function sDot = getSDot(t, s_comp, mat_models_hard, mat_models_soft, strain_model)
    % Propagate both hard and soft values
    s_hard = s_comp(1); s_soft = s_comp(2);
    % Implement hard viscoelastic 
    sigAFun = mat_models_hard{1};
    sigDotAFun = mat_models_hard{2};
    sigBFun = mat_models_hard{3};
    sigDotBFun = mat_models_hard{4};
    sigCFun = mat_models_hard{5};
    
    [e, eDot] = strain_model(t);
    
    sDotA = sigDotAFun(e,eDot);
    sA = sigAFun(e);
    sC = s_hard - sA;
    eDotC = fzero(@(eDotC)sigCFun(eDotC) - sC, 0.25);
    eDotB = eDot - eDotC;
    sDotB = sigDotBFun(eDotB);
    sDot_hard = sDotA + sDotB;
    
    % Implement soft viscoelastic 
    sigAFun = mat_models_soft{1};
    sigDotAFun = mat_models_soft{2};
    sigBFun = mat_models_soft{3};
    sigDotBFun = mat_models_soft{4};
    sigCFun = mat_models_soft{5};
    
    [e, eDot] = strain_model(t);
    
    sDotA = sigDotAFun(e,eDot);
    sA = sigAFun(e);
    sC = s_soft - sA;
    eDotC = fzero(@(eDotC)sigCFun(eDotC) - sC, 0.25);
    eDotB = eDot - eDotC;
    sDotB = sigDotBFun(eDotB);
    sDot_soft = sDotA + sDotB;
    
    sDot = [sDot_hard, sDot_soft]';
end

function sigA = sigAModel(lam, mu, Jm)
    % Hyperelastic response of viscoelastic model
    I1 = lam.^2 + 2./lam;
    sigA = mu*Jm./(Jm-(I1-3)).*(lam.^2 - 1./lam)./lam;
    
end

function sigDotA = sigDotAModel(lam, lamDot, mu, Jm)
    % Derivative of stress in hyperelastic spring
    I1 = lam.^2 + 2./lam;
    
    DsigDlam = mu*Jm./(Jm-(I1-3)).*((2*lam-2./lam.^2)*(lam-1./lam.^2)/(Jm-(I1-3)) + (1+2./lam.^3));
    sigDotA = DsigDlam .* lamDot;
    
end

function sigB = sigBModel(strainB, Eb)
    % Elastic spring in Maxwell element of viscoelastic model

    sigB = Eb*strainB;
end

function sigDotB = sigDotBModel(strainDotB, Eb)
    % Derivative of elastic spring
    
    sigDotB = Eb*strainDotB;
end

function sigC = sigCModel(strainDotC, eta)
    % Viscous response

    sigC = eta*strainDotC;
end

% Given a time (scalar or vector), a maximum strain, and a strain rate,
% return the current strain and strain rate using a cyclic (triangular)
% loading
function [e,eDot] = cyclic_loading(t, e_max, e_rate)
    t_half = e_max/e_rate; % half period
    
    e = zeros(1,length(t)); eDot = e; % Initialize
    % Get indices for rising and falling
    rising = find(mod(floor(t/t_half),2) == 0);
    falling = find(mod(floor(t/t_half),2) == 1);
    % Strains in rising and falling portions
    e(rising) = e_rate*mod(t(rising),t_half);
    e(falling) = e_rate*(t_half - mod(t(falling),t_half));
    % Strain rates in rising and falling portions
    eDot(rising) = e_rate;
    eDot(falling) = -e_rate;

end

% Given a time (scalar or vector), a strain amplitude and a frequency,
% return the current strain and strain rate using a sinusoidal loading
function [e,eDot] = sin_loading(t, e_max, e_freq)
    e = e_max*sin(e_freq*t);
    eDot = e_max*e_freq*cos(e_freq*t);
end

function [value,isterminal,direction] = cycle_event(t,s,strain_model)
    % End ode45 integration when stress reaches 0 again
    value = strain_model(t);
    isterminal = 0; % set to 1 to stop simulation when strain completes one period
    direction = 1;
end



