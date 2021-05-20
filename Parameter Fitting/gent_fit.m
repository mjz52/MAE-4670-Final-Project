% Michael Zakoworotny
% Obtain the Gent material parameters from the data in Wang et al

clear all; clc; close all;
% formatMike(); % plot formatting
set(0,'DefaultLineLineWidth',2);

%%%%%%%% Plot representative model of Gent and Neo Hookean %%%%%%%%%%%%%%%
mu = 1;
Jm = 1;
stretch_lock = fzero(@(x)x^3 - (Jm+3)*x + 2, 10);

L1 = linspace(1,stretch_lock,1000);
L2 = 1./sqrt(L1); L3 = L2;
L = [L1; L2; L3];

I1 = sum(L.^2);
sig = mu*Jm./(Jm-(I1-3)) .* (L.^2 - 1./L1)./L; % nominal stress
s1 = sig(1,:);

figure; box on; hold on;
plot(L1, s1, '-k');
plot(L1, mu*(L1-1./L1.^2), '-r'); % Neo-Hookean
plot(L1, 3*mu*(L1-1), '--b'); % Initial slope
xlim([1,stretch_lock]); ylim([0,(stretch_lock-1)*mu*10]);
legend({'Gent','NH','$3\mu$ slope'})

%%%%%%%%%%%%%%%% EXPERIMENTAL FITTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import hard experimental data
hard_data = xlsread('stress_strain.xlsx','HardPDMS');
L_hard_exp = hard_data(:,1); s_hard_exp = hard_data(:,2);

% Produce fit to experimental data - hard
diff_fun = @(x) sum(abs(getS1(L_hard_exp,x(1),x(2)) - s_hard_exp));
[mat_fit] = fminsearch(diff_fun, [1,2]);
% [mat_fit] = fmincon(diff_fun,[1,1],[],[],[],[],[0,0],[20,20]);
mu_hard_fit = mat_fit(1); Jm_hard_fit = mat_fit(2);
L_hard_fit = linspace(1,max(L_hard_exp),100);
s_hard_fit = getS1(L_hard_fit,mu_hard_fit,Jm_hard_fit);

% Import hard experimental data
soft_data = xlsread('stress_strain.xlsx','SoftPDMS');
L_soft_exp = soft_data(:,1); s_soft_exp = soft_data(:,2);

% Produce fit to experimental data - soft
diff_fun = @(x) sum(abs(getS1(L_soft_exp,x(1),x(2)) - s_soft_exp));
[mat_fit] = fminsearch(diff_fun, [0.05,5]);
% [mat_fit] = fmincon(diff_fun,[1,1],[],[],[],[],[0,0],[20,20]);
mu_soft_fit = mat_fit(1); Jm_soft_fit = mat_fit(2);
L_soft_fit = linspace(1,max(L_soft_exp),100);
s_soft_fit = getS1(L_soft_fit,mu_soft_fit,Jm_soft_fit);

% Plot curves for experimentla data and git
fig = figure; box on; hold on;
plot(L_hard_fit, s_hard_fit, '-r');
plot(L_hard_exp, s_hard_exp, 'x--k','LineWidth',1,'MarkerSize',5);
plot(L_soft_fit, s_soft_fit, '-b');
plot(L_soft_exp, s_soft_exp, 'x--c','LineWidth',1,'MarkerSize',5);
text(1.9,2.25,{['$\mu_h=$',num2str(mu_hard_fit,4),' MPa'],[''],['$J_{m,h}=$',num2str(Jm_hard_fit,4)]});
text(1.75,0.5,{['$\mu_s=$',num2str(mu_soft_fit,3),' MPa'],[''],['$J_{m,s}=$',num2str(Jm_soft_fit,4)]});
xlim([1,max([L_hard_exp;L_soft_exp])]); ylim([0,1.1*max(s_hard_exp)]);
legend({'Hard Fit','Hard Exp.','Soft Fit','Soft Exp.'},'Location','nw');
xlabel('$\lambda$'); ylabel('Nominal Stress (MPa)');
% saveas(fig,'gent_fit.png');

% Get work done
W_fit = trapz(L_hard_fit, s_hard_fit);

% Paper values
mu_hard_pa = 0.59;
Jm_hard_pa = 2.9;
mu_soft_pa = 0.042;
Jm_soft_pa = 30;

%%%%%%%%%%%%%% COMPOSITIVE CONSTITUTIVE RESPONSES %%%%%%%%%%%%%%%%%%%%%%%%%
% 10% composite (10% hard to 90% soft)
data_10 = xlsread('stress_strain.xlsx','10p_PDMS');
L_10_exp = data_10(:,1); s_10_exp = data_10(:,2);
Vf = 0.1;
L_10 = L_hard_fit;
s_10 = Vf*getS1(L_10,mu_hard_fit,Jm_hard_fit) + (1-Vf)*getS1(L_10,mu_soft_pa,Jm_soft_pa);

% 50% composite (50% hard to 50% soft)
data_50 = xlsread('stress_strain.xlsx','50p_PDMS');
L_50_exp = data_50(:,1); s_50_exp = data_50(:,2);
Vf = 0.5;
L_50 = L_hard_fit;
s_50 = Vf*getS1(L_50,mu_hard_fit,Jm_hard_fit) + (1-Vf)*getS1(L_50,mu_soft_pa,Jm_soft_pa);

% Plot composite responses
fig = figure; box on; hold on;
plot(L_50_exp,s_50_exp,'--','color',getColor('LimeGreen'));
plot(L_50,s_50,'-g');
plot(L_10_exp,s_10_exp,'--','color',getColor('salmon'));
plot(L_10,s_10,'-m');
legend({'50\% Exp.','50\% Fit','10\% Exp.','10\% Fit'},'Location','nw');
xlabel('$\lambda$'); ylabel('Nominal Stress (MPa)');
saveas(fig,'composite_fit.png');



function s1 = getS1(L,mu,Jm)
    % Given lambda_1, reutrn sigma 1, using the Gent nominal stress
    
    s1 = mu*Jm./(Jm-(L.^2+2./L-3)).*(L - 1./L.^2);
end






