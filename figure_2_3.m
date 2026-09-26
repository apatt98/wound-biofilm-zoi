%% Figure 2.3: comparison of antimicrobial MIC level curves
% This script compares the antimicrobial inhibition boundary for four
% combinations of initial condition and spatial domain:
%
%   1. Dirac delta initial condition on an infinite domain
%   2. Heaviside initial condition on an infinite domain
%   3. Dirac delta initial condition on a finite domain
%   4. Heaviside initial condition on a finite domain
%
% In each case, the plotted curve is the level set a(r,t) = a_MIC. The
% maximum inhibition radius and corresponding stopping time are also
% identified for the two infinite-domain solutions.
%
% Running this script produces:
%   outputs/Figure_2_3.pdf

clear; close all; clc;

%% Parameters
aMIC = 0.02;
Rd = 7;
Nr = 320;
Nt = 360;
Nterms = 4000;
Nl = 240;
tMin = 1e-3;
tMax = 5;

r = linspace(0,Rd,Nr);
% Extra resolution at early times without making the plotted grid irregular.
t = unique([logspace(log10(tMin),log10(0.25),140), ...
            linspace(0.25,tMax,Nt)]);
[R,T] = meshgrid(r,t);

%% Neumann eigenvalues
alpha = besselj1_zeros(Nterms);

%% Case 1: Dirac delta initial condition, infinite domain
A1 = (1./(4*pi*T)).*exp(-(R.^2)./(4*T));

%% Case 2: Heaviside initial condition, infinite domain
lq = linspace(0,1,Nl);
A2 = zeros(size(R));

for it = 1:numel(t)
    tt = t(it);
    RR = r(:);
    LL = lq(:).';
    z = (RR.*LL)./(2*tt);

    % besseli(0,z,1) = exp(-abs(z))*besseli(0,z). Since z >= 0,
    % this is algebraically identical to the original integrand but stable.
    integrand = (1/pi).*(LL./(2*tt)) ...
        .* exp(-((RR-LL).^2)./(4*tt)) ...
        .* besseli(0,z,1);

    A2(it,:) = trapz(lq,integrand,2);
end

%% Case 3: Dirac delta initial condition, finite domain
A3 = (1/(pi*Rd^2))*ones(size(R));
for n = 1:Nterms
    an = alpha(n);
    modeR = besselj(0,an*r/Rd);
    modeT = exp(-(an^2)*t(:)/Rd^2);
    A3 = A3 + (1/(pi*Rd^2))*(1/besselj(0,an)^2) ...
        *(modeT*modeR);
end

%% Case 4: Heaviside initial condition, finite domain
A4 = (1/(pi*Rd^2))*ones(size(R));
for n = 1:Nterms
    an = alpha(n);
    coefficient = (1/(pi*Rd))*(2/(an*besselj(0,an)^2)) ...
        *besselj(1,an/Rd);
    modeR = besselj(0,an*r/Rd);
    modeT = exp(-(an^2)*t(:)/Rd^2);
    A4 = A4 + coefficient*(modeT*modeR);
end

%% Infinite-domain turning points
rZDirac = 1/sqrt(exp(1)*pi*aMIC);
tZDirac = 1/(4*exp(1)*pi*aMIC);

% Solve the Heaviside level and turning-point conditions simultaneously.
% Log variables enforce r_Z > 0 and t_Z > 0.
options = optimoptions('fsolve','Display','none', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12);
logRT = fsolve(@(x) heaviside_turning_system(x,aMIC), ...
    log([rZDirac;tZDirac]),options);
rZHeaviside = exp(logRT(1));
tZHeaviside = exp(logRT(2));

fprintf('Dirac infinite domain:     r_Z = %.9f, t_Z = %.9f\n', ...
    rZDirac,tZDirac);
fprintf('Heaviside infinite domain: r_Z = %.9f, t_Z = %.9f\n', ...
    rZHeaviside,tZHeaviside);

%% Plot
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 12]);
ax = axes(fig);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontName','Times','FontSize',12,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');

[~,h1] = contour(ax,R,T,A1,[aMIC aMIC], ...
    'LineWidth',2,'Color','blue','LineStyle','-');
[~,h2] = contour(ax,R,T,A2,[aMIC aMIC], ...
    'LineWidth',2,'Color','magenta','LineStyle','-');
[~,h3] = contour(ax,R,T,A3,[aMIC aMIC], ...
    'LineWidth',2,'Color','red','LineStyle','--');
[~,h4] = contour(ax,R,T,A4,[aMIC aMIC], ...
    'LineWidth',2,'Color','green','LineStyle','--');

p1 = plot(ax,rZDirac,tZDirac,'o','MarkerSize',9, ...
    'MarkerFaceColor','blue','MarkerEdgeColor','k','LineWidth',1.5);
p2 = plot(ax,rZHeaviside,tZHeaviside,'d','MarkerSize',9, ...
    'MarkerFaceColor','magenta','MarkerEdgeColor','k','LineWidth',1.5);

text(ax,rZDirac+0.08,tZDirac+0.08, ...
    sprintf('$r_Z=%.4f,\\;t_Z=%.4f$',rZDirac,tZDirac), ...
    'Interpreter','latex','FontSize',10,'BackgroundColor','w','Margin',2);
text(ax,rZHeaviside+0.08,tZHeaviside-0.18, ...
    sprintf('$r_Z=%.4f,\\;t_Z=%.4f$',rZHeaviside,tZHeaviside), ...
    'Interpreter','latex','FontSize',10,'BackgroundColor','w','Margin',2);

xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$t$','Interpreter','latex');
xlim(ax,[0 3]);
ylim(ax,[0 tMax]);
legend(ax,[h1 h2 h3 h4 p1 p2], ...
    {'Case 1: Dirac delta IC on infinite domain', ...
     'Case 2: Heaviside IC on infinite domain', ...
     'Case 3: Dirac delta IC on finite domain', ...
     'Case 4: Heaviside IC on finite domain', ...
     'Case 1: $(r_Z,t_Z)$', ...
     'Case 2: $(r_Z,t_Z)$'}, ...
    'Location','southwest','FontSize',10,'Interpreter','latex','Box','on');

%% Export
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
exportgraphics(fig,fullfile(outputDir,'Figure_2_3.pdf'), ...
    'ContentType','vector','Resolution',300);

%% Local functions
function F = heaviside_turning_system(logRT,aMIC)
    rZ = exp(logRT(1));
    tZ = exp(logRT(2));
    z = rZ/(2*tZ);

    % a_H(r,t) = pi^{-1}[1-Q_1(r/sqrt(2t),1/sqrt(2t))].
    levelResidual = (1/pi)*(1-marcumq( ...
        rZ/sqrt(2*tZ),1/sqrt(2*tZ))) - aMIC;
    turningResidual = rZ*besseli(1,z,1)-besseli(0,z,1);
    F = [levelResidual;turningResidual];
end

function z = besselj1_zeros(N)
    z = zeros(N,1);
    for k = 1:N
        x0 = (k+0.25)*pi;
        z(k) = fzero(@(x)besselj(1,x),[x0-pi/3,x0+pi/3]);
    end
end
