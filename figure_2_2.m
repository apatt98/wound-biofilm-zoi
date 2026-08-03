%% Figure 2.2: comparison of the four antimicrobial level curves
% The level a(r,t) = a_MIC is plotted for both initial conditions on finite
% and infinite domains.
%
% This script is self-contained. Running it produces Figure_2_2.pdf in the
% outputs folder alongside this file.

clear; close all; clc;

%% Parameters

aMIC = 0.02;
Rd = 7;
Nr = 220;
Nt = 220;
Nterms = 4000;
Nl = 180;
tMin = 1e-3;
tMax = 5;

r = linspace(0,Rd,Nr);
t = linspace(tMin,tMax,Nt);
[R,T] = meshgrid(r,t);

%% Neumann eigenvalues

alpha = besselj1_zeros(Nterms);

%% Case 1: Dirac delta initial condition, infinite domain

A1 = (1./(4*pi*T)).*exp(-(R.^2)./(4*T));

%% Case 2: Heaviside initial condition, infinite domain

lq = linspace(0,1,Nl);
A2 = zeros(size(R));

for it = 1:Nt
    tt = t(it);
    RR = r(:);
    LL = lq(:).';

    integrand = (1/pi).*(LL./(2*tt)) ...
        .* exp(-(RR.^2 + LL.^2)./(4*tt)) ...
        .* besseli(0,(RR.*LL)./(2*tt));

    A2(it,:) = trapz(lq,integrand,2);
end

%% Case 3: Dirac delta initial condition, finite domain

A3 = (1/(pi*Rd^2))*ones(size(R));
for n = 1:Nterms
    an = alpha(n);
    modeR = besselj(0,an*r/Rd);
    modeT = exp(-(an^2)*t(:)/Rd^2);
    A3 = A3 + (1/(pi*Rd^2))*(1/besselj(0,an)^2)*(modeT*modeR);
end

%% Case 4: Heaviside initial condition, finite domain

A4 = (1/(pi*Rd^2))*ones(size(R));
for n = 1:Nterms
    an = alpha(n);
    coefficient = (1/(pi*Rd))*(2/(an*besselj(0,an)^2)) ...
        * besselj(1,an/Rd);
    modeR = besselj(0,an*r/Rd);
    modeT = exp(-(an^2)*t(:)/Rd^2);
    A4 = A4 + coefficient*(modeT*modeR);
end

%% Critical points for the infinite-domain cases

rZDirac = 1/sqrt(exp(1)*pi*aMIC);
tZDirac = 1/(4*exp(1)*pi*aMIC);

% The original calculation locates the Heaviside turning point directly
% from the computed contour. This is retained here.
CHeaviside = contourc(r,t,A2,[aMIC aMIC]);
[rZHeaviside,tZHeaviside] = max_r_from_contour(CHeaviside);

fprintf('Dirac infinite domain:     r_Z = %.6f, t_Z = %.6f\n', ...
    rZDirac,tZDirac);
fprintf('Heaviside infinite domain: r_Z = %.6f, t_Z = %.6f\n', ...
    rZHeaviside,tZHeaviside);

%% Plot

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 12]);
ax = axes(fig);
hold(ax,'on');
box(ax,'on');
grid(ax,'on');
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
    sprintf('$r_Z=%.3f,\\ t_Z=%.3f$',rZDirac,tZDirac), ...
    'Interpreter','latex','FontSize',12,'BackgroundColor','w','Margin',2);
text(ax,rZHeaviside+0.08,tZHeaviside-0.18, ...
    sprintf('$r_Z=%.3f,\\ t_Z=%.3f$',rZHeaviside,tZHeaviside), ...
    'Interpreter','latex','FontSize',12,'BackgroundColor','w','Margin',2);

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
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(fig,fullfile(outputDir,'Figure_2_2.pdf'), ...
    'ContentType','vector','Resolution',300);

%% Local functions

function z = besselj1_zeros(N)
% First N positive zeros of J_1.

z = zeros(N,1);
for k = 1:N
    x0 = (k+0.25)*pi;
    z(k) = fzero(@(x)besselj(1,x),[x0-pi/3,x0+pi/3]);
end
end

function [rMax,tAtMax] = max_r_from_contour(C)
% Return the contour point having the largest radial coordinate.

if isempty(C)
    error('No contour was found at the requested level.');
end

column = 1;
rMax = -Inf;
tAtMax = NaN;

while column < size(C,2)
    nPoints = C(2,column);
    points = C(:,column+1:column+nPoints);
    [candidate,index] = max(points(1,:));

    if candidate > rMax
        rMax = candidate;
        tAtMax = points(2,index);
    end

    column = column+nPoints+1;
end

if ~isfinite(rMax)
    error('The maximum-radius contour point could not be extracted.');
end
end
