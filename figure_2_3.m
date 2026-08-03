%% Figure 2.3: finite-domain level curves for different dish radii
% Panel group (a) uses the Dirac delta initial condition and panel group (b)
% uses the Heaviside initial condition.
%
% This script is self-contained. Running it produces Figure_2_3.pdf in the
% outputs folder alongside this file.

clear; close all; clc;

%% Parameters

aMIC = 0.02;
RdValues = [3 4 5 6];
RdCritical = 1/sqrt(pi*aMIC);

Nterms = 4000;
Nr = 320;
NtEarly = 180;
NtLate = 260;
tMin = 1e-4;
tMid = 0.25;
tMax = 5;
chunkSize = 250;

t = unique([linspace(tMin,tMid,NtEarly),linspace(tMid,tMax,NtLate)]);

%% Bessel data used by both initial conditions

alpha = besselj1_zeros(Nterms);
J0alpha = besselj(0,alpha);

%% Plot

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 23]);
tl = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(RdValues)
    Rd = RdValues(j);
    r = linspace(0,Rd,Nr);
    A = dirac_finite_grid(r,t,Rd,alpha,J0alpha,chunkSize);

    ax = nexttile(tl,j);
    contour(ax,r,t,A.',[aMIC aMIC], ...
        'LineWidth',1.5,'Color','red');
    format_panel(ax,Rd,tMax);
end

for j = 1:numel(RdValues)
    Rd = RdValues(j);
    r = linspace(0,Rd,Nr);
    A = heaviside_finite_grid(r,t,Rd,alpha,J0alpha,chunkSize);

    ax = nexttile(tl,j+4);
    contour(ax,r,t,A.',[aMIC aMIC], ...
        'LineWidth',1.5,'Color','green');
    format_panel(ax,Rd,tMax);
end

annotation(fig,'textbox',[0.48 0.505 0.04 0.03], ...
    'String','(a)','EdgeColor','none','HorizontalAlignment','center', ...
    'Interpreter','latex');
annotation(fig,'textbox',[0.48 0.012 0.04 0.03], ...
    'String','(b)','EdgeColor','none','HorizontalAlignment','center', ...
    'Interpreter','latex');

fprintf('Figure 2.3: 1/sqrt(pi*a_MIC) = %.6f\n',RdCritical);

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(fig,fullfile(outputDir,'Figure_2_3.pdf'), ...
    'ContentType','vector','Resolution',300);

%% Local functions

function A = dirac_finite_grid(r,t,Rd,alpha,J0alpha,chunkSize)
% Fourier--Bessel series for the finite-domain Dirac solution.
% The sum is evaluated in blocks to avoid forming a large 3-D array.

A = ones(numel(r),numel(t))/(pi*Rd^2);

for first = 1:chunkSize:numel(alpha)
    last = min(first+chunkSize-1,numel(alpha));
    a = alpha(first:last).';
    weights = 1./(J0alpha(first:last).'.^2);

    radialModes = besselj(0,r(:)*a/Rd);
    timeModes = exp(-t(:)*(a.^2)/Rd^2);

    A = A + ((radialModes.*weights)*timeModes.')/(pi*Rd^2);
end
end

function A = heaviside_finite_grid(r,t,Rd,alpha,J0alpha,chunkSize)
% Fourier--Bessel series for the finite-domain Heaviside solution.

A = ones(numel(r),numel(t))/(pi*Rd^2);

for first = 1:chunkSize:numel(alpha)
    last = min(first+chunkSize-1,numel(alpha));
    a = alpha(first:last).';
    J0 = J0alpha(first:last).';

    coefficients = (2./(a.*J0.^2)).*besselj(1,a/Rd);
    radialModes = besselj(0,r(:)*a/Rd);
    timeModes = exp(-t(:)*(a.^2)/Rd^2);

    A = A + ((radialModes.*coefficients)*timeModes.')/(pi*Rd);
end
end

function z = besselj1_zeros(N)
% First N positive zeros of J_1.

z = zeros(N,1);
for k = 1:N
    x0 = (k+0.25)*pi;
    left = x0-pi/3;
    right = x0+pi/3;

    % Expand the bracket in the unlikely event that the asymptotic guess
    % does not straddle a zero.
    attempts = 0;
    while sign(besselj(1,left)) == sign(besselj(1,right)) && attempts < 8
        left = left-pi/6;
        right = right+pi/6;
        attempts = attempts+1;
    end

    z(k) = fzero(@(x)besselj(1,x),[left right]);
end
end

function format_panel(ax,Rd,tMax)

box(ax,'on');
grid(ax,'on');
set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');
title(ax,sprintf('$R_d=%g$',Rd),'Interpreter','latex');
xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$t$','Interpreter','latex');
xlim(ax,[0 Rd]);
ylim(ax,[0 tMax]);
end
