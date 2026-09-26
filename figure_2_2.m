%% Figure 2.2: finite-domain level curves for different domain radii
% This script examines how the boundary of the antimicrobial inhibition
% region depends on the size of a finite radial domain.
%
% A Dirac delta initial condition is used for the antimicrobial
% concentration. For each domain radius R_d, the diffusion equation is
% evaluated using a Fourier--Bessel expansion with a no-flux boundary
% condition at r = R_d.
%
% The plotted curves show the locations in the (r,t) plane where the
% antimicrobial concentration is equal to the minimum inhibitory
% concentration a_MIC.
%
% Running this script produces:
%   outputs/Figure_2_2.pdf

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

%% Bessel data

alpha = besselj1_zeros(Nterms);
J0alpha = besselj(0,alpha);

%% Plot

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 12]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

for j = 1:numel(RdValues)
    Rd = RdValues(j);
    r = linspace(0,Rd,Nr);
    A = dirac_finite_grid(r,t,Rd,alpha,J0alpha,chunkSize);

    ax = nexttile(tl);
    contour(ax,r,t,A.',[aMIC aMIC], ...
        'LineWidth',1.5,'Color','red');
    format_panel(ax,Rd,tMax);
end

fprintf('Figure 2.2: 1/sqrt(pi*a_MIC) = %.6f\n',RdCritical);

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(fig,fullfile(outputDir,'Figure_2_2.pdf'), ...
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
