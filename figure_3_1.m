%% Figure 3.1: bacterial density profiles
% This script examines the bacterial density predicted by the model after
% exposure to a diffusing antimicrobial.
%
% Figures 3.1a and 3.1b show bacterial density profiles at selected times
% for a fixed bacterial death-rate parameter. The two panels correspond to
% Dirac delta and Heaviside antimicrobial initial conditions, respectively.
%
% Figures 3.1c and 3.1d show bacterial density profiles at a fixed time for
% several values of the bacterial death-rate parameter k. Again, the two
% panels correspond to Dirac delta and Heaviside antimicrobial initial
% conditions.
%
% Bacterial death occurs only where the antimicrobial concentration exceeds
% the minimum inhibitory concentration a_MIC. The density at each radius is
% calculated from the accumulated antimicrobial exposure.
%
% Running this script produces:
%   outputs/Figure_3_1a.pdf
%   outputs/Figure_3_1b.pdf
%   outputs/Figure_3_1c.pdf
%   outputs/Figure_3_1d.pdf

clear; close all; clc;

%% Figure 3.1a and Figure 3.1b: profiles at selected times
%% Parameters

rMax = 10;
rZ = 4;
Nr = 800;
r = linspace(0,rMax,Nr);
r(1) = 1e-6;
b0 = double(r >= 1);

k = 50;
timeFractions = [0 0.5 1 5];

%% Stopping times and threshold concentrations

aDirac = @(r,t) (1./(4*pi*t)).*exp(-(r.^2)./(4*t));

tZDirac = rZ^2/4;

g = @(y) rZ*besseli(1,rZ./(2*exp(y))) ...
    - besseli(0,rZ./(2*exp(y)));
options = optimoptions('fsolve','Display','none', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);
tZHeaviside = exp(fsolve(g,log(tZDirac),options));

aMICDirac = aDirac(rZ,tZDirac);
aMICHeaviside = (1/pi)*(1-marcumq( ...
    rZ/sqrt(2*tZHeaviside),1/sqrt(2*tZHeaviside)));

timesDirac = timeFractions*tZDirac;
timesHeaviside = timeFractions*tZHeaviside;

fprintf('Dirac:     t_Z = %.6f, a_MIC = %.6e\n',tZDirac,aMICDirac);
fprintf('Heaviside: t_Z = %.6f, a_MIC = %.6e\n', ...
    tZHeaviside,aMICHeaviside);

%% Plot styles

styles = { ...
    struct('Color',[0.00 0.45 0.74],'LineStyle','-'), ...
    struct('Color',[0.85 0.33 0.10],'LineStyle','--'), ...
    struct('Color',[0.47 0.67 0.19],'LineStyle','-.'), ...
    struct('Color',[0.49 0.18 0.56],'LineStyle',':')};

%% Figure 3.1a: Dirac initial condition

figA = figure('Color','w','Units','centimeters','Position',[2 2 11 9]);
ax = axes(figA);
hold(ax,'on');
for j = 1:numel(timesDirac)
    if timesDirac(j) == 0
        bj = b0;
    else
        bj = b_snapshot_time(1,aMICDirac,k,r,timesDirac(j),b0);
    end

    plot(ax,r,bj,'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_density_axes_31(ax,rMax);
legend(ax,time_labels_31(timeFractions), ...
    'Interpreter','latex','Location','southeast','FontSize',10,'Box','on');

%% Figure 3.1b: Heaviside initial condition

figB = figure('Color','w','Units','centimeters','Position',[14 2 11 9]);
ax = axes(figB);
hold(ax,'on');
for j = 1:numel(timesHeaviside)
    if timesHeaviside(j) == 0
        bj = b0;
    else
        bj = b_snapshot_time(2,aMICHeaviside,k,r,timesHeaviside(j),b0);
    end

    plot(ax,r,bj,'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_density_axes_31(ax,rMax);
legend(ax,time_labels_31(timeFractions), ...
    'Interpreter','latex','Location','southeast','FontSize',10,'Box','on');

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(figA,fullfile(outputDir,'Figure_3_1a.pdf'), ...
    'ContentType','vector','Resolution',300);
exportgraphics(figB,fullfile(outputDir,'Figure_3_1b.pdf'), ...
    'ContentType','vector','Resolution',300);

%% Figure 3.1c and Figure 3.1d: profiles for different death rates
%% Parameters
rMax = 10;
rZ = 4;
Nr = 800;
r = linspace(0,rMax,Nr);
r(1) = 1e-6;
b0 = double(r >= 1);

kValues = [1 10 50 100 1000];
timeFactor = 2;

%% Stopping times and threshold concentrations
aDirac = @(x,t) (1./(4*pi*t)).*exp(-(x.^2)./(4*t));
tZDirac = rZ^2/4;

g = @(y) rZ*besseli(1,rZ./(2*exp(y)),1) ...
    - besseli(0,rZ./(2*exp(y)),1);

options = optimoptions('fsolve','Display','none', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);

tZHeaviside = exp(fsolve(g,log(tZDirac),options));

aMICDirac = aDirac(rZ,tZDirac);

aMICHeaviside = (1/pi)*(1-marcumq( ...
    rZ/sqrt(2*tZHeaviside),1/sqrt(2*tZHeaviside)));

tFixedDirac = timeFactor*tZDirac;
tFixedHeaviside = timeFactor*tZHeaviside;

fprintf('Dirac:     t_Z = %.9f, a_MIC = %.9e\n', ...
    tZDirac,aMICDirac);

fprintf('Heaviside: t_Z = %.9f, a_MIC = %.9e\n', ...
    tZHeaviside,aMICHeaviside);

%% Line styles
styles = { ...
    struct('Color',[0.00 0.45 0.74],'LineStyle','-'), ...
    struct('Color',[0.85 0.33 0.10],'LineStyle','--'), ...
    struct('Color',[0.47 0.67 0.19],'LineStyle','-.'), ...
    struct('Color',[0.49 0.18 0.56],'LineStyle',':'), ...
    struct('Color',[0.30 0.30 0.30],'LineStyle','-')};

% Plot black curve first so that the green and purple curves are drawn
% above it wherever the curves overlap.
%
% Original k-order:
%   1, 10, 50, 100, 1000
%
% Plotting order:
%   1000, 1, 10, 50, 100
plotOrder = [5 1 2 3 4];

%% Figure 3.1c: Dirac kernel
figDirac = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 11 9]);

ax = axes(figDirac);
hold(ax,'on');

% Store handles according to the ORIGINAL legend order
hDirac = gobjects(numel(kValues),1);

for j = plotOrder
    bj = b_snapshot_rates( ...
        1,aMICDirac,kValues(j),r,tFixedDirac,b0);

    hDirac(j) = plot(ax,r,bj, ...
        'LineWidth',1.8, ...
        'Color',styles{j}.Color, ...
        'LineStyle',styles{j}.LineStyle);
end

format_density_axes_31_rates(ax,rMax);

legend(ax,hDirac,k_labels_31(kValues), ...
    'Interpreter','latex', ...
    'Location','southeast', ...
    'FontSize',9, ...
    'Box','on');

%% Figure 3.1d: Heaviside kernel
figHeaviside = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[14 2 11 9]);

ax = axes(figHeaviside);
hold(ax,'on');

% Again store handles according to the ORIGINAL legend order
hHeaviside = gobjects(numel(kValues),1);

for j = plotOrder
    bj = b_snapshot_rates( ...
        2,aMICHeaviside,kValues(j),r,tFixedHeaviside,b0);

    hHeaviside(j) = plot(ax,r,bj, ...
        'LineWidth',1.8, ...
        'Color',styles{j}.Color, ...
        'LineStyle',styles{j}.LineStyle);
end

format_density_axes_31_rates(ax,rMax);

legend(ax,hHeaviside,k_labels_31(kValues), ...
    'Interpreter','latex', ...
    'Location','southeast', ...
    'FontSize',9, ...
    'Box','on');

%% Export
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics( ...
    figDirac, ...
    fullfile(outputDir,'Figure_3_1c.pdf'), ...
    'ContentType','vector', ...
    'Resolution',300);

exportgraphics( ...
    figHeaviside, ...
    fullfile(outputDir,'Figure_3_1d.pdf'), ...
    'ContentType','vector', ...
    'Resolution',300);

%% Local functions for Figure 3.1a and Figure 3.1b

function bj = b_snapshot_time(kernelIndex,a_min,k,x,tj,b0)
% Numerical evaluation of the accumulated antimicrobial exposure.
% kernelIndex = 1: radial Dirac; kernelIndex = 2: radial Heaviside.

    t0 = 1e-6;

    if tj <= 0
        bj = b0;
        return;
    end

    % First identify the time range where the threshold is exceeded.
    tc = linspace(t0,tj,100);
    Ac = a_eval_time(kernelIndex,x,tc);

    if all(Ac <= a_min,'all')
        bj = b0;
        return;
    end

    tEnd = tc(find(any(Ac > a_min,1),1,'last'));

    if kernelIndex == 2
        Nq = 160;
    else
        Nq = 300;
    end

    tq = linspace(t0,tEnd,Nq);
    A = a_eval_time(kernelIndex,x,tq);

    exposure = trapz(tq,A.*(A > a_min),2);
    bj = (b0(:).*exp(-k*exposure)).';
end

function A = a_eval_time(kernelIndex,x,t)
% Antimicrobial concentration on the tensor product of x and t.

    X = x(:);
    T = t(:).';

    switch kernelIndex
        case 1
            A = (1./(4*pi*T)).*exp(-(X.^2)./(4*T));

        case 2
            % Keep this Marcum-Q call in the same vectorised form as the
            % original script. Replacing it with repeated scalar calls is
            % substantially slower for Figure 3.1.
            A = (1/pi)*(1-marcumq( ...
                X./sqrt(2*T),ones(numel(X),1)./sqrt(2*T)));

        otherwise
            error('Unknown kernelIndex.');
    end
end

function labels = time_labels_31(fractions)

labels = cell(size(fractions));
for j = 1:numel(fractions)
    labels{j} = sprintf('$t=%.0f\\%%\\,t_Z$',100*fractions(j));
end
end

function format_density_axes_31(ax,rMax)

box(ax,'on');
grid(ax,'on');
set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');
xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$b(r,t)$','Interpreter','latex');
xlim(ax,[0 rMax]);
ylim(ax,[-0.05 1.05]);
end

%% Local functions for Figure 3.1c and Figure 3.1d
function bj = b_snapshot_rates(kernelIndex,a_min,k,x,tj,b0)

    t0 = 1e-6;

    if tj <= 0
        bj = b0;
        return;
    end

    tc = linspace(t0,tj,80);
    Ac = a_eval_rates(kernelIndex,x,tc);

    if all(Ac <= a_min,'all')
        bj = b0;
        return;
    end

    tEnd = tc(find(any(Ac > a_min,1),1,'last'));

    if kernelIndex == 2
        Nq = 160;
    else
        Nq = 300;
    end

    tq = linspace(t0,tEnd,Nq);
    A = a_eval_rates(kernelIndex,x,tq);

    exposure = trapz(tq,A.*(A > a_min),2);

    bj = (b0(:).*exp(-k*exposure)).';
end

function A = a_eval_rates(kernelIndex,x,t)

    X = x(:);
    T = t(:).';

    switch kernelIndex

        case 1
            A = (1./(4*pi*T)).*exp(-(X.^2)./(4*T));

        case 2
            A = (1/pi)*(1-marcumq( ...
                X./sqrt(2*T), ...
                ones(numel(X),1)./sqrt(2*T)));

        otherwise
            error('Unknown kernelIndex.');
    end
end

function labels = k_labels_31(kValues)

    labels = arrayfun( ...
        @(x)sprintf('$k=%g$',x), ...
        kValues, ...
        'UniformOutput',false);
end

function format_density_axes_31_rates(ax,rMax)

    box(ax,'on');
    grid(ax,'on');

    set(ax, ...
        'FontName','Times', ...
        'FontSize',11, ...
        'LineWidth',0.9, ...
        'TickLabelInterpreter','latex');

    xlabel(ax,'$r$','Interpreter','latex');
    ylabel(ax,'$b(r,t)$','Interpreter','latex');

    xlim(ax,[0 rMax]);
    ylim(ax,[-0.05 1.05]);
end
