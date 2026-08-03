%% Figure 3.2: bacterial density profiles for different death rates
% Both panels are evaluated at 2t_Z. The left panel uses the Dirac delta
% initial condition and the right panel uses the Heaviside initial condition.
%
% This script is self-contained. Running it produces Figure_3_2.pdf in the
% outputs folder alongside this file.

clear; close all; clc;

%% Parameters

rMax = 10;
rZ = 4;
Nr = 800;
r = linspace(0,rMax,Nr);
r(1) = 1e-6;
b0 = double(r >= 1);

kValues = [1 10 50 100];
timeFactor = 2;

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

tFixedDirac = timeFactor*tZDirac;
tFixedHeaviside = timeFactor*tZHeaviside;

fprintf('Dirac:     t_Z = %.6f, a_MIC = %.6e\n',tZDirac,aMICDirac);
fprintf('Heaviside: t_Z = %.6f, a_MIC = %.6e\n', ...
    tZHeaviside,aMICHeaviside);

%% Plot

styles = { ...
    struct('Color',[0.00 0.45 0.74],'LineStyle','-'), ...
    struct('Color',[0.85 0.33 0.10],'LineStyle','--'), ...
    struct('Color',[0.47 0.67 0.19],'LineStyle','-.'), ...
    struct('Color',[0.49 0.18 0.56],'LineStyle',':')};

fig = figure('Color','w','Units','centimeters','Position',[2 2 22 9]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl);
hold(ax,'on');
for j = 1:numel(kValues)
    bj = b_snapshot(1,aMICDirac,kValues(j),r,tFixedDirac,b0);
    plot(ax,r,bj,'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_density_axes(ax,rMax);
legend(ax,k_labels(kValues),'Interpreter','latex', ...
    'Location','southeast','FontSize',10,'Box','on');

ax = nexttile(tl);
hold(ax,'on');
for j = 1:numel(kValues)
    bj = b_snapshot(2,aMICHeaviside,kValues(j),r,tFixedHeaviside,b0);
    plot(ax,r,bj,'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_density_axes(ax,rMax);
legend(ax,k_labels(kValues),'Interpreter','latex', ...
    'Location','southeast','FontSize',10,'Box','on');

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(fig,fullfile(outputDir,'Figure_3_2.pdf'), ...
    'ContentType','vector','Resolution',300);

%% Local functions

function bj = b_snapshot(kernelIndex,a_min,k,x,tj,b0)
% Numerical evaluation of the accumulated antimicrobial exposure.

    t0 = 1e-6;

    if tj <= 0
        bj = b0;
        return;
    end

    tc = linspace(t0,tj,80);
    Ac = a_eval(kernelIndex,x,tc);

    if all(Ac <= a_min,'all')
        bj = b0;
        return;
    end

    tEnd = tc(find(any(Ac > a_min,1),1,'last'));
    if kernelIndex == 2
        Nq = 120;
    else
        Nq = 250;
    end

    tq = linspace(t0,tEnd,Nq);
    A = a_eval(kernelIndex,x,tq);

    exposure = trapz(tq,A.*(A > a_min),2);
    bj = (b0(:).*exp(-k*exposure)).';
end

function A = a_eval(kernelIndex,x,t)

    X = x(:);
    T = t(:).';

    switch kernelIndex
        case 1
            A = (1./(4*pi*T)).*exp(-(X.^2)./(4*T));

        case 2
            A = (1/pi)*(1-marcumq( ...
                X./sqrt(2*T),ones(numel(X),1)./sqrt(2*T)));

        otherwise
            error('Unknown kernelIndex.');
    end
end

function labels = k_labels(kValues)

labels = cell(size(kValues));
for j = 1:numel(kValues)
    labels{j} = sprintf('$k=%g$',kValues(j));
end
end

function format_density_axes(ax,rMax)

box(ax,'on');
grid(ax,'on');
set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');
xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$b(r,t)$','Interpreter','latex');
xlim(ax,[0 rMax]);
ylim(ax,[-0.05 1.05]);
end
