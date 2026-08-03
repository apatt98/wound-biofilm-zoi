%% Figure 2.1: antimicrobial concentration profiles
% Dirac delta initial condition on an infinite domain.
%
% This script is self-contained. Running it produces Figure_2_1.pdf in the
% outputs folder alongside this file.

clear; close all; clc;

%% Parameters

aMIC = 0.02;
rMax = 3;
Nr = 1000;
r = linspace(0, rMax, Nr);
r(1) = 1e-6;

% Closed-form maximum radius and stopping time.
rZ = 1/sqrt(exp(1)*pi*aMIC);
tZ = 1/(4*exp(1)*pi*aMIC);

timeFractions = [0, 0.25, 0.5, 1, 1.5, 2];
times = timeFractions*tZ;

aDirac = @(r,t) (1./(4*pi*t)).*exp(-(r.^2)./(4*t));

%% Plot

fig = figure('Color','w','Units','centimeters','Position',[2 2 22 12]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

curveColour = [0 0 1];
lineWidth = 1.8;
yLimits = [0 0.12];

for j = 1:numel(times)
    ax = nexttile(tl);
    hold(ax,'on');
    box(ax,'on');
    grid(ax,'on');

    set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
        'TickLabelInterpreter','latex');

    tj = times(j);

    if tj == 0
        % The Dirac initial condition is represented schematically at r = 0.
        plot(ax,[0 0],yLimits,'Color',curveColour, ...
            'LineWidth',lineWidth);
        title(ax,'$t=0$','Interpreter','latex');
    else
        aj = aDirac(r,tj);
        plot(ax,r,aj,'Color',curveColour,'LineWidth',lineWidth);
        title(ax,sprintf('$t=%.0f\\%%\\,t_Z$',100*timeFractions(j)), ...
            'Interpreter','latex');

        % Radius at which a(r,t) = a_MIC at the chosen time.
        logArgument = 4*pi*tj*aMIC;
        if logArgument > 0 && logArgument <= 1
            rLevel = sqrt(-4*tj*log(logArgument));
            if isreal(rLevel) && rLevel <= rMax
                plot(ax,rLevel,aMIC,'o','MarkerSize',7, ...
                    'MarkerFaceColor',curveColour, ...
                    'MarkerEdgeColor',curveColour);
            end
        end
    end

    yline(ax,aMIC,'k--','LineWidth',1.1,'HandleVisibility','off');
    xline(ax,rZ,'k--','LineWidth',1.1,'HandleVisibility','off');

    text(ax,0.5,0.03,'$a_{\mathrm{MIC}}$','Interpreter','latex', ...
        'FontSize',10,'BackgroundColor','w','Margin',1.5);
    text(ax,2.5,0.105,'$r_Z$','Interpreter','latex', ...
        'FontSize',10,'BackgroundColor','w','Margin',1.5);

    xlabel(ax,'$r$','Interpreter','latex');
    ylabel(ax,'$a(r,t)$','Interpreter','latex');
    xlim(ax,[0 rMax]);
    ylim(ax,yLimits);
end

fprintf('Figure 2.1: r_Z = %.6f, t_Z = %.6f\n',rZ,tZ);

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(fig,fullfile(outputDir,'Figure_2_1.pdf'), ...
    'ContentType','vector','Resolution',300);
