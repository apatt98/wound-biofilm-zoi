%% Figure 2.1: antimicrobial concentration profiles and MIC level curve
% This script illustrates the diffusion of an antimicrobial from a Dirac
% delta initial condition on an infinite radial domain.
%
% Figure 2.1a shows the antimicrobial concentration profile a(r,t) at
% selected times relative to the stopping time t_Z. The minimum inhibitory
% concentration a_MIC and the maximum inhibition radius r_Z are also shown.
%
% Figure 2.1b shows the corresponding level curve a(r,t) = a_MIC in the
% (r,t) plane, with the times used in Figure 2.1a marked on the curve.
%
% Running this script produces:
%   outputs/Figure_2_1a.pdf
%   outputs/Figure_2_1b.pdf

clear; close all; clc;

%% Parameters

aMIC = 0.02;
rMax = 3;
Nr = 1000;
r = linspace(0,rMax,Nr);
r(1) = 1e-6;

% Closed-form maximum radius and stopping time.
rZ = 1/sqrt(exp(1)*pi*aMIC);
tZ = 1/(4*exp(1)*pi*aMIC);

timeFractions = [0, 0.25, 0.5, 1, 1.5, 2];
times = timeFractions*tZ;

aDirac = @(rr,tt) (1./(4*pi*tt)).*exp(-(rr.^2)./(4*tt));

curveColour = [0 0 1];
lineWidth = 1.8;

%% Figure 2.1a: concentration profiles

figA = figure('Color','w','Units','centimeters','Position',[2 2 22 12]);
tl = tiledlayout(figA,2,3,'TileSpacing','compact','Padding','compact');

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

        % Radius at which a(r,t) = a_MIC at this time.
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

%% Figure 2.1b: infinite-domain Dirac level curve
% This is the same Case 1 level curve used in Figure 2.2:
% A1(r,t) = (1/(4*pi*t))*exp(-r^2/(4*t)), with A1 = a_MIC.

% Use the same time range as Figure 2.2 and a fine grid for a smooth contour.
tMinLevel = 1e-3;
tMaxLevel = 5;
NrLevel = 700;
NtLevel = 700;

rLevelGrid = linspace(0,rMax,NrLevel);
tLevelGrid = linspace(tMinLevel,tMaxLevel,NtLevel);
[RLevel,TLevel] = meshgrid(rLevelGrid,tLevelGrid);
A1 = (1./(4*pi*TLevel)).*exp(-(RLevel.^2)./(4*TLevel));

% Exact radius corresponding to each Figure 2.1a time on a(r,t)=a_MIC.
% The t=0 point is the limiting point (r,t)=(0,0).
rAtTimes = NaN(size(times));
rAtTimes(1) = 0;
for j = 2:numel(times)
    tj = times(j);
    logArgument = 4*pi*tj*aMIC;
    if logArgument > 0 && logArgument <= 1
        rAtTimes(j) = sqrt(-4*tj*log(logArgument));
    end
end

figB = figure('Color','w','Units','centimeters','Position',[2 2 16 12]);
axB = axes(figB);
hold(axB,'on');
box(axB,'on');
grid(axB,'on');
set(axB,'FontName','Times','FontSize',12,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');

% Plot the infinite-domain Dirac MIC level curve used for Figure 2.1b.
contour(axB,RLevel,TLevel,A1,[aMIC aMIC], ...
    'LineWidth',2,'Color',curveColour,'LineStyle','-');

% Add the Figure 2.1a times at their corresponding radii.
valid = isfinite(rAtTimes);
plot(axB,rAtTimes(valid),times(valid),'o', ...
    'LineStyle','none','MarkerSize',7.5, ...
    'MarkerFaceColor',curveColour,'MarkerEdgeColor','k','LineWidth',1.0);

% Label each marked point with the same time fraction used in Figure 2.1a.
% Each label is offset approximately normal to the local level curve so that
% it stays close to its marker without crossing the curve.
labelDx = [0.05, -0.05, -0.05, 0.07, 0.06, 0.06];
labelDy = [0.06,  0.06,  0.06, 0.00, 0.06, 0.06];
labelHorizontal = {'left','right','right','left','left','left'};
labelVertical   = {'bottom','bottom','bottom','middle','bottom','bottom'};

for j = 1:numel(times)
    if ~isfinite(rAtTimes(j))
        continue;
    end

    if j == 1
        labelText = '$t=0$';
    else
        labelText = sprintf('$t=%.0f\\%%\\,t_Z$',100*timeFractions(j));
    end


    text(axB,rAtTimes(j)+labelDx(j),times(j)+labelDy(j),labelText, ...
        'Interpreter','latex','FontSize',10, ...
        'HorizontalAlignment',labelHorizontal{j}, ...
        'VerticalAlignment',labelVertical{j});
end



xlabel(axB,'$r$','Interpreter','latex');
ylabel(axB,'$t$','Interpreter','latex');
xlim(axB,[0 rMax]);
ylim(axB,[0 tMaxLevel]);

%% Report coordinates

fprintf('Figure 2.1: r_Z = %.6f, t_Z = %.6f\n',rZ,tZ);
fprintf('\nFigure 2.1b marked level-curve coordinates:\n');
for j = 1:numel(times)
    if isfinite(rAtTimes(j))
        fprintf('  %5.2f t_Z:  t = %.6f,  r = %.6f\n', ...
            timeFractions(j),times(j),rAtTimes(j));
    else
        fprintf('  %5.2f t_Z:  no real a_MIC level radius\n',timeFractions(j));
    end
end

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

exportgraphics(figA,fullfile(outputDir,'Figure_2_1a.pdf'), ...
    'ContentType','vector','Resolution',300);

exportgraphics(figB,fullfile(outputDir,'Figure_2_1b.pdf'), ...
    'ContentType','vector','Resolution',300);
