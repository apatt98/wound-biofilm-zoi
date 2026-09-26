%% Figure 3.2: error in the measured zone-of-inhibition radius
% This script quantifies the difference between the ideal antimicrobial
% inhibition boundary and the radius that would be inferred from a chosen
% bacterial-density threshold.
%
% The relative error is expressed as a percentage of the ideal
% zone-of-inhibition radius.
%
% Figure 3.2a shows the error over time for several values of the bacterial
% death-rate parameter k at a fixed bacterial-density threshold beta.
%
% Figure 3.2b shows the error over time for several values of beta at a
% fixed value of k.
%
% Running this script produces:
%   outputs/Figure_3_2a.pdf
%   outputs/Figure_3_2b.pdf

clear; close all; clc;

%% Parameters
rZ = 4;
aMIC = 1/(exp(1)*pi*rZ^2);
tZ = rZ^2/4;
timeFractions = linspace(0,3,800);
times = timeFractions*tZ;

betaFixed = 0.5;
kVals = [1 10 50 100 1000];
kFixed = 50;
betaVals = [0.25 0.50 0.75 0.90];

%% Ideal beta = 1 inhibition boundary
r1 = arrayfun(@(tj)r_beta(tj,1,1,aMIC,rZ,tZ),times);

%% Relative percentage error for varying k
errorK = zeros(numel(kVals),numel(times));
for i = 1:numel(kVals)
    for j = 1:numel(times)
        rb = r_beta(times(j),betaFixed,kVals(i),aMIC,rZ,tZ);
        errorK(i,j) = 100*(r1(j)-rb)/r1(j);
    end
end

%% Relative percentage error for varying beta
errorBeta = zeros(numel(betaVals),numel(times));
for i = 1:numel(betaVals)
    for j = 1:numel(times)
        rb = r_beta(times(j),betaVals(i),kFixed,aMIC,rZ,tZ);
        errorBeta(i,j) = 100*(r1(j)-rb)/r1(j);
    end
end

for i = 1:numel(kVals)
    [peakError,idx] = max(errorK(i,:));
    fprintf(['k = %g: E(3t_Z) = %.6f%%, max E = %.6f%% ', ...
        'at t/t_Z = %.4f\n'],kVals(i),errorK(i,end), ...
        peakError,timeFractions(idx));
end

%% Plot styles
styles = { ...
    struct('Color',[0.00 0.45 0.74],'LineStyle','-'), ...
    struct('Color',[0.85 0.33 0.10],'LineStyle','--'), ...
    struct('Color',[0.47 0.67 0.19],'LineStyle','-.'), ...
    struct('Color',[0.49 0.18 0.56],'LineStyle',':'), ...
    struct('Color',[0.30 0.30 0.30],'LineStyle','-')};

errorMax = max([errorK(:);errorBeta(:)]);
yMax = max(5,1.05*errorMax);

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Figure 3.2a: varying k
figA = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 10.5 8.5],'Renderer','painters');
ax = axes(figA);
hold(ax,'on');
for j = 1:numel(kVals)
    plot(ax,timeFractions,errorK(j,:),'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_error_axes(ax,yMax);
legend(ax,k_labels(kVals),'Interpreter','latex', ...
    'Location','northeast','FontSize',9,'Box','on');
exportgraphics(figA,fullfile(outputDir,'Figure_3_2a.pdf'), ...
    'ContentType','vector');

%% Figure 3.2b: varying beta
figB = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 10.5 8.5],'Renderer','painters');
ax = axes(figB);
hold(ax,'on');
for j = 1:numel(betaVals)
    plot(ax,timeFractions,errorBeta(j,:),'LineWidth',1.8, ...
        'Color',styles{j}.Color,'LineStyle',styles{j}.LineStyle);
end
format_error_axes(ax,yMax);
legend(ax,beta_labels(betaVals),'Interpreter','latex', ...
    'Location','northeast','FontSize',10,'Box','on');
exportgraphics(figB,fullfile(outputDir,'Figure_3_2b.pdf'), ...
    'ContentType','vector');

%% Local functions
function rb = r_beta(t,beta,k,aMIC,rZ,tZ)
    tol = 1e-10;
    if t <= 0
        rb = 1;
        return;
    end

    if abs(beta-1) < tol
        if t < tZ
            q = 4*pi*aMIC*t;
            if q < 1
                rFront = sqrt(-4*t*log(q));
                rb = max(1,min(rFront,rZ));
            else
                rb = 1;
            end
        else
            rb = rZ;
        end
        return;
    end

    rLeft = 1+1e-8;
    bLeft = b_dirac(rLeft,t,k,aMIC,rZ);
    if bLeft >= beta
        rb = 1;
        return;
    end

    rRight = rZ;
    f = @(x)b_dirac(x,t,k,aMIC,rZ)-beta;
    fLeft = f(rLeft);
    fRight = f(rRight);
    if abs(fLeft) < tol
        rb = rLeft;
    elseif abs(fRight) < tol
        rb = rRight;
    elseif fLeft*fRight > 0
        error(['Unable to bracket r_beta at t = %.6g, k = %.6g, ', ...
            'beta = %.6g. f(1+) = %.6e, f(r_Z) = %.6e.'], ...
            t,k,beta,fLeft,fRight);
    else
        rb = fzero(f,[rLeft rRight]);
    end
end

function b = b_dirac(r,t,k,aMIC,rZ)
    if t <= 0 || r >= rZ
        b = 1;
        return;
    end

    z = -pi*aMIC*r^2;
    Wm1 = double(lambertw(-1,z));
    W0 = double(lambertw(0,z));
    tMinus = -r^2/(4*Wm1);
    tPlus = -r^2/(4*W0);
    Ei = @(x)-expint(-x);

    if t < tMinus
        b = 1;
    elseif t <= tPlus
        exposure = Ei(Wm1)-Ei(-r^2/(4*t));
        b = exp(-k*exposure/(4*pi));
    else
        exposure = Ei(Wm1)-Ei(W0);
        b = exp(-k*exposure/(4*pi));
    end
end

function labels = k_labels(kVals)
    labels = arrayfun(@(x)sprintf('$k=%g$',x),kVals, ...
        'UniformOutput',false);
end

function labels = beta_labels(betaVals)
    labels = arrayfun(@(x)sprintf('$\\beta=%.2f$',x),betaVals, ...
        'UniformOutput',false);
end

function format_error_axes(ax,yMax)
    box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
        'TickLabelInterpreter','latex');
    xlabel(ax,'$t/t_Z$','Interpreter','latex');
    ylabel(ax,'Relative ZOI error (\%)','Interpreter','latex');
    xlim(ax,[0 3]);
    ylim(ax,[0 yMax]);
end
