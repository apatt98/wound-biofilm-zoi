%% Figure 3.3: measurement-error dependence on k and beta
% This script examines how the error in the measured zone-of-inhibition
% radius depends jointly on the bacterial death-rate parameter k and the
% bacterial-density threshold beta.
%
% The measured boundary is defined by b(r,t) = beta and is compared with
% the ideal inhibition boundary obtained in the limit beta = 1.
%
% Figure 3.3a shows the relative measurement error at the final observation
% time across the (k,beta) parameter space.
%
% Figure 3.3b shows the maximum relative measurement error attained over
% the full observation interval.
%
% Running this script produces:
%   outputs/Figure_3_3a.pdf
%   outputs/Figure_3_3b.pdf

clear; close all; clc;

%% ------------------------------------------------------------------------
%  Model parameters
%  ------------------------------------------------------------------------

rZ = 4;

% From r_Z = 1/sqrt(e*pi*a_MIC)
aMIC = 1/(exp(1)*pi*rZ^2);

% From t_Z = r_Z^2/4
tZ = rZ^2/4;

% End time
tEnd = 3*tZ;

%% ------------------------------------------------------------------------
%  Parameter ranges
%  ------------------------------------------------------------------------

% These grids are more than sufficient for the displayed heatmap
% resolution while being substantially cheaper than a 1000 x 999 grid.

Nk = 350;
Nbeta = 350;

kVals = logspace(0,3,Nk);
betaVals = linspace(0.001,0.999,Nbeta);

%% ------------------------------------------------------------------------
%  Numerical grids for exposure calculation
%  ------------------------------------------------------------------------

Nr = 1000;
Nt = 700;

r = linspace(1,rZ,Nr);

% Use a mildly non-uniform time grid to retain greater resolution at early
% times, where the measurement error changes most rapidly.
s = linspace(0,1,Nt);
t = tEnd*s.^1.5;

R = r(:);

%% ------------------------------------------------------------------------
%  Antimicrobial concentration
%  ------------------------------------------------------------------------

A = zeros(Nr,Nt);

T = t(2:end);

A(:,2:end) = ...
    (1./(4*pi*T)) .* exp(-(R.^2)./(4*T));

%% ------------------------------------------------------------------------
%  Accumulated antimicrobial exposure
%  ------------------------------------------------------------------------

% Bacterial death occurs only while a >= a_MIC.

Akill = A;
Akill(A < aMIC) = 0;

% I(r,t) = integral_0^t a(r,s) H(a-a_MIC) ds
I = cumtrapz(t,Akill,2);

%% ------------------------------------------------------------------------
%  Ideal beta = 1 inhibition boundary
%  ------------------------------------------------------------------------

r1 = ones(1,Nt);

indGrowing = (t > 0) & (t < tZ);

qFront = 4*pi*aMIC*t(indGrowing);

rFront = sqrt(-4*t(indGrowing).*log(qFront));

r1(indGrowing) = max(1,min(rFront,rZ));

r1(t >= tZ) = rZ;

%% ------------------------------------------------------------------------
%  Collapse (k,beta) dependence to q = log(1/beta)/k
%  ------------------------------------------------------------------------

qMap = log(1./betaVals)./kVals(:);

qMin = min(qMap(:));
qMax = max(qMap(:));

% Use logarithmic spacing because q spans several orders of magnitude.
Nq = 1400;
qVals = logspace(log10(qMin),log10(qMax),Nq);

%% ------------------------------------------------------------------------
%  Calculate error as a function of q
%  ------------------------------------------------------------------------

errorEndQ = zeros(1,Nq);
errorMaxQ = zeros(1,Nq);

fprintf('Calculating measurement-error curves...\n');

for n = 2:Nt

    exposureProfile = I(:,n);

    % No killing has occurred anywhere outside the application region.
    if exposureProfile(1) <= 0
        continue;
    end

    % Exposure decreases with increasing r. Reverse the vectors so that
    % exposure is increasing for interp1.
    exposureIncreasing = flipud(exposureProfile);
    radiusIncreasingExposure = flipud(R);

    % Numerical integration creates regions with identical exposure
    % (especially I = 0). interp1 requires unique x values.
    [exposureUnique,indUnique] = unique( ...
        exposureIncreasing,'stable');

    radiusUnique = radiusIncreasingExposure(indUnique);

    rBetaQ = ones(1,Nq);

    if numel(exposureUnique) > 1

        Imin = exposureUnique(1);
        Imax = exposureUnique(end);

        % If q is smaller than the minimum resolved exposure, beta is
        % extremely close to one and r_beta lies at the outer boundary.
        indOuter = qVals <= Imin;

        rBetaQ(indOuter) = r1(n);

        % Values between the minimum and maximum exposure have an
        % interior r_beta determined by interpolation.
        indInterior = ...
            (qVals > Imin) & ...
            (qVals < Imax);

        if any(indInterior)

            rBetaQ(indInterior) = interp1( ...
                exposureUnique, ...
                radiusUnique, ...
                qVals(indInterior), ...
                'linear');

        end

        % For q >= Imax, even the bacteria at r = 1 have not fallen
        % sufficiently for the beta threshold to be attained.
        rBetaQ(qVals >= Imax) = 1;

    end

    % Restrict to the physically admissible interval.
    rBetaQ = max(1,min(rBetaQ,r1(n)));

    % Relative error in percent.
    errorNow = 100*(r1(n)-rBetaQ)./r1(n);

    errorNow = max(0,errorNow);

    % Update maximum error.
    errorMaxQ = max(errorMaxQ,errorNow);

    % Retain final-time error.
    if n == Nt
        errorEndQ = errorNow;
    end

end

fprintf('Exposure calculation complete.\n');

%% ------------------------------------------------------------------------
%  Map q-results back into (k,beta) space
%  ------------------------------------------------------------------------

% Interpolate in log(q), which is more accurate because q spans several
% orders of magnitude.

logQ = log(qVals);
logQMap = log(qMap);

errorEnd = interp1( ...
    logQ,errorEndQ,logQMap,'linear');

errorMax = interp1( ...
    logQ,errorMaxQ,logQMap,'linear');

%% ------------------------------------------------------------------------
%  Common plotting properties
%  ------------------------------------------------------------------------

maxColour = max([errorEnd(:);errorMax(:)]);

contourLevels = [1 5 10 25 50];

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% ------------------------------------------------------------------------
%  Figure 3.3a: error at t = 3 t_Z
%  ------------------------------------------------------------------------

figA = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 10.5 8.5]);

ax = axes(figA);
hold(ax,'on');

surf(ax,betaVals,kVals,errorEnd, ...
    'EdgeColor','none');

view(ax,2);

set(ax, ...
    'YScale','log', ...
    'FontName','Times', ...
    'FontSize',11, ...
    'LineWidth',0.9, ...
    'Layer','top', ...
    'TickLabelInterpreter','latex');

xlabel(ax,'$\beta$', ...
    'Interpreter','latex');

ylabel(ax,'$k$', ...
    'Interpreter','latex');

title(ax,'$E(t=3t_Z)$', ...
    'Interpreter','latex', ...
    'FontSize',12);

xlim(ax,[min(betaVals) max(betaVals)]);
ylim(ax,[min(kVals) max(kVals)]);

clim(ax,[0 maxColour]);

box(ax,'on');
grid(ax,'on');

contour(ax,betaVals,kVals,errorEnd, ...
    contourLevels, ...
    'k', ...
    'LineWidth',0.65, ...
    'ShowText','on');

cb = colorbar(ax);

cb.Label.String = 'Relative ZOI error (\%)';
cb.Label.Interpreter = 'latex';

set(cb, ...
    'FontName','Times', ...
    'FontSize',10, ...
    'TickLabelInterpreter','latex');

colormap(ax,turbo);

% Raster export is preferable for a dense heatmap. At 600 dpi this is
% substantially higher resolution than required at the final LaTeX size.
exportgraphics(figA, ...
    fullfile(outputDir,'Figure_3_3a.pdf'), ...
    'ContentType','image', ...
    'Resolution',600);

%% ------------------------------------------------------------------------
%  Figure 3.3b: maximum error over 0 <= t <= 3 t_Z
%  ------------------------------------------------------------------------

figB = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 10.5 8.5]);

ax = axes(figB);
hold(ax,'on');

surf(ax,betaVals,kVals,errorMax, ...
    'EdgeColor','none');

view(ax,2);

set(ax, ...
    'YScale','log', ...
    'FontName','Times', ...
    'FontSize',11, ...
    'LineWidth',0.9, ...
    'Layer','top', ...
    'TickLabelInterpreter','latex');

xlabel(ax,'$\beta$', ...
    'Interpreter','latex');

ylabel(ax,'$k$', ...
    'Interpreter','latex');

title(ax,'$\displaystyle\max_{0\leq t\leq3t_Z}E(t)$', ...
    'Interpreter','latex', ...
    'FontSize',12);

xlim(ax,[min(betaVals) max(betaVals)]);
ylim(ax,[min(kVals) max(kVals)]);

clim(ax,[0 maxColour]);

box(ax,'on');
grid(ax,'on');

contour(ax,betaVals,kVals,errorMax, ...
    contourLevels, ...
    'k', ...
    'LineWidth',0.65, ...
    'ShowText','on');

cb = colorbar(ax);

cb.Label.String = 'Relative ZOI error (\%)';
cb.Label.Interpreter = 'latex';

set(cb, ...
    'FontName','Times', ...
    'FontSize',10, ...
    'TickLabelInterpreter','latex');

colormap(ax,turbo);

exportgraphics(figB, ...
    fullfile(outputDir,'Figure_3_3b.pdf'), ...
    'ContentType','image', ...
    'Resolution',600);