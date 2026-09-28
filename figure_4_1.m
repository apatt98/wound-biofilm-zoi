%% Figure 4.1: weighted least-squares fits to experimental data
% This script generates the parameter-estimation figures for ciprofloxacin
% and meropenem described in Section 4.3.
%
% Running this script produces:
%   outputs/Figure_4_1a.pdf
%   outputs/Figure_4_1b.pdf

clear; close all; clc;

%% Figure 4.1a: ciprofloxacin

rng(10);

%% Parameters and data

cfg.agar_depth_mm = 4;
cfg.tobs_hr = 16;
cfg.Nboot = 4000;
cfg.alpha = 0.05;
cfg.Ncurve = 600;
cfg.xlim_MIC = [0.01 0.5];
cfg.ylim_R2 = [150 520];

% Columns: [MIC (ug/mL), ZOI radius (mm), number of isolates]
cip_compact = [
0.016 16.5 1;
0.016 17.5 3;
0.016 18.0 7;
0.016 18.5 6;
0.016 19.0 6;
0.016 19.5 6;
0.016 20.0 6;
0.016 20.5 3;
0.016 21.0 4;
0.016 21.5 1;
0.016 22.0 1;

0.032 17.0 1;
0.032 17.5 1;
0.032 18.0 2;
0.032 18.5 6;
0.032 19.0 3;
0.032 20.0 5;
0.032 20.5 2;
0.032 21.0 2;

0.064 14.5 1;
0.064 15.0 1;
0.064 15.5 2;
0.064 16.5 2;
0.064 17.0 1;
0.064 17.5 2;
0.064 18.0 1;
0.064 18.5 1;
0.064 19.0 1;
0.064 20.0 1;

0.125 14.5 3;
0.125 15.5 8;
0.125 16.0 3;
0.125 16.5 2;
0.125 17.0 2;
0.125 18.5 1;

0.25 13.0 1;
0.25 13.5 2;
0.25 14.0 4;
0.25 14.5 9;
0.25 15.0 7;
0.25 15.5 4;
0.25 16.0 1;
0.25 16.5 1;

0.5 13.0 1;
0.5 14.5 1
];

cip_data.MIC = cip_compact(:,1);
cip_data.R = cip_compact(:,2);
cip_data.N = cip_compact(:,3);

%% Fit

fitC = fit_section_43( ...
    cip_data.MIC,cip_data.R,cip_data.N, ...
    cfg.agar_depth_mm,cfg.tobs_hr,cfg.Nboot,cfg.alpha);

MICgrid = logspace(log10(cfg.xlim_MIC(1)), ...
    log10(cfg.xlim_MIC(2)),cfg.Ncurve).';

[r2_fit,r2_lo,r2_hi] = mean_response_band( ...
    fitC,MICgrid,cfg.alpha);

%% Plot

figA = figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
ax = axes(figA);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');

set(ax,'FontName','Times','FontSize',12,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex','XScale','log');

fill(ax,[MICgrid; flipud(MICgrid)], ...
    [r2_lo; flipud(r2_hi)],[0.75 0.85 0.95], ...
    'FaceAlpha',0.45,'EdgeColor','none', ...
    'DisplayName','95\% confidence interval');

scatter(ax,fitC.MIC_center,fitC.R.^2,70,fitC.N,'s','filled', ...
    'MarkerEdgeColor','k','DisplayName','Salmonella isolates');

colourBar = colorbar(ax);
colourBar.Label.String = 'Number of isolates';
colourBar.Label.Interpreter = 'latex';
colourBar.TickLabelInterpreter = 'latex';

plot(ax,MICgrid,r2_fit, ...
    'Color',[0 0.45 0.74],'LineWidth',2.1, ...
    'DisplayName','Weighted least-squares fit');

xlabel(ax,'MIC ($\mu\mathrm{g/mL}$), dilution-bin centre', ...
    'Interpreter','latex');
ylabel(ax,'$r_{\mathrm{obs}}^{*2}~(\mathrm{mm}^2)$', ...
    'Interpreter','latex');

xlim(ax,cfg.xlim_MIC);
ylim(ax,cfg.ylim_R2);

text(ax,0.97,0.96, ...
    sprintf(['$\\Lambda = %.1f \\pm %.1f~\\mathrm{mm}^2$' newline ...
             '$R^2 = %.2f$' newline ...
             '$M^* = %.2f~\\mu\\mathrm{g}$' newline ...
             '95\\%% CI $[%.2f,\\ %.2f]~\\mu\\mathrm{g}$' newline ...
             'nominal $M^*=5~\\mu\\mathrm{g}$'], ...
    fitC.Lambda,fitC.SE_Lambda,fitC.R2,fitC.M_ug, ...
    fitC.M_CI_ug(1),fitC.M_CI_ug(2)), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',10,'Interpreter','latex');

legend(ax,'Interpreter','latex','Location','southwest', ...
    'FontSize',10,'Box','on');

fprintf('Figure 4.1a: Lambda = %.4f +/- %.4f mm^2\n', ...
    fitC.Lambda,fitC.SE_Lambda);
fprintf('Figure 4.1a: fitted D_a* = %.7f mm^2/min\n',fitC.D_mm2_min);
fprintf('Figure 4.1a: M* = %.6g ug\n',fitC.M_ug);
fprintf('Figure 4.1a: R^2 = %.6f\n',fitC.R2);

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

safe_export_pdf(figA,fullfile(outputDir,'Figure_4_1a.pdf'));

%% Figure 4.1b: meropenem

rng(10);

clear cfg MICgrid r2_fit r2_lo r2_hi ax colourBar;

%% Parameters and data

cfg.agar_depth_mm = 4;
cfg.Nboot = 4000;
cfg.alpha = 0.05;
cfg.Ncurve = 600;
cfg.xlim_MIC = [0.04 4];
cfg.ylim_R2 = [300 1250];

% Columns: [MIC (ug/mL), ZOI radius (mm)]
mer_pseudomonas = [
    2       29.0;
    2       24.5;
    2       24.5;
    2       25.5;
    2       24.5;
    2       24.5;
    2       25.0;
    2       27.5;
    2       24.0;
    2       25.0;
    0.0625  30.5;
    0.5     29.0;
    0.25    29.5;
    0.25    30.0;
    0.125   29.5;
    4       18.5;
    0.25    29.0;
    0.0625  30.0;
    0.25    29.5;
    0.25    30.0;
    0.0625  30.5;
    1       29.0;
    2       24.5;
    1       29.5;
    2       25.5;
    0.125   30.0;
    0.125   29.5;
    0.5     30.0;
    2       24.5;
    0.5     29.5;
    0.125   29.5;
    0.5     29.5;
    0.5     21.0
];

mer_acinetobacter = [
    2     18.0;
    2     20.0;
    2     19.0;
    0.5   20.0;
    0.5   20.5;
    0.5   20.0;
    1     22.0;
    1     21.5;
    1     20.5;
    2     18.5;
    1     21.0;
    2     19.0;
    1     21.0;
    2     19.0;
    1     21.0;
    1     22.0;
    1     21.5;
    1     21.0;
    1     22.0;
    2     18.0;
    0.5   22.0;
    1     21.0;
    0.5   22.5
];

%% Fit

fitP = fit_section_43( ...
    mer_pseudomonas(:,1),mer_pseudomonas(:,2), ...
    ones(size(mer_pseudomonas,1),1), ...
    cfg.agar_depth_mm,18,cfg.Nboot,cfg.alpha);

fitA = fit_section_43( ...
    mer_acinetobacter(:,1),mer_acinetobacter(:,2), ...
    ones(size(mer_acinetobacter,1),1), ...
    cfg.agar_depth_mm,24,cfg.Nboot,cfg.alpha);

MICgrid = logspace(log10(cfg.xlim_MIC(1)), ...
    log10(cfg.xlim_MIC(2)),cfg.Ncurve).';

[r2_P,r2_P_lo,r2_P_hi] = mean_response_band( ...
    fitP,MICgrid,cfg.alpha);

[r2_A,r2_A_lo,r2_A_hi] = mean_response_band( ...
    fitA,MICgrid,cfg.alpha);

%% Plot

figB = figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
ax = axes(figB);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');

set(ax,'FontName','Times','FontSize',12,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex','XScale','log');

fill(ax,[MICgrid; flipud(MICgrid)], ...
    [r2_P_lo; flipud(r2_P_hi)],[0.75 0.85 0.95], ...
    'FaceAlpha',0.35,'EdgeColor','none', ...
    'HandleVisibility','off');

fill(ax,[MICgrid; flipud(MICgrid)], ...
    [r2_A_lo; flipud(r2_A_hi)],[0.95 0.82 0.74], ...
    'FaceAlpha',0.35,'EdgeColor','none', ...
    'HandleVisibility','off');

scatter(ax,fitP.MIC_center,fitP.R.^2,30,'o', ...
    'MarkerFaceColor',[0.82 0.82 0.82],'MarkerEdgeColor','k', ...
    'DisplayName','Pseudomonas');

scatter(ax,fitA.MIC_center,fitA.R.^2,34,'s', ...
    'MarkerFaceColor',[0.92 0.92 0.92],'MarkerEdgeColor','k', ...
    'DisplayName','Acinetobacter');

plot(ax,MICgrid,r2_P, ...
    'Color',[0 0.45 0.74],'LineWidth',2.1, ...
    'DisplayName','Pseudomonas fit');

plot(ax,MICgrid,r2_A, ...
    'Color',[0.85 0.33 0.10],'LineWidth',2.1, ...
    'DisplayName','Acinetobacter fit');

xlabel(ax,'MIC ($\mu\mathrm{g/mL}$), dilution-bin centre', ...
    'Interpreter','latex');
ylabel(ax,'$r_{\mathrm{obs}}^{*2}~(\mathrm{mm}^2)$', ...
    'Interpreter','latex');

xlim(ax,cfg.xlim_MIC);
ylim(ax,cfg.ylim_R2);

text(ax,0.97,0.96, ...
    sprintf(['\\textit{Pseudomonas}: ' ...
             '$\\Lambda=%.1f\\pm%.1f~\\mathrm{mm}^2$, $R^2=%.2f$' newline ...
             '$M^*=%.0f~\\mu\\mathrm{g}$ (95\\%% CI %.0f--%.0f)' ...
             newline newline ...
             '\\textit{Acinetobacter}: ' ...
             '$\\Lambda=%.1f\\pm%.1f~\\mathrm{mm}^2$, $R^2=%.2f$' newline ...
             '$M^*=%.0f~\\mu\\mathrm{g}$ (95\\%% CI %.0f--%.0f)' ...
             newline newline ...
             'nominal $M^*=1000~\\mu\\mathrm{g}$'], ...
    fitP.Lambda,fitP.SE_Lambda,fitP.R2,fitP.M_ug, ...
    fitP.M_CI_ug(1),fitP.M_CI_ug(2), ...
    fitA.Lambda,fitA.SE_Lambda,fitA.R2,fitA.M_ug, ...
    fitA.M_CI_ug(1),fitA.M_CI_ug(2)), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',9.5,'Interpreter','latex');

legend(ax,'Interpreter','latex','Location','southwest', ...
    'FontSize',10,'Box','on');

fprintf('Figure 4.1b: Pseudomonas Lambda = %.4f +/- %.4f mm^2\n', ...
    fitP.Lambda,fitP.SE_Lambda);
fprintf('Figure 4.1b: Pseudomonas fitted D_a* = %.7f mm^2/min\n', ...
    fitP.D_mm2_min);
fprintf('Figure 4.1b: Pseudomonas M* = %.6g ug\n',fitP.M_ug);
fprintf('Figure 4.1b: Pseudomonas R^2 = %.6f\n',fitP.R2);

fprintf('Figure 4.1b: Acinetobacter Lambda = %.4f +/- %.4f mm^2\n', ...
    fitA.Lambda,fitA.SE_Lambda);
fprintf('Figure 4.1b: Acinetobacter fitted D_a* = %.7f mm^2/min\n', ...
    fitA.D_mm2_min);
fprintf('Figure 4.1b: Acinetobacter M* = %.6g ug\n',fitA.M_ug);
fprintf('Figure 4.1b: Acinetobacter R^2 = %.6f\n',fitA.R2);

%% Export

safe_export_pdf(figB,fullfile(outputDir,'Figure_4_1b.pdf'));

%% Local functions

function fit = fit_section_43(MIC,R,N,agar_depth_mm,tobs_hr,Nboot,alpha)
% Fit r_obs^2 against ln(MIC) using the procedure described in Section 4.3.

    MIC = MIC(:);
    R = R(:);
    N = N(:);

    MIC_center = MIC/sqrt(2);
    x = log(MIC_center);
    y = R.^2;
    weights = N./R.^2;

    X = [ones(size(x)) x];

    XtWX = X'*(weights.*X);
    beta = XtWX\(X'*(weights.*y));

    intercept = beta(1);
    Lambda = -beta(2);

    residuals = y-X*beta;
    n = sum(N);
    p = 2;

    weighted_sse = sum(weights.*residuals.^2);
    sigma2 = weighted_sse/(n-p);
    cov_beta = sigma2*inv(XtWX);

    SE_Lambda = sqrt(cov_beta(2,2));

    ybar = sum(weights.*y)/sum(weights);
    weighted_sst = sum(weights.*(y-ybar).^2);
    R2 = 1-weighted_sse/weighted_sst;

    M_g = pi*Lambda*agar_depth_mm*1e-9*exp(intercept/Lambda);
    M_ug = M_g*1e6;

    D_mm2_min = Lambda/(4*60*tobs_hr);

    expanded = expand_isolates(MIC,R,N);
    Mboot = nan(Nboot,1);

    for j = 1:Nboot
        ind = randi(size(expanded,1),size(expanded,1),1);
        sample = expanded(ind,:);
        Mboot(j) = fit_loading(sample(:,1),sample(:,2),agar_depth_mm);
    end

    M_CI_ug = percentile_local(Mboot, ...
        [100*alpha/2 100*(1-alpha/2)]);

    fit.MIC_center = MIC_center;
    fit.R = R;
    fit.N = N;
    fit.beta = beta;
    fit.cov_beta = cov_beta;
    fit.df = n-p;
    fit.Lambda = Lambda;
    fit.SE_Lambda = SE_Lambda;
    fit.D_mm2_min = D_mm2_min;
    fit.M_ug = M_ug;
    fit.M_CI_ug = M_CI_ug;
    fit.R2 = R2;
end

function [yfit,ylo,yhi] = mean_response_band(fit,MICgrid,alpha)

    Xgrid = [ones(numel(MICgrid),1) log(MICgrid(:))];

    yfit = Xgrid*fit.beta;

    variance = sum((Xgrid*fit.cov_beta).*Xgrid,2);
    se = sqrt(max(variance,0));

    if exist('tinv','file') == 2
        tcrit = tinv(1-alpha/2,fit.df);
    else
        tcrit = 1.96;
    end

    ylo = yfit-tcrit*se;
    yhi = yfit+tcrit*se;
end

function expanded = expand_isolates(MIC,R,N)

    expanded = nan(sum(N),2);
    row = 1;

    for j = 1:numel(MIC)
        expanded(row:row+N(j)-1,:) = repmat([MIC(j) R(j)],N(j),1);
        row = row+N(j);
    end
end

function M_ug = fit_loading(MIC,R,agar_depth_mm)

    x = log(MIC(:)/sqrt(2));
    y = R(:).^2;
    weights = 1./R(:).^2;

    X = [ones(size(x)) x];
    beta = (X'*(weights.*X))\(X'*(weights.*y));

    Lambda = -beta(2);

    if ~isfinite(Lambda) || Lambda <= 0
        M_ug = NaN;
        return;
    end

    intercept = beta(1);
    M_g = pi*Lambda*agar_depth_mm*1e-9*exp(intercept/Lambda);
    M_ug = M_g*1e6;
end

function q = percentile_local(x,p)

    x = sort(x(isfinite(x)));
    q = nan(size(p));

    for j = 1:numel(p)
        pos = 1+(numel(x)-1)*p(j)/100;
        lo = floor(pos);
        hi = ceil(pos);

        if lo == hi
            q(j) = x(lo);
        else
            q(j) = x(lo)+(pos-lo)*(x(hi)-x(lo));
        end
    end
end

function safe_export_pdf(fig,pdf_name)
% Transparent patches are exported as a rasterised PDF for reliability.

    try
        exportgraphics(fig,pdf_name,'ContentType','image','Resolution',300);
    catch ME
        warning('exportgraphics failed: %s',ME.message);
        set(fig,'PaperPositionMode','auto');
        print(fig,pdf_name,'-dpdf','-r300');
    end
end