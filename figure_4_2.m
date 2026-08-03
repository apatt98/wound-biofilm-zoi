%% Figure 4.2: meropenem data and fitted model prediction
% This script contains the data, fitting, combined uncertainty sampling and
% plotting code needed for Figure 4.2.

clear; close all; clc;

% rng(10) for reproducibility
rng(10);

%% Parameters and data

cfg.drug_name = 'Meropenem';
cfg.dose_ug = 1000;
cfg.agar_depth_mm = 4;
cfg.D0_mm2_min = 2.0e-2;
cfg.D_sensitivity_range_mm2_min = [0.5 1.5]*cfg.D0_mm2_min;
cfg.Nsamp = 5000;
cfg.ribbon_lo_pct = 5;
cfg.ribbon_hi_pct = 95;
cfg.Ncurve = 600;
cfg.xlim_mm = [0 35];
cfg.effective_loading_factor_range = [0.5 1.5];
cfg.reference_loading_factor = mean(cfg.effective_loading_factor_range);
cfg.tobs_mode = 'discrete';
cfg.tobs_values_hr = [18 24];
cfg.tobs0_hr = mean(cfg.tobs_values_hr);

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
    0.25      29.0;
    0.0625  30.0;
    0.25    29.5;
    0.25    30.0;
    0.625   30.5;
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

mer_data.MIC = [mer_pseudomonas(:,1); mer_acinetobacter(:,1)];
mer_data.R = [mer_pseudomonas(:,2); mer_acinetobacter(:,2)];
mer_data.group = [repmat("Pseudomonas",size(mer_pseudomonas,1),1); ...
                  repmat("Acinetobacter",size(mer_acinetobacter,1),1)];

%% Fit the reference diffusivity

[D_fit,sse_fit] = fit_diffusivity_to_radius( ...
    mer_data.MIC,mer_data.R,ones(size(mer_data.R)), ...
    cfg.dose_ug,cfg.agar_depth_mm,cfg.reference_loading_factor, ...
    cfg.tobs0_hr,cfg.D_sensitivity_range_mm2_min);

cfg.D0_mm2_min = D_fit;
cfg.D_sensitivity_range_mm2_min = [0.5 1.5]*D_fit;

%% Model curve and uncertainty band

amin_grid = make_amin_grid(mer_data.MIC,cfg.Ncurve);
ygrid = log10(amin_grid);
yexp = log10(mer_data.MIC);
dose_ref_ug = cfg.dose_ug*cfg.reference_loading_factor;

[~,r_fix,r_stop,use_fix,use_stop,amincrit] = ...
    combined_curve_physical(amin_grid,dose_ref_ug,cfg.agar_depth_mm, ...
    cfg.D0_mm2_min,60*cfg.tobs0_hr);

R_all = sample_prediction_curves(cfg,amin_grid);
r_lo = prctile(R_all,cfg.ribbon_lo_pct,2);
r_hi = prctile(R_all,cfg.ribbon_hi_pct,2);

%% Plot

fig = figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
ax = axes(fig);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontName','Times','FontSize',12,'LineWidth',0.9, ...
    'TickLabelInterpreter','latex');

fill_band(ax,r_lo,r_hi,ygrid,cfg);

isP = mer_data.group == "Pseudomonas";
isA = mer_data.group == "Acinetobacter";
scatter(ax,mer_data.R(isP),yexp(isP),30,'o', ...
    'MarkerFaceColor',[0.82 0.82 0.82],'MarkerEdgeColor','k', ...
    'DisplayName','Pseudomonas');
scatter(ax,mer_data.R(isA),yexp(isA),34,'s', ...
    'MarkerFaceColor',[0.92 0.92 0.92],'MarkerEdgeColor','k', ...
    'DisplayName','Acinetobacter');

plot_branches(ax,r_fix,r_stop,use_fix,use_stop,ygrid);

joinRadius = sqrt(4*cfg.D0_mm2_min*60*cfg.tobs0_hr);
xline(ax,joinRadius,':','Join','Interpreter','latex', ...
    'HandleVisibility','off');

xlabel(ax,'ZOI radius (mm)','Interpreter','latex');
ylabel(ax,'$\log_{10}(MIC)~(\mu\mathrm{g/mL})$','Interpreter','latex');
xlim(ax,cfg.xlim_mm);
legend(ax,'Interpreter','latex','Location','southwest', ...
    'FontSize',10,'Box','on');

fprintf('Figure 4.2: fitted meropenem D_a* = %.7f mm^2/min\n',D_fit);
fprintf('Figure 4.2: SSE = %.6g\n',sse_fit);
fprintf('Figure 4.2: critical MIC = %.6g ug/mL\n',amincrit);

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
safe_export_pdf(fig,fullfile(outputDir,'Figure_4_2.pdf'));

%% Local functions

function R_all = sample_prediction_curves(cfg,amin_grid)
% Sample diffusivity, observation time and effective loading together.

    R_all = nan(numel(amin_grid),cfg.Nsamp);

    for j = 1:cfg.Nsamp
        D_j = 10.^unifrnd(log10(cfg.D_sensitivity_range_mm2_min(1)), ...
            log10(cfg.D_sensitivity_range_mm2_min(2)));

        if strcmp(cfg.tobs_mode,'discrete')
            tobs_hr = cfg.tobs_values_hr(randi(numel(cfg.tobs_values_hr)));
        elseif strcmp(cfg.tobs_mode,'uniform')
            tobs_hr = unifrnd(cfg.tobs_range_hr(1),cfg.tobs_range_hr(2));
        else
            error('Unknown tobs_mode: %s',cfg.tobs_mode);
        end

        loading_factor = unifrnd(cfg.effective_loading_factor_range(1), ...
            cfg.effective_loading_factor_range(2));
        dose_j = cfg.dose_ug*loading_factor;

        R_all(:,j) = combined_curve_physical(amin_grid,dose_j, ...
            cfg.agar_depth_mm,D_j,60*tobs_hr);
    end
end


function [D_fit,sse_fit] = fit_diffusivity_to_radius(MIC,R_obs,weights, ...
        dose_ug,agar_depth_mm,loading_factor,tobs_hr,D_bounds)
% Fit D on a logarithmic scale so that the diffusivity remains positive.

    if nargin < 3 || isempty(weights)
        weights = ones(size(MIC));
    end

    dose_eff_ug = dose_ug*loading_factor;
    tobs_min = 60*tobs_hr;

    objective = @(logD) weighted_sse_for_D(exp(logD),MIC,R_obs,weights, ...
        dose_eff_ug,agar_depth_mm,tobs_min);

    logD_fit = fminbnd(objective,log(D_bounds(1)),log(D_bounds(2)));
    D_fit = exp(logD_fit);
    sse_fit = objective(logD_fit);
end

function sse = weighted_sse_for_D(D,MIC,R_obs,weights, ...
        dose_eff_ug,agar_depth_mm,tobs_min)

    R_pred = combined_curve_physical(MIC,dose_eff_ug, ...
        agar_depth_mm,D,tobs_min);
    residuals = R_obs-R_pred;
    valid = isfinite(residuals) & isfinite(weights);
    sse = sum(weights(valid).*residuals(valid).^2);
end

function [r_comb,r_fix,r_stop,use_fix,use_stop,amincrit_ugmL] = ...
    combined_curve_physical(amin_ugmL,dose_ug,agar_depth_mm,D_mm2_min,tobs_min)
% Piecewise radius obtained from the dimensional Dirac solution.

    % 1 ug/mL = 1e-9 g/mm^3. Multiplication by agar depth gives an
    % effective surface concentration in g/mm^2.
    amin_surface = amin_ugmL*1e-9*agar_depth_mm;
    M_g = dose_ug*1e-6;

    r_stop = sqrt(M_g./(exp(1)*pi.*amin_surface));
    t_stop = M_g./(4*exp(1)*pi*D_mm2_min.*amin_surface);

    log_arg = 4*pi*D_mm2_min*tobs_min.*amin_surface./M_g;
    r_fix = nan(size(amin_ugmL));
    valid_fix = log_arg > 0 & log_arg < 1;
    r_fix(valid_fix) = sqrt(-4*D_mm2_min*tobs_min.*log(log_arg(valid_fix)));

    use_fix = t_stop > tobs_min & valid_fix;
    use_stop = ~use_fix;

    r_comb = nan(size(amin_ugmL));
    r_comb(use_fix) = r_fix(use_fix);
    r_comb(use_stop) = r_stop(use_stop);

    amincrit_surface = M_g/(4*exp(1)*pi*D_mm2_min*tobs_min);
    amincrit_ugmL = amincrit_surface/(1e-9*agar_depth_mm);
end

function amin_grid = make_amin_grid(MIC_values,Ncurve)

    amin_lo = max(min(MIC_values(MIC_values > 0))*0.8,1e-4);
    amin_hi = max(MIC_values)*1.2;
    amin_grid = logspace(log10(amin_lo),log10(amin_hi),Ncurve).';
end

function fill_band(ax,r_lo,r_hi,ygrid,cfg)

    valid = isfinite(r_lo) & isfinite(r_hi);
    fill(ax,[r_lo(valid); flipud(r_hi(valid))], ...
        [ygrid(valid); flipud(ygrid(valid))],[0.75 0.85 0.95], ...
        'FaceAlpha',0.45,'EdgeColor','none', ...
        'DisplayName',sprintf('%d--%d%% sensitivity band', ...
        cfg.ribbon_lo_pct,cfg.ribbon_hi_pct));
end

function plot_branches(ax,r_fix,r_stop,use_fix,use_stop,ygrid)

    plot(ax,r_stop(use_stop),ygrid(use_stop), ...
        'Color',[0 0.45 0.74],'LineWidth',2.1, ...
        'DisplayName','Stopping-radius branch');
    plot(ax,r_fix(use_fix),ygrid(use_fix), ...
        'Color',[0.85 0.33 0.10],'LineWidth',2.1, ...
        'DisplayName','Fixed-time branch');
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

