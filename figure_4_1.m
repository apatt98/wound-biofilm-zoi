%% Figure 4.1: meropenem one-at-a-time sensitivity analysis
% This script contains the data, fit, sensitivity calculation and plotting
% code needed for Figure 4.1. It does not depend on any other project file.

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

fprintf('Figure 4.1: fitted meropenem D_a* = %.7f mm^2/min\n',D_fit);
fprintf('Figure 4.1: SSE = %.6g\n',sse_fit);

%% Sensitivity panels

amin_grid = make_amin_grid(mer_data.MIC,cfg.Ncurve);
ygrid = log10(amin_grid);
yexp = log10(mer_data.MIC);
D0 = cfg.D0_mm2_min;
t0_min = 60*cfg.tobs0_hr;
sensitivity = make_sensitivity_definitions(cfg);

fig = figure('Color','w','Units','centimeters','Position',[2 2 27 9]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

for p = 1:numel(sensitivity)
    [R_all,r_fix,r_stop,use_fix,use_stop] = ...
        sample_one_parameter_curves(cfg,amin_grid,sensitivity(p),D0,t0_min);

    r_lo = prctile(R_all,cfg.ribbon_lo_pct,2);
    r_hi = prctile(R_all,cfg.ribbon_hi_pct,2);

    ax = nexttile(tl);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax,'FontName','Times','FontSize',11,'LineWidth',0.9, ...
        'TickLabelInterpreter','latex');

    fill_band(ax,r_lo,r_hi,ygrid,cfg);

    isP = mer_data.group == "Pseudomonas";
    isA = mer_data.group == "Acinetobacter";
    scatter(ax,mer_data.R(isP),yexp(isP),24,'o', ...
        'MarkerFaceColor',[0.82 0.82 0.82],'MarkerEdgeColor','k', ...
        'HandleVisibility','off');
    scatter(ax,mer_data.R(isA),yexp(isA),28,'s', ...
        'MarkerFaceColor',[0.92 0.92 0.92],'MarkerEdgeColor','k', ...
        'HandleVisibility','off');

    plot_branches(ax,r_fix,r_stop,use_fix,use_stop,ygrid);

    title(ax,['Varying ',sensitivity(p).name],'Interpreter','latex');
    xlabel(ax,'ZOI radius (mm)','Interpreter','latex');
    ylabel(ax,'$\log_{10}(MIC)~(\mu\mathrm{g/mL})$','Interpreter','latex');
    xlim(ax,cfg.xlim_mm);

    if p == 1
        legend(ax,'Interpreter','latex','Location','southwest', ...
            'FontSize',8,'Box','on');
    end
end

%% Export

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'outputs');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
safe_export_pdf(fig,fullfile(outputDir,'Figure_4_1.pdf'));

%% Local functions

function sensitivity = make_sensitivity_definitions(cfg)
% Define the one-at-a-time sampling used in the three panels.

    sensitivity(1).name = '$D_a^*$';
    sensitivity(1).type = 'D';
    sensitivity(1).sampling = 'loguniform';
    sensitivity(1).range = cfg.D_sensitivity_range_mm2_min;

    sensitivity(2).name = '$t_{\mathrm{obs}}^*$';
    sensitivity(2).type = 't';
    if strcmp(cfg.tobs_mode,'discrete')
        sensitivity(2).sampling = 'discrete';
        sensitivity(2).values = 60*cfg.tobs_values_hr;
    else
        sensitivity(2).sampling = 'uniform';
        sensitivity(2).range = 60*cfg.tobs_range_hr;
    end

    sensitivity(3).name = 'effective loading';
    sensitivity(3).type = 'dosefactor';
    sensitivity(3).sampling = 'uniform';
    sensitivity(3).range = cfg.effective_loading_factor_range;
end

function [R_all,r_fix_base,r_stop_base,use_fix_base,use_stop_base] = ...
    sample_one_parameter_curves(cfg,amin_grid,sensitivity,D0,t0_min)

    R_all = nan(numel(amin_grid),cfg.Nsamp);
    dose_ref_ug = cfg.dose_ug*cfg.reference_loading_factor;

    [~,r_fix_base,r_stop_base,use_fix_base,use_stop_base] = ...
        combined_curve_physical(amin_grid,dose_ref_ug, ...
        cfg.agar_depth_mm,D0,t0_min);

    for j = 1:cfg.Nsamp
        dose_j = dose_ref_ug;
        D_j = D0;
        t_j = t0_min;
        sampled_value = sample_sensitivity_value(sensitivity);

        switch sensitivity.type
            case 'D'
                D_j = sampled_value;
            case 't'
                t_j = sampled_value;
            case 'dosefactor'
                dose_j = cfg.dose_ug*sampled_value;
            otherwise
                error('Unknown sensitivity type: %s',sensitivity.type);
        end

        R_all(:,j) = combined_curve_physical(amin_grid,dose_j, ...
            cfg.agar_depth_mm,D_j,t_j);
    end
end

function value = sample_sensitivity_value(sensitivity)

    switch sensitivity.sampling
        case 'loguniform'
            value = 10.^unifrnd(log10(sensitivity.range(1)), ...
                log10(sensitivity.range(2)));
        case 'uniform'
            value = unifrnd(sensitivity.range(1),sensitivity.range(2));
        case 'discrete'
            value = sensitivity.values(randi(numel(sensitivity.values)));
        otherwise
            error('Unknown sampling mode: %s',sensitivity.sampling);
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

