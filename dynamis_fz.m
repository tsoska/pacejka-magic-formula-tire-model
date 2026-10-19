%% Dynamis: Fx–kappa data maps per 50 N (RR) — Pure MF fit (Massimo)
clc; clear; close all;

%% -------------------- Φόρτωμα δεδομένων --------------------
% === Pacejka params από σημειώσεις ===
use_manual_P = true;   % αν το κάνεις false θα ξαναψάξει .mat

P_manual = [ ...
   -1.43878473984818,  2.16894381997033,  3.28400630195105e-10, ...
   -0.378358308003816, -0.313035477022556, 0.208175235232244, ...
    30.9795890258771,  0.000118526283949617, -0.192430442880484, ...
   -0.669821923670443, -0.348918405484979, -1.97844080465441, ...
    3.51718376067039,  9.48646938200230,  0.0173339783986732 ];

load('DynamisPRC_Acc (1).mat');    
S = run01;

t  = S.VxEst(:,1);                
V  = S.VxEst(:,2);                
ax = S.Ax_DAIMURRAW(:,2);         
V  = fillmissing(V,'previous');
ax = fillmissing(ax,'previous');
if numel(ax) > 5, ax = movmedian(ax,5); end

%% -------------------- Παράμετροι οχήματος / aero --------------------
p.m    = 230;    p.g = 9.81;
p.wb   = 1.59;   p.x_cg = 0.779;  p.a = p.x_cg; p.b = p.wb - p.a;
p.h_cg = 0.275;  p.Rw = 0.21;
rho = 1.225;
p.rho = rho;
p.ClA = 5.3; 
p.CdA = 1.8;
p.x_cp_down = 0.731;   
p.h_cp_drag = 0;       % θεωρώ μηδέν από κέντρο βάρους για να μήν έχω ροπή λόγω drag

%% -------------------- Υπολογισμός Fz --------------------
[FzF, FzR, FzFL, FzFR, FzRL, FzRR, Fdown, Fdrag] = compute_Fz_timeseries(V, ax, p);

%% -------------------- Slip ratio από RPM --------------------
wFL = S.nMotorFL(:,2) * 2*pi/60;
wFR = S.nMotorFR(:,2) * 2*pi/60;
wRL = S.nMotorRL(:,2) * 2*pi/60;
wRR = S.nMotorRR(:,2) * 2*pi/60;
Vden = max(V, 1.0);                % ασφαλής παρονομαστής
kFL = (wFL*p.Rw - V) ./ Vden;
kFR = (wFR*p.Rw - V) ./ Vden;
kRL = (wRL*p.Rw - V) ./ Vden;
kRR = (wRR*p.Rw - V) ./ Vden;

%% -------------------- Fx_tot & επιμερισμός στο Fz --------------------
Crr   = 0.010;
q     = 0.5*p.rho.*V.^2;
Fdrag = q*p.CdA;
Fdown = q*p.ClA;
Fx_tot = p.m.*ax + Fdrag + Crr*(p.m*p.g + Fdown);

Fz_sum = FzFL + FzFR + FzRL + FzRR;  epsN = 1e-9;
w_FL = FzFL ./ max(Fz_sum, epsN);
w_FR = FzFR ./ max(Fz_sum, epsN);
w_RL = FzRL ./ max(Fz_sum, epsN);
w_RR = FzRR ./ max(Fz_sum, epsN);

Fx_FL = Fx_tot .* w_FL; 
Fx_FR = Fx_tot .* w_FR;
Fx_RL = Fx_tot .* w_RL;
Fx_RR = Fx_tot .* w_RR;

%% -------------------- time plots (Fz & Aero) --------------------
figure('Color','w');
plot(t,FzFL, t,FzFR, t,FzRL, t,FzRR, 'LineWidth',1.2); grid on
xlabel('t [s]'); ylabel('F_z per wheel [N]');
legend('FL','FR','RL','RR','Location','best');
title('Normal loads per wheel');

figure('Color','w');
plot(t, FzF, t, FzR, 'LineWidth',1.2); grid on
xlabel('t [s]'); ylabel('F_z axle sum [N]');
legend('Front axle sum','Rear axle sum','Location','best');
title('Axle normal loads');

figure('Color','w');
plot(t, Fdown, t, Fdrag, 'LineWidth',1.2); grid on
xlabel('t [s]'); ylabel('Force [N]');
legend('Downforce total','Drag total','Location','best');
title('Aero forces');

%% Ταχύτητα vs Χρόνος (από raw data)
figure('Color','w');
plot(t, V, 'k', 'LineWidth',1.5); grid on
xlabel('t [s]'); ylabel('V_x [m/s]');
title('Vehicle speed V_x vs time');

%% Επιτάχυνση vs Χρόνος
figure('Color','w');
plot(t, ax, 'b', 'LineWidth',1.5); grid on
xlabel('t [s]'); ylabel('a_x [m/s^2]');
title('Longitudinal acceleration vs time');

%% Slip ratio vs Χρόνος (τροχοί & axle averages)
% axle averages (για context)
kFront = ((wFL + wFR)/2*p.Rw - V) ./ max(((wFL + wFR)/2*p.Rw), 1.0);
kRear  = ((wRL + wRR)/2*p.Rw - V) ./ max(((wRL + wRR)/2*p.Rw), 1.0);

figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; % Rear wheels
plot(t, kRL, 'b', 'LineWidth',1.2); hold on;
plot(t, kRR, 'r--', 'LineWidth',1.2);
plot(t, kRear, 'k-', 'LineWidth',1.4);
yline(0,'k:');
grid on; xlabel('t [s]'); ylabel('\kappa');
title('Slip ratio vs time — Rear');
legend('RL','RR','Rear axle avg','Location','best');
ylim([-0.5 0.5]);

nexttile; % Front wheels
plot(t, kFL, 'g', 'LineWidth',1.2); hold on;
plot(t, kFR, 'm--', 'LineWidth',1.2);
plot(t, kFront, 'k-', 'LineWidth',1.4);
yline(0,'k:');
grid on; xlabel('t [s]'); ylabel('\kappa');
title('Slip ratio vs time — Front');
legend('FL','FR','Front axle avg','Location','best');
ylim([-0.5 0.5]);

%% -------------------- Fx–kappa maps ανά 50 N (±25) για RR, 4/σελίδα --------------------
Vmin_keep   = 15/3.6;
dt          = max(median(diff(t)), eps);
jerk        = [0; diff(ax)/dt];     % m/s^3
jerk_max    = 4;
Mbase       = (V >= Vmin_keep) & abs(jerk) <= jerk_max;

% RR μόνο
Fz_i = FzRR;  k_i = kRR;  Fx_i = Fx_RR;

% ΜΟΝΟ ΘΕΤΙΚΟ SLIP & ΜΟΝΟ ΘΕΤΙΚΕΣ ΔΥΝΑΜΕΙΣ
kPosMin = 0;                 
PosMask = (k_i >= kPosMin) & (Fx_i >= 0);

% bins 400..Fz_top κάθε 50 N (πλάτος ±25)
Fz_max   = max(Fz_i);      if isempty(Fz_max) || ~isfinite(Fz_max), Fz_max = 900; end
Fz_top   = 25*ceil(Fz_max/25);
Fz_all   = 400:50:Fz_top;   hw = 25;

% --- Κράτα μόνο όσα bins έχουν Ν>0 σημεία ---
keep_centers = [];
for cFz = Fz_all
    inBin = Mbase & PosMask & ...
            (Fz_i >= cFz - hw) & (Fz_i <  cFz + hw) & ...
            isfinite(k_i) & isfinite(Fx_i);
    if any(inBin)
        keep_centers(end+1) = cFz; 
    end
end
    errors_summary = [];  % [Fz, RMSE, MAE, MAXE]
if isempty(keep_centers)
    warning('RR: Δεν υπάρχουν bins με δεδομένα.');
else
    % 4 subplots ανά figure (2x2)
    pageSize = 4;
    nBins    = numel(keep_centers);
    nPages   = ceil(nBins / pageSize);

    for pg = 1:nPages
        i1 = (pg-1)*pageSize + 1;
        i2 = min(pg*pageSize, nBins);
        nThis = i2 - i1 + 1;

        % δυναμική διάταξη ώστε να μην υπάρχουν κενά tiles
        nCol = min(2, nThis);
        nRow = ceil(nThis / nCol);

        figure('Color','w','Name',sprintf('RR | Fx-\\kappa maps (bins %d–%d)', i1, i2));
        tl = tiledlayout(nRow, nCol,'TileSpacing','compact','Padding','compact');
        title(tl, sprintf('RR | Fx vs \\kappa per F_z bin (±%d N) | page %d/%d', hw, pg, nPages));

        for j = 1:nThis
            cFz = keep_centers(i1 + j - 1);

            inBin = Mbase & PosMask & ...
                    (Fz_i >= cFz - hw) & (Fz_i <  cFz + hw) & ...
                    isfinite(k_i) & isfinite(Fx_i);

            k_plot  = k_i(inBin);
            Fx_plot = Fx_i(inBin);

            nexttile;
            scatter(k_plot, Fx_plot, 8, 'filled'); grid on
            xlabel('\kappa'); ylabel('F_x [N]');
            title(sprintf('F_z \\approx %d\\pm%d N (N=%d)', cFz, hw, numel(Fx_plot)));
            xlim([max(0,kPosMin) 0.30]);
            ylim([0, max(Fx_plot)*1.05]);
        end
    end

%% -------- Export SR–FX pairs (RR) --------
targets=[700, 850, 900, 1050];   % N (bin centers)
ratios_BCDE = [];   % [Fz  B_ratio  C_ratio  D_ratio  E_ratio]

for cFz = targets
    inBin = Mbase & PosMask & ...
            (Fz_i >= cFz - hw) & (Fz_i < cFz + hw) & ...
            isfinite(k_i) & isfinite(Fx_i);

    SR = double(k_i(inBin));     
    FX = double(Fx_i(inBin));    

    [SR, ord] = sort(SR);  FX = FX(ord);

    if isempty(SR)
        warning('RR @ Fz≈%d±%d N: no points to export.', cFz, hw);
    else
        T = table(SR, FX);
        meta.cFz = cFz; meta.hw = hw; meta.N = height(T);
    end
end

%% -------- Pure Magic Formula (Massimo) (RR) --------
tolFz = hw;  % ±25 N
err_rows = [];                 % θα έχει γραμμές [Fz, RMSE, MAE]
err_mode = strings(0,1);       % 'scaled','original','missing','no-overlap','mf-failed','too-few'

for cFz = targets
    % ---- επιλογή σημείων για το bin ----
    mask = Mbase & PosMask & ...
           (Fz_i >= cFz - tolFz) & (Fz_i < cFz + tolFz) & ...
           isfinite(k_i) & isfinite(Fx_i);

    SR = double(k_i(mask));
    FX = double(Fx_i(mask));
    FZ = double(Fz_i(mask));

    if numel(SR) < 6
        fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (too few points for MF fit)\n', cFz);
        continue;
    end

    % ===== OUTER TRY: MF fit (Dynamis) =====
    try
        [B,C,D,E,xm,ym] = mf_massimo_from_arrays( ...
            SR, FX, FZ, cFz, tolFz, ...
            sprintf('RR @ F_z\\approx %d\\pm%d N', cFz, tolFz), true);

        fprintf('RR @ Fz≈%d: B=%.3f, C=%.3f, D=%.1f, E=%.3f, xm=%.4f, ym=%.1f (N=%d)\n', ...
            cFz, B, C, D, E, xm, ym, numel(SR));

        % === Smooth MF curve (Dynamis) στο κοινό grid ===
        kq_pos  = linspace(0, 0.30, 800);
        MF_fun  = @(x) D .* sin( C .* atan( B.*x - E.*(B.*x - atan(B.*x)) ) );
        Fxq_pos = MF_fun(kq_pos);

        % για συμβατότητα με τα υπόλοιπα scripts (προαιρετικό save)
        SR_dyn_MF = kq_pos(:);
        FX_dyn_MF = Fxq_pos(:);
        save('dyn_MF_curve_Fz850_pm25.mat', 'SR_dyn_MF', 'FX_dyn_MF', 'cFz', 'tolFz');

        % ====== defaults για τις rig καμπύλες πριν το INNER TRY ======
        Fx_rig_MF62    = nan(size(kq_pos));
        Fx_rig_scaledP = [];

        % ===== INNER TRY: Rig MF62 + scaling =====
        try
            % --- P του rig ---
            if exist('use_manual_P','var') && use_manual_P
                P_opt = P_manual;
            else
                tmp = load('P_opt_rig.mat');
                P_opt = tmp.P_opt;
            end

            % --- συνθήκες αναφοράς ---
            Fz_ref = cFz;
            IA_ref = deg2rad(0);
            P_ref  = 10*6.89476e3;

            % --- αρχική rig καμπύλη στο ίδιο grid ---
            kq = kq_pos(:);
            Fx_rig_MF62 = MF62_FX_vector(P_opt, kq, Fz_ref, IA_ref, P_ref);

            % --- υπολογισμός B,C,D,E του rig από 15P ---
            ex  = 1.0e-3; Fz0 = 1112.0554070627252; p0 = 82737.090;
            dfz = (Fz_ref - Fz0) / Fz0;
            dp  = (P_ref  - p0)  / p0;
            IA_sin = sin(IA_ref);

            C_rig = abs(P_opt(1));
            D_rig = Fz_ref * (P_opt(2) - abs(P_opt(3))*dfz) * ...
                    (1 + P_opt(4)*dp + P_opt(5)*dp^2) * (1 - P_opt(6)*IA_sin^2);
            Kx_rig = Fz_ref * (P_opt(7) + P_opt(8)*dfz) * exp(P_opt(9)*dfz) * ...
                     (1 + P_opt(10)*dp + P_opt(11)*dp^2);
            B_rig = Kx_rig / (C_rig * D_rig + ex);
            E_rig = (P_opt(12) + P_opt(13)*dfz + P_opt(14)*dfz^2) * (1 - P_opt(15));

            % --- λόγοι Dyn/Rig και στόχοι ---
            ratio_B = B / B_rig;  ratio_C = C / C_rig;
            ratio_D = D / D_rig;  ratio_E = E / E_rig;
            ratios_BCDE = [ratios_BCDE; cFz, ratio_B, ratio_C, ratio_D, ratio_E];

            B_target = B_rig * ratio_B;
            C_target = C_rig * ratio_C;
            D_target = D_rig * ratio_D;
            E_target = E_rig * ratio_E;

            fprintf('Dyn/Rig @ Fz≈%d:  B=%.3f  C=%.3f  D=%.3f  E=%.3f\n', ...
                    cFz, ratio_B, ratio_C, ratio_D, ratio_E);

            % --- overlay plots ---
            figure('Color','w','Name', ...
                sprintf('FX–\\kappa overlay @ F_z \\approx %d \\pm %d N', cFz, tolFz));
            hold on; grid on;
            plot(kq_pos, Fxq_pos,      'b-', 'LineWidth', 2, 'DisplayName','Dynamis MF (B,C,D,E)');
            plot(kq_pos, Fx_rig_MF62,  'r-', 'LineWidth', 2, 'DisplayName','Rig MF62 (15P)');
            xlabel('Slip Ratio \kappa'); ylabel('F_x [N]');
            title(sprintf('FX–\\kappa overlay @ F_z \\approx %d \\pm %d N', cFz, tolFz));
            legend('Location','best');

            % --- fit P_scaled από B,C,D,E_target ---
            P_scaled = fit_MF62_from_BCDE(B_target, C_target, D_target, E_target, ...
                                          cFz, 0, 10*6.89476e3, P_opt(:));
            Fx_rig_scaledP = MF62_FX_vector(P_scaled, kq_pos, cFz, IA_ref, P_ref);

            figure('Color','w','Name', ...
                sprintf('FX–\\kappa overlay (NEW P) @ F_z \\approx %d \\pm %d N', cFz, tolFz));
            hold on; grid on;
            plot(kq_pos, Fxq_pos,       'b-', 'LineWidth', 2, 'DisplayName','Dynamis MF');
            plot(kq_pos, Fx_rig_MF62,   'r--','LineWidth', 1.6,'DisplayName','Rig MF62 (original P)');
            plot(kq_pos, Fx_rig_scaledP,'m-', 'LineWidth', 2, 'DisplayName','Rig MF62 (P from scaled B,C,D,E)');
            xlabel('Slip Ratio \kappa'); ylabel('F_x [N]');
            title(sprintf('FX–\\kappa (original vs NEW P) @ F_z \\approx %d \\pm %d N', cFz, tolFz));
            legend('Location','best');

           if ~isempty(Fx_rig_scaledP) && numel(Fx_rig_scaledP)==numel(Fxq_pos)
    valid = isfinite(Fxq_pos(:)) & isfinite(Fx_rig_scaledP(:));
    if any(valid)
        e    = Fxq_pos(valid) - Fx_rig_scaledP(valid);
        RMSE = sqrt(mean(e.^2,'omitnan'));
        MAE  = mean(abs(e),'omitnan');

        % --- Εκτύπωση στο Command Window ---
        fprintf('RMSE=%.2f | MAE=%.2f @ Fz≈%d N\n', RMSE, MAE, cFz);

        % --- ⬇️ ΝΕΟ: Καταγραφή στο table σφαλμάτων ---
        err_rows(end+1,:) = [cFz, RMSE, MAE];
        err_mode(end+1,1) = "scaled";

        % --- Σχολιασμός στο plot χωρίς range() ---
        xl = xlim; yl = ylim;
        dx = xl(2) - xl(1); dy = yl(2) - yl(1);
        txt = sprintf('RMSE=%.1f  MAE=%.1f N', RMSE, MAE);
        text(xl(1)+0.05*dx, yl(2)-0.12*dy, txt, ...
             'BackgroundColor',[1 1 1 0.7], 'EdgeColor',[0.7 0.7 0.7], ...
             'FontSize',10);
    else
        fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (no overlapping valid points)\n', cFz);
        err_rows(end+1,:) = [cFz, NaN, NaN];
        err_mode(end+1,1) = "no-overlap";
    end
else
    fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (scaled rig curve not available)\n', cFz);
    err_rows(end+1,:) = [cFz, NaN, NaN];
    err_mode(end+1,1) = "missing";
           end

        catch ME  % ===== end INNER TRY =====
            warning('Rig MF62 overlay failed @ Fz≈%d: %s', cFz, ME.message);

            % Fallback: σύγκριση με την αρχική rig καμπύλη
            if ~isempty(Fx_rig_MF62) && numel(Fx_rig_MF62)==numel(Fxq_pos)
                valid_fb = isfinite(Fxq_pos(:)) & isfinite(Fx_rig_MF62(:));
                if any(valid_fb)
                    e_fb = Fxq_pos(valid_fb) - Fx_rig_MF62(valid_fb);
                    RMSE = sqrt(mean(e_fb.^2,'omitnan'));
                    MAE  = mean(abs(e_fb),'omitnan');
                    fprintf('RMSE=%.2f | MAE=%.2f @ Fz≈%d N (vs ORIGINAL rig)\n', RMSE, MAE, cFz);
                    err_rows(end+1,:) = [cFz, RMSE, MAE];
err_mode(end+1,1) = "original";

                else
                    fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (no fallback points)\n', cFz);
                    err_rows(end+1,:) = [cFz, NaN, NaN];
err_mode(end+1,1) = "missing";

                end
            else
                fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (original rig curve not available)\n', cFz);
                err_rows(end+1,:) = [cFz, NaN, NaN];
err_mode(end+1,1) = "missing";

            end
        end  % <-- END inner try/catch

    catch ME  % ===== end OUTER TRY =====
        warning('MF (Massimo) fit failed at RR @ %d±%d: %s', cFz, tolFz, ME.message);
        fprintf('RMSE=NaN | MAE=NaN @ Fz≈%d N (MF fit failed)\n', cFz);
    end  % <-- END outer try/catch

end  % <-- END for cFz loop

% ---- ΕΚΤΥΠΩΣΗ συνοπτικών λόγων Dyn/Rig ΜΕΤΑ το loop (προαιρετικό) ----
if exist('ratios_BCDE','var') && ~isempty(ratios_BCDE)
    fprintf('\n===== Ratios Dynamis / Rig =====\n');
    disp(array2table(ratios_BCDE, ...
        'VariableNames', {'Fz','B_ratio','C_ratio','D_ratio','E_ratio'}));
    mean_B = mean(ratios_BCDE(:,2),'omitnan');
    mean_C = mean(ratios_BCDE(:,3),'omitnan');
    mean_D = mean(ratios_BCDE(:,4),'omitnan');
    mean_E = mean(ratios_BCDE(:,5),'omitnan');
    fprintf('\n===== Mean Ratios (Dynamis / Rig) =====\n');
    fprintf('B_mean = %.3f\nC_mean = %.3f\nD_mean = %.3f\nE_mean = %.3f\n', ...
        mean_B, mean_C, mean_D, mean_E);
end
end
% ====== Τελικός πίνακας σφαλμάτων για όλα τα Fz ======
if ~isempty(err_rows)
    % Τακτοποίηση/στρογγυλοποίηση για όμορφη εκτύπωση
    Fz_col   = err_rows(:,1);
    RMSE_col = round(err_rows(:,2), 2);
    MAE_col  = round(err_rows(:,3), 2);
    Mode_col = err_mode;  % string array

    % Προαιρετικά: ταξινόμηση κατά Fz
    [Fz_col, order] = sort(Fz_col);
    RMSE_col = RMSE_col(order);
    MAE_col  = MAE_col(order);
    Mode_col = Mode_col(order);

    errors_table = table(Fz_col, RMSE_col, MAE_col, Mode_col, ...
        'VariableNames', {'Fz_N','RMSE_N','MAE_N','ComparedTo'});

    disp('===== Error table (Dynamis vs Rig) =====');
    disp(errors_table);
else
    disp('===== Error table (Dynamis vs Rig) =====');
    disp('No error rows collected.');
end



%% ============================== ΒΟΗΘΗΤΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ ==============================
function [B, C, D, E, xm, ym] = mf_massimo_from_arrays(SR_in, FX_in, FZ_in, Fz_target, tolerance, plotTitle, doPlot)
% Απλή εκδοχή Massimo:
% B = (κλίση) από (0,0) στο πρώτο σημείο
% C από alpha, D & xm από peak (εδώ: μέσος των 5 υψηλότερων), E από τύπο.

    if nargin < 6 || isempty(plotTitle)
        plotTitle = 'MF fit';
    end
    if nargin < 7
        doPlot = true;      % default: να κάνει plot
    end

    % --- Φίλτρο ---
    M = isfinite(SR_in) & isfinite(FX_in) & isfinite(FZ_in) & ...
        (SR_in > 0) & (FX_in >= 0) & (abs(FZ_in - Fz_target) <= tolerance);
    SR = double(SR_in(M));
    FX = double(FX_in(M));

    if numel(SR) < 2
        B = 0; C = 1.46;
        [D,ix] = max(FX);
        if isempty(D), D=0; xm=0.08; else, xm = SR(ix); end
        ym = D; E = 0;
        return
    end

    % --- Sort & unique-ish ---
    [SR, ord] = sort(SR);  FX = FX(ord);
    keep = [true; diff(SR) > 1e-9];
    SR = SR(keep); FX = FX(keep);

    % --- D & xm από μέσο όρο των 5 υψηλότερων ---
    [FX_sorted, idx_sorted] = sort(FX, 'descend');
    nTop  = min(5, numel(FX_sorted));
    topIdx = idx_sorted(1:nTop);

    D  = mean(FX(topIdx));
    xm = mean(SR(topIdx));
    ym = D;

    % --- C presupposed ---
    alpha = 0.75;
    C = 2 - (2/pi)*asin(alpha);

    % --- B από τα πρώτα σημεία ---
    nB = min(1, numel(SR));
    xB = SR(1:nB); yB = FX(1:nB);
    B  = sum(yB) / ( (sum(xB)*C*D) + eps );

    % --- E από τον τύπο ---
    Bxm = B * xm;
    denom = (Bxm - atan(Bxm));
    if abs(denom) < 1e-12
        E = 0;
    else
        E = (Bxm - tan(pi/(2*C))) / denom;
    end
    E = min(max(E, -10), 1);

    % --- Plot (μόνο αν ζητηθεί) ---
    if doPlot
        kmax = max(SR)*1.1;
        kq = linspace(-kmax, kmax, 1000);
        MF = @(x) D .* sin( C .* atan( B.*x - E.*(B.*x - atan(B.*x)) ) );
        Fxq = MF(kq);

        assignin('base','MF_kq',kq);
        assignin('base','MF_Fxq',Fxq);
        assignin('base','MF_Fz',Fz_target);

        figure('Color','w','Name',plotTitle);
        hold on; grid on;
        scatter(SR, FX, 15, 'filled', 'DisplayName','measured (+ only)');
        plot(kq, Fxq, 'r-', 'LineWidth', 2, 'DisplayName','Magic Formula (sym.)');
        xlabel('Slip Ratio \kappa'); ylabel('F_x [N]');
        title(sprintf('%s | B=%.3f, C=%.3f, D=%.0f, E=%.3f, x_m=%.3f (N=%d)', ...
            plotTitle, B, C, D, E, xm, numel(SR)));
        legend('Location','best');
        xlim([-kmax, kmax]);
        ylim([-D*1.2, D*1.2]);
    end
end

function [FzF, FzR, FzFL, FzFR, FzRL, FzRR, Fdown, Fdrag] = compute_Fz_timeseries(V, ax, p)
% V, ax: διανύσματα στήλης ίδιου μήκους
m  = p.m;   g  = p.g;   wb = p.wb;
a  = p.a;   b  = p.b;   h  = p.h_cg;

rho = p.rho; ClA = p.ClA; CdA = p.CdA;
xcpd = p.x_cp_down;  hcpd = p.h_cp_drag;

% Στατικό
Fz0F = m*g*(b/wb);
Fz0R = m*g*(a/wb);

% Aero
q      = 0.5*rho.*V.^2;
Fdown  = q*ClA;
Fdrag  = q*CdA;

% Κατανομή downforce εμπρός/πίσω
FdownF = Fdown.*(wb - xcpd)/wb;
FdownR = Fdown.*(xcpd)/wb;

% Ροπή από drag (αν έχει offset)
Mdrag  = Fdrag.*hcpd;
dFdrag = Mdrag / wb;
FdragF = -dFdrag;  FdragR = +dFdrag;

% Διαμήκης μεταφορά βάρους
dFlong = m.*ax.*h / wb;
FxfF   = -dFlong;
FxfR   = +dFlong;

% Αθροίσματα αξόνων
FzF = Fz0F + FdownF + FdragF + FxfF;
FzR = Fz0R + FdownR + FdragR + FxfR;

% L/R ίσα (ay=0)
FzFL = 0.5*FzF; FzFR = 0.5*FzF;
FzRL = 0.5*FzR; FzRR = 0.5*FzR;

% clip ασφαλείας
FzF  = max(FzF,  0);  FzR  = max(FzR,  0);
FzFL = max(FzFL, 0);  FzFR = max(FzFR, 0);
FzRL = max(FzRL, 0);  FzRR = max(FzRR, 0);
end
function Fx = MF62_FX_vector(P, SL_vec, Fz_vec, IA_vec, P_vec)
    ex = 1.0e-3;
    Fz0 = 1112.0554070627252;
    p0  = 82737.090;

    IA  = sin(IA_vec);
    dfz = (Fz_vec - Fz0) ./ Fz0;
    dp  = (P_vec - p0)  ./ p0;

    C = abs(P(1));
    D = Fz_vec .* (P(2) - abs(P(3)) .* dfz) .* (1 + P(4)*dp + P(5)*dp.^2) .* (1 - P(6)*IA.^2);
    Kx = Fz_vec .* (P(7) + P(8)*dfz) .* exp(P(9)*dfz) .* (1 + P(10)*dp + P(11)*dp.^2);
    B  = Kx ./ (C .* D + ex);
    E  = (P(12) + P(13)*dfz + P(14)*dfz.^2) .* (1 - P(15) * sign(SL_vec));
    Fx = D .* sin( C .* atan( B .* SL_vec - E .* (B .* SL_vec - atan(B .* SL_vec)) ) );
end
function [B_mf, C_mf, D_mf, E_mf] = MF62_BCDE_from_P(P, Fz_ref, IA_ref, P_ref)
    % ίδιοι ορισμοί με αυτούς που έχεις πιο πάνω
    ex  = 1.0e-3;
    Fz0 = 1112.0554070627252;
    p0  = 82737.090;

    dfz    = (Fz_ref - Fz0) / Fz0;
    dp     = (P_ref  - p0)  / p0;
    IA_sin = sin(IA_ref);

    % C
    C_mf = abs(P(1));

    % D
    D_mf = Fz_ref * (P(2) - abs(P(3)) * dfz) * ...
           (1 + P(4)*dp + P(5)*dp^2) * ...
           (1 - P(6)*IA_sin^2);

    % Kx
    Kx = Fz_ref * (P(7) + P(8)*dfz) * exp(P(9)*dfz) * ...
         (1 + P(10)*dp + P(11)*dp^2);

    % B
    B_mf = Kx / (C_mf * D_mf + ex);

    % E (για + slip)
    E_mf = (P(12) + P(13)*dfz + P(14)*dfz^2) * (1 - P(15));
end
function J = MF62_BCDE_cost(P, B_t, C_t, D_t, E_t, Fz_ref, IA_ref, P_ref, P_anchor)
    % P_anchor = το "παλιό / base" P (π.χ. P_opt ή P_manual)
    [B_mf, C_mf, D_mf, E_mf] = MF62_BCDE_from_P(P, Fz_ref, IA_ref, P_ref);

    % σφάλμα σχετικό (για να μη βαραίνει μόνο το D)
    eB = (B_mf - B_t) / max(abs(B_t), 1e-6);
    eC = (C_mf - C_t) / max(abs(C_t), 1e-6);
    eD = (D_mf - D_t) / max(abs(D_t), 1e-3);
    eE = (E_mf - E_t) / max(abs(E_t), 1e-6);

    wB = 1; wC = 1; wD = 2; wE = 1;   % βάρη, παίξε αν θες

    err_term = wB*eB^2 + wC*eC^2 + wD*eD^2 + wE*eE^2;

    % regularization: να μην απομακρυνθεί πολύ από το αρχικό P
    if nargin < 9 || isempty(P_anchor)
        reg_term = 0;
    else
        reg_term = 1e-2 * sum( (P - P_anchor).^2 );  % 1e-2 δοκιμαστικό
    end

    J = err_term + reg_term;
end
function P_new = fit_MF62_from_BCDE(B_t, C_t, D_t, E_t, ...
                                    Fz_ref, IA_ref, P_ref, ...
                                    P_start)
    % Αν δεν δώσεις αρχικό, βάλε μηδενικά ή κάποιο χειροκίνητο
    if nargin < 8 || isempty(P_start)
        P_start = zeros(15,1);
    end

    lb = -10 * ones(15,1);
    ub =  10 * ones(15,1);

    cost = @(P) MF62_BCDE_cost(P, B_t, C_t, D_t, E_t, ...
                               Fz_ref, IA_ref, P_ref, P_start);

    opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','interior-point', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 5e4);

    P_new = fmincon(cost, P_start, [], [], [], [], lb, ub, [], opts);
end

























