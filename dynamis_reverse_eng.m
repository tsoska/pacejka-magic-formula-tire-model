%% Dynamis: Fx–kappa data maps per 50 N (RR) — Pure MF fit via reverse-engineering (GA/fmincon)
clc; clear; close all;

%% -------------------- Φόρτωμα δεδομένων --------------------
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
end  

%% -------- Pure Magic Formula via REVERSE-ENGINEERING (GA/fmincon) για 800 (RR) --------
%mfDir = fullfile(pwd, 'mf_opt_rr');
%if ~exist(mfDir, 'dir'), mkdir(mfDir); end

tolFz = hw;   % ±25 N
targets = 850;   

for cFz = targets
    mask = Mbase & PosMask & ...
           (Fz_i >= cFz - tolFz) & (Fz_i < cFz + tolFz) & ...
           isfinite(k_i) & isfinite(Fx_i);

    SR = double(k_i(mask));
    FX = double(Fx_i(mask));
    FZ = double(Fz_i(mask));

    if numel(SR) < 6
        warning('RR @ Fz≈%d±%d N: too few points for MF fit.', cFz, tolFz);
        continue;
    end

    try
        [B,C,D,E,xm,ym,score] = mf_optimize_from_arrays( ...
            SR, FX, FZ, cFz, tolFz, ...
            sprintf('RR @ F_z\\approx %d\\pm%d N (RE optimize)', cFz, tolFz));

        %save(fullfile(mfDir, sprintf('RR_MF_opt_Fz%d.mat', cFz)), ...
        %    'B','C','D','E','xm','ym','cFz','tolFz','score');

        fprintf(['RR @ Fz≈%d: B=%.3f, C=%.3f, D=%.1f, E=%.3f, ' ...
                 'xm=%.4f, ym=%.1f | RMSE=%.2f (N=%d)\n'], ...
                 cFz, B, C, D, E, xm, ym, score, numel(SR));

    catch ME
        warning('MF reverse-engineering failed at RR @ %d±%d: %s', ...
            cFz, tolFz, ME.message);
    end
end

%% ============================== ΒΟΗΘΗΤΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ ==============================
function [B, C, D, E, xm, ym, score] = mf_optimize_from_arrays(SR_in, FX_in, FZ_in, Fz_target, tolFz, plotTitle)
% Reverse-engineering Magic Formula παραμέτρων με GA -> fmincon.
% Κρατάμε C σταθερό (=1.46) και βελτιστοποιούμε μόνο B, D, E.
% MF(x) = D * sin( Cfix * atan( B*x - E*(B*x - atan(B*x)) ) )

    if nargin < 6, plotTitle = 'MF RE optimize'; end
    Cfix = 1.46;   % <<< ΣΤΑΘΕΡΟ C

    % ----------- ΦΙΛΤΡΟ & ΤΑΞΙΝΟΜΗΣΗ -----------
    M = isfinite(SR_in) & isfinite(FX_in) & isfinite(FZ_in) & ...
        (SR_in > 0) & (FX_in >= 0) & (abs(FZ_in - Fz_target) <= tolFz);
    SR = double(SR_in(M));
    FX = double(FX_in(M));

    if numel(SR) < 8
        error('Not enough points after filtering (need >= 8).');
    end

    [SR, ord] = sort(SR);
    FX = FX(ord);

    % Peak από ΔΕΔΟΜΕΝΑ:
    [ym, ipk] = max(FX);
    xm        = SR(ipk);

    % ----------- WEIGHTING γύρω από το peak -----------
    SR_peak = xm;
    sigma   = 0.05;
    w       = 1 + 2*exp(-((SR - SR_peak).^2)/(2*sigma^2));  % 1..3x

    % ----------- ΜΟΝΤΕΛΟ (C FIXED) -----------
    % P = [B, D, E]
    MF = @(P,x) P(2) .* sin( Cfix .* atan( P(1).*x - P(3).*(P(1).*x - atan(P(1).*x)) ) );

    % ----------- ΚΟΣΤΟΣ (weighted RMSE) -----------
    costFun = @(P) sqrt( mean( w .* ( MF(P,SR) - FX ).^2 ) );

    % ----------- ΟΡΙΑ / INIT -----------
    D0  = ym;
    lb  = [   0.01,  0.50*max(1,D0),  -2.0 ];
    ub  = [  80.00,  1.50*max(1,D0),   1.0 ];

    % αρχικοποίηση B από αρχική κλίση (χρησιμοποιεί Cfix)
    x1 = SR(1); y1 = FX(1);
    m0 = max(1e-6, y1/max(x1,1e-6));
    B0 = min( ub(1), max(lb(1), m0/(Cfix*max(D0,1)) ) );
    P0 = [B0, D0, 0.0];

    % ----------- GA -----------
    gaOpts = optimoptions('ga', ...
        'PopulationSize',          160, ...
        'MaxGenerations',          200, ...
        'MaxStallGenerations',     80,  ...
        'FunctionTolerance',       1e-7, ...
        'UseParallel',             false, ...
        'Display',                 'off');

    [P_ga, ~] = ga(costFun, 3, [], [], [], [], lb, ub, [], gaOpts);

    % ----------- fmincon refinement -----------
    fminOpts = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...
        'Display','off', ...
        'FunctionTolerance',1e-8, ...
        'StepTolerance',1e-10, ...
        'MaxIterations',2000, ...
        'MaxFunctionEvaluations',2e5);

    [P_opt, score] = fmincon(costFun, P_ga, [], [], [], [], lb, ub, [], fminOpts);

    B = P_opt(1);  C = Cfix;  D = P_opt(2);  E = P_opt(3);

    % ----------- Plot (με αντι-συμμετρική επέκταση) -----------
    kq_pos  = linspace(0, max(SR)*1.1, 800);
    Fxq_pos = MF(P_opt, kq_pos);
    kq_neg  = -fliplr(kq_pos);
    Fxq_neg = -fliplr(Fxq_pos);
    kq_full = [kq_neg, kq_pos];
    Fxq_full= [Fxq_neg, Fxq_pos];

    figure('Color','w','Name',plotTitle); hold on; grid on;
    scatter(SR, FX, 18, 'filled', 'DisplayName','measured (+ only)');
    plot(kq_full, Fxq_full, 'r-', 'LineWidth', 2, 'DisplayName','MF (sym. extended, C=1.46 fixed)');
    xlabel('Slip Ratio \kappa'); ylabel('F_x [N]');
    title(sprintf('F_z \\approx %d\\pm%d N (RE optimize) | B=%.3f, C=%.2f (fixed), D=%.1f, E=%.3f | RMSE=%.2f', ...
        Fz_target, tolFz, B, C, D, E, score));
    legend('Location','best');
end

function Fx = mf_eval(k, B, C, D, E)
    Bx = B.*k;
    Fx = D .* sin( C .* atan( Bx - E .* (Bx - atan(Bx)) ) );
end

function x = clamp(x, lb, ub)
    x = min(max(x, lb), ub);
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

