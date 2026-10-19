clc;
clear;
close all;
addpath(genpath('./'));

% === Clear Output Directory ===
outRoot = 'Data Plots and Mats';
if isfolder(outRoot)
    try
        delete(fullfile(outRoot, '**', '*.*'));
        pause(0.1);
        rmdir(outRoot, 's');
    catch ME
        warning('Δεν ήταν δυνατή η διαγραφή του φακέλου: %s\n%s', outRoot, ME.message);
    end
end

%% === Step 1: Load and Combine Data ===
load('B2356run71.mat');
T1 = table(FZ, V, IA, P, FX, FY, SA, SL);
load('B2356run72.mat');
T2 = table(FZ, V, IA, P, FX, FY, SA, SL);
load('B2356run73.mat');
T3 = table(FZ, V, IA, P, FX, FY, SA, SL);
inputTable = [T1; T2; T3];

dataFolder = fullfile(outRoot, 'Data_Fx');
plotFolder = fullfile(outRoot, 'Plots_Fx');
mkdir(dataFolder);
mkdir(plotFolder);
save(fullfile(outRoot, 'Longitudinal_Data.mat'), 'inputTable');

%% === Step 2: Categorization ===
nRows = height(inputTable);
FZ_lb = inputTable.FZ * 0.2248089431;

catLoad = strings(nRows,1);
catP    = strings(nRows,1);
catIA   = strings(nRows,1);
catV    = strings(nRows,1);

for i = 1:nRows
    FZ_ = FZ_lb(i);
    P_  = inputTable.P(i);
    IA_ = inputTable.IA(i);
    V_  = inputTable.V(i);

    if  -55 < FZ_ && FZ_ <  -45, catLoad(i)="50";
    elseif -109< FZ_ && FZ_ < -91, catLoad(i)="100";
    elseif -160< FZ_ && FZ_ < -140,catLoad(i)="150";
    elseif -212< FZ_ && FZ_ < -188,catLoad(i)="200";
    elseif -265< FZ_ && FZ_ < -235,catLoad(i)="250";
    else, catLoad(i)="NA";
    end

    if   49 < P_ && P_ <  61, catP(i)="8";
    elseif 62 < P_ && P_ <  76, catP(i)="10";
    elseif 77 < P_ && P_ <  90, catP(i)="12";
    elseif 91 < P_ && P_ < 102, catP(i)="14";
    else, catP(i)="NA";
    end

    if   -0.7 < IA_ && IA_ < 0.7, catIA(i)="0";
    elseif 1.2 < IA_ && IA_ < 2.8, catIA(i)="2";
    elseif 3.2 < IA_ && IA_ < 4.8, catIA(i)="4";
    else, catIA(i)="NA";
    end

    if   19 < V_ && V_ < 30, catV(i)="15";
    elseif 34 < V_ && V_ < 48, catV(i)="25";
    elseif 62 < V_ && V_ < 79, catV(i)="45";
    else, catV(i)="NA";
    end
end

inputTable.catLoad = catLoad;
inputTable.catP    = catP;
inputTable.catIA   = catIA;
inputTable.catV    = catV;
save(fullfile(outRoot, 'Longitudinal_Data.mat'), 'inputTable');

%% === Step 3: Save Categorized Files and Plot FX vs SL (and SA distribution) ===
validRows = catLoad ~= "NA" & catP ~= "NA" & catIA ~= "NA" & catV ~= "NA";
filteredTable = inputTable(validRows, :);
combinations = unique(filteredTable(:, ["catLoad", "catP", "catIA", "catV"]));

h = waitbar(0, 'Generating plots...');
for i = 1:height(combinations)
    waitbar(i / height(combinations), h);
    c = combinations(i, :);
    subTable = filteredTable(...
        filteredTable.catLoad == c.catLoad & ...
        filteredTable.catP == c.catP & ...
        filteredTable.catIA == c.catIA & ...
        filteredTable.catV == c.catV, :);

    nameTag = sprintf('L%s-P%s-IA%s-V%s', ...
        c.catLoad{1}, c.catP{1}, c.catIA{1}, c.catV{1});

    save(fullfile(dataFolder, [nameTag '.mat']), 'subTable');

    % Plot FX vs SL
    figure('Visible','off');
    scatter(subTable.SL, subTable.FX, 10, 'filled');
    xlabel('Slip Ratio');
    ylabel('Longitudinal Force FX (N)');
    title(['FX vs SL: ' nameTag]);
    grid on;
    saveas(gcf, fullfile(plotFolder, [nameTag '.png']));
    close;

    % Optional: plot SA distribution for the group
    figure('Visible','off');
    histogram(subTable.SA, 20);
    xlabel('Slip Angle (deg)');
    ylabel('Count');
    title(['SA distribution: ' nameTag]);
    saveas(gcf, fullfile(plotFolder, [nameTag '_SAhist.png']));
    close;
end
close(h);

%% === Step 4: Interpolation for FX (και συλλογή SA) ===
L_indx = [50; 100; 150; 200; 250];
P_indx = [8; 10; 12; 14];
IA_indx = [0;2;4];
V_indx = [15; 25; 45];

allSL = [];
allFX = [];
allLoad = [];
allP = [];
allIA = [];
allSA = []; % συλλογή slip angle

for l = 1:length(L_indx)
    for p = 1:length(P_indx)
        for ia = 1:length(IA_indx)
            for v = 1:length(V_indx)
                nameTag = sprintf('L%d-P%d-IA%d-V%d', ...
                    L_indx(l), P_indx(p), IA_indx(ia), V_indx(v));
                matFilePath = fullfile(dataFolder, [nameTag '.mat']);

                if exist(matFilePath, 'file')
                    subTable = load(matFilePath).subTable;
                    subTable = rmmissing(subTable);
                    if height(subTable) < 10, continue; end

                    % Sort by SL and take average FX per unique SL
                    [SL_sorted, sortIdx] = sort(subTable.SL);
                    FX_sorted = subTable.FX(sortIdx);
                    SA_sorted = subTable.SA(sortIdx); % σταθμισμένο SA

                    % Make SL unique with mean FX and mean SA per SL bin
                    [SL_unique, ~, ic] = unique(SL_sorted);
                    FX_mean = accumarray(ic, FX_sorted, [], @mean);
                    SA_mean = accumarray(ic, SA_sorted, [], @mean);

                    if numel(SL_unique) < 2
                        continue;
                    end

                    % Interpolation grid on SL
                    SLi = linspace(min(SL_unique), max(SL_unique), 100)';
                    FXi = interp1(SL_unique, FX_mean, SLi, 'pchip');
                    SAi = interp1(SL_unique, SA_mean, SLi, 'pchip'); % interpolate SA vs SL

                    % Add to pooled dataset
                    allSL = [allSL; SLi];
                    allFX = [allFX; FXi];
                    allLoad = [allLoad; repmat(L_indx(l), length(SLi), 1)];
                    allP = [allP; repmat(P_indx(p), length(SLi), 1)];
                    allIA = [allIA; repmat(IA_indx(ia), length(SLi), 1)];
                    allSA = [allSA; SAi];
                end
            end
        end
    end
end

%% === Step 5: Optimization with GA + fmincon (combined model) ===
% Μετατροπές μονάδων / ρύθμιση εισόδων
Fz_N = allLoad / 0.2248089431;      % Fz σε N (lb->N)
Pinput = allP * 6.89476 * 1000;     % pressure σε Pa
IA_rad = deg2rad(allIA);            % camber σε rad
SL = allSL;                         % slip ratio
FX_target = allFX;                  % πραγματικά FX
SA_rad = deg2rad(allSA);            % slip angle σε rad

% --- Παράμετροι P ---
% 4 νέοι παράμετροι για το G_xα (C_xa, B_xa, E_xa, S_Hxa)
% -> συνολικά D = 19.
D = 19;
lb = -10 * ones(1, D);
ub = 10 * ones(1, D);

fitnessFcn = @(P) MF62_RMS_FX_Cost_combined(P, SL, Fz_N, IA_rad, Pinput, SA_rad, FX_target);

gaOpts = optimoptions('ga', ...
    'PopulationSize', 250, ...
    'MaxGenerations', 1000, ...
    'MaxStallGenerations', 150, ...
    'FunctionTolerance', 1e-8, ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'PlotFcn', {@gaplotbestf});

[P_ga, bestScore] = ga(fitnessFcn, D, [], [], [], [], lb, ub, [], gaOpts); %#ok<ASGLU>

fminOpts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'FunctionTolerance', 1e-6, ...
    'StepTolerance', 1e-8, ...
    'MaxIterations', 5000, ...
    'MaxFunctionEvaluations', 5e5);

problem = createOptimProblem('fmincon', ...
    'objective', fitnessFcn, ...
    'x0', P_ga, ...
    'options', fminOpts);

[P_opt, finalScore] = fmincon(problem); %#ok<ASGLU>

%% === Step 6: Validation Plot for a sample group ===
% ΓΝΩΣΤΕΣ pure-parameters (15)
PURE_FX_PARAMS = [ ...
   -1.43878473984818,  2.16894381997033,  3.28400630195105e-10, ...
   -0.378358308003816, -0.313035477022556, 0.208175235232244, ...
    30.9795890258771,  0.000118526283949617, -0.192430442880484, ...
   -0.669821923670443, -0.348918405484979, -1.97844080465441, ...
    3.51718376067039,  9.48646938200230,  0.0173339783986732 ];

% Φόρτωσε ένα subTable από το disk για να κάνεις σύγκριση
load(fullfile(dataFolder, 'L50-P8-IA2-V25.mat'));
Slip_Ratio = linspace(-0.3, 0.3, 1000);

% Μέση SA και IA του φακέλου για combined
meanSA = mean(subTable.SA, 'omitnan');
meanIA = deg2rad(mean(subTable.IA, 'omitnan'));
press  = 8*6.89476*1000;  % ενδεικτική πίεση
Fz_val = 222.411;         % ενδεικτικό Fz

% Combined prediction (βασισμένο στα 19 βελτιστοποιημένα)
Fx_comb_pred = MF62_FX_vector_combined(P_opt, Slip_Ratio, Fz_val, meanIA, press, deg2rad(meanSA));

% Pure prediction (με τις 15 γνωστές παραμέτρους)
Fx_pure_pred = MF62_FX_vector(PURE_FX_PARAMS, Slip_Ratio, Fz_val, meanIA, press);

figure('Color','w');
scatter(subTable.SL, subTable.FX, 15, 'b', 'filled'); hold on;
plot(Slip_Ratio, Fx_comb_pred, 'k',  'LineWidth', 2);
plot(Slip_Ratio, Fx_pure_pred, 'r--','LineWidth', 1.6);
xlabel('Slip Ratio');
ylabel('Longitudinal Force FX (N)');
legend('Measured FX', 'Predicted FX (combined)', 'Predicted FX (pure)','Location','best');
title('FX vs Slip Ratio: Real vs Predicted (combined vs pure)');
grid on;

%% === Functions: pure FX + combined + cost ===
function Fx = MF62_FX_vector_combined(P, SL_vec, Fz_vec, IA_vec, P_vec, SA_vec)
    % υπολογίζει πρώτα Fx_pure και μετά πολλαπλασιάζει με G_xa(SA).
    % P(1:15) -> pure longitudinal parameters
    % P(16:19) -> [C_xa, B_xa, E_xa, S_Hxa] παράμετροι για combined effect
    Ppure = P(1:15);
    Fx0 = MF62_FX_vector(Ppure, SL_vec, Fz_vec, IA_vec, P_vec);

    C_xa = abs(P(16));     % ensure C positive
    B_xa = abs(P(17));     % ensure B positive
    E_xa = P(18);
    S_Hxa = P(19);

    alpha_s = SA_vec + S_Hxa;  % SA σε rad
    Gxa = cos( C_xa .* atan( B_xa .* alpha_s - E_xa .* (B_xa .* alpha_s - atan(B_xa .* alpha_s)) ) );
    Gxa = max(Gxa, 0);         % μην πάει αρνητικό

    Fx = Gxa .* Fx0;
end

function Fx = MF62_FX_vector(P, SL_vec, Fz_vec, IA_vec, P_vec)
    % MF62 pure longitudinal
    ex = 0.0010;
    Fz0 = 1112.05540706272520;
    p0 = 82737.090;
    IA = sin(IA_vec);
    dfz = (Fz_vec - Fz0) ./ Fz0;
    dp = (P_vec - p0) ./ p0;

    C = abs(P(1));
    D = Fz_vec .* (P(2) - abs(P(3)) .* dfz) .* (1 + P(4)*dp + P(5)*dp.^2) .* (1 - P(6)*IA.^2);
    Kx = Fz_vec .* (P(7) + P(8)*dfz) .* exp(P(9)*dfz) .* (1 + P(10)*dp + P(11)*dp.^2);
    B = Kx ./ (C .* D + ex);
    E = (P(12) + P(13)*dfz + P(14)*dfz.^2) .* (1 - P(15) .* sign(SL_vec));

    Fx = D .* sin(C .* atan(B .* SL_vec - E .* (B .* SL_vec - atan(B .* SL_vec))));
end

function rmsError = MF62_RMS_FX_Cost_combined(P, SL, Fz, IA, Pinput, SA, realFx)
    Fx_pred = MF62_FX_vector_combined(P, SL, Fz, IA, Pinput, SA);
    rmsError = sqrt(mean((realFx - Fx_pred).^2));
end


