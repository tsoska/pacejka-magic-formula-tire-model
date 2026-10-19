% - Fitting στα RAW δεδομένα
% - Combined με G_yk = cos(Ck*atan(Byk*SR - Ek*(...)))

clc; clear; close all; addpath(genpath('./'));

%% === Output folders ===
outRoot = 'Data Plots and Mats';
if isfolder(outRoot)
    try
        delete(fullfile(outRoot, '**', '*.*')); pause(0.1); rmdir(outRoot,'s');
    catch, end
end
dataFolder = fullfile(outRoot,'Data_Fy');
plotFolder = fullfile(outRoot,'Plots_Fy');
mkdir(dataFolder); mkdir(plotFolder);

%% === Step 1: Load raw runs & combine ===
load('B2356run2.mat'); T1 = table(FZ, V, IA, P, FX, FY, SA, SR);
load('B2356run4.mat'); T2 = table(FZ, V, IA, P, FX, FY, SA, SR);
load('B2356run5.mat'); T3 = table(FZ, V, IA, P, FX, FY, SA, SR);
load('B2356run6.mat'); T4 = table(FZ, V, IA, P, FX, FY, SA, SR);
inputTable = [T1; T2; T3; T4];

%% === Step 2: Κατηγοριοποίηση ===
nRows = height(inputTable);
FZ_lb = inputTable.FZ * 0.2248089431;

catLoad = strings(nRows,1); catP = strings(nRows,1); catIA = strings(nRows,1); catV = strings(nRows,1);
for i = 1:nRows
    FZ_ = FZ_lb(i); P_ = inputTable.P(i); IA_ = inputTable.IA(i); V_ = inputTable.V(i);
    if  -55 < FZ_ && FZ_ <  -45, catLoad(i)="50";
    elseif -109< FZ_ && FZ_ < -91, catLoad(i)="100";
    elseif -160< FZ_ && FZ_ < -140,catLoad(i)="150";
    elseif -212< FZ_ && FZ_ < -188,catLoad(i)="200";
    elseif -265< FZ_ && FZ_ < -235,catLoad(i)="250";
    else, catLoad(i)="NA"; end

    if   49 < P_ && P_ <  61, catP(i)="8";
    elseif 62 < P_ && P_ <  76, catP(i)="10";
    elseif 77 < P_ && P_ <  90, catP(i)="12";
    elseif 91 < P_ && P_ < 102, catP(i)="14";
    else, catP(i)="NA"; end

    if   -0.7 < IA_ && IA_ < 0.7, catIA(i)="0";
    elseif 1.2 < IA_ && IA_ < 2.8, catIA(i)="2";
    elseif 3.2 < IA_ && IA_ < 4.8, catIA(i)="4";
    else, catIA(i)="NA"; end

    if   19 < V_ && V_ < 30, catV(i)="15";
    elseif 34 < V_ && V_ < 48, catV(i)="25";
    elseif 62 < V_ && V_ < 79, catV(i)="45";
    else, catV(i)="NA"; end
end
inputTable.catLoad = catLoad; inputTable.catP = catP; inputTable.catIA = catIA; inputTable.catV = catV;
save(fullfile(outRoot,'Lateral_Data.mat'),'inputTable');

%% === Step 3: Save per-combination & plots ===
validRows = catLoad ~= "NA" & catP ~= "NA" & catIA ~= "NA" & catV ~= "NA";
filteredTable = inputTable(validRows, :);

filteredTable = filteredTable(filteredTable.SA < 0 & filteredTable.FY > 0, :);

combinations = unique(filteredTable(:, ["catLoad","catP","catIA","catV"]));
for i = 1:height(combinations)
    c = combinations(i,:);
    subTable = filteredTable( ...
        filteredTable.catLoad==c.catLoad & ...
        filteredTable.catP==c.catP & ...
        filteredTable.catIA==c.catIA & ...
        filteredTable.catV==c.catV, :);
    nameTag = sprintf('L%s-P%s-IA%s-V%s', c.catLoad{1}, c.catP{1}, c.catIA{1}, c.catV{1});
    if ~isempty(subTable) && height(subTable)>=6
        save(fullfile(dataFolder,[nameTag '.mat']),'subTable');
        hfig = figure('Visible','off');
        scatter(subTable.SA, subTable.FY, 10, 'filled');
        xlabel('Slip Angle (deg)'); ylabel('Lateral Force FY (N)');
        title(['FY vs SA: ' nameTag]); grid on;
        saveas(hfig, fullfile(plotFolder,[nameTag '.png']));
        close(hfig);
    end
end

%% === Step 4: Συλλογή RAW για fitting ===
L_indx = [50;100;150;200;250]; P_indx=[8;10;12;14]; IA_indx=[0]; V_indx=[15;25;45];

allFY=[]; allSA=[]; allLoad=[]; allP=[]; allIA=[]; allSR=[];
for l=1:length(L_indx)
    for p=1:length(P_indx)
        for ia=1:length(IA_indx)
            for v=1:length(V_indx)
                nameTag = sprintf('L%d-P%d-IA%d-V%d', L_indx(l), P_indx(p), IA_indx(ia), V_indx(v));
                matFilePath = fullfile(dataFolder,[nameTag '.mat']);
                if exist(matFilePath,'file')
                    S = load(matFilePath);
                    sub = rmmissing(S.subTable);
                    sub = sub(sub.SA<0 & sub.FY>0,:);
                    if height(sub) < 8, continue; end

                    allSA   = [allSA;   sub.SA];
                    allFY   = [allFY;   sub.FY];
                    allLoad = [allLoad; repmat(L_indx(l),height(sub),1)];
                    allP    = [allP;    repmat(P_indx(p),height(sub),1)];
                    allIA   = [allIA;   repmat(IA_indx(ia),height(sub),1)];
                    allSR   = [allSR;   sub.SR];   % <-- χρησιμοποιούμε SR
                end
            end
        end
    end
end
fprintf('Collected %d raw sample points. SR nonzero? %d\n', numel(allFY), any(abs(allSR)>1e-6));

%% === Step 5: Inputs ===
Fz_N  = allLoad ./ 0.2248089431;        % lb->N
Pinput = allP * 6.89476 * 1000;         % psi->Pa
SA_rad = deg2rad(allSA);
SR     = allSR;

% Robust weights + ισορροπημένη έμφαση ανά |SA|
iq  = iqr(allFY); 
medFY = median(allFY);

w = ones(size(allFY));
w(abs(allFY - medFY) > 3*iq) = 0.25;      % outliers downweight

SAdeg = abs(allSA);

% ➊ ΔΩΣΕ ΜΕΓΑΛΥΤΕΡΟ ΒΑΡΟΣ ΣΤΗ ΓΡΑΜΜΙΚΗ ΖΩΝΗ (<=3°)
w = w .* (1 + 4.0*(SAdeg <= 3));

% ➋ ΗΠΙΟ ramp στην ουρά
w = w .* (1 + 0.02*SAdeg.^2);   % μικρότερη έμφαση καθώς μεγαλώνει |SA|


%% === Combined on/off ===
doCombined = any(abs(SR) > 1e-6);
fprintf('%s\n', tern(doCombined,'Running combined fit (with SR)...','No SR variation: only pure fit.'));

%% === GA / fmincon options ===
UseParallel = true;
gaOpts = optimoptions('ga','PopulationSize',220,'MaxGenerations',160,...
    'MaxStallGenerations',60,'FunctionTolerance',1e-6,'Display','iter','UseParallel',UseParallel);
fminOpts = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'FunctionTolerance',1e-6,'MaxIterations',2000,'MaxFunctionEvaluations',2e5);

%% === Pure fit ===
Dpure = 23;
lb_pure = -5*ones(1,Dpure);  ub_pure =  5*ones(1,Dpure);

% C πιο “ρεαλιστικό”
lb_pure(1)  = 1.2;  ub_pure(1)  = 1.8;

lb_pure(12) = -1.0;   ub_pure(12) = 3.0;   % Kya_scale ∈ [e^-1, e^3] ≈ [0.37, 20]

lb_pure(19:23) = -2.0;   % επιτρέπει αρκετά αρνητικό E (να πέφτει)
ub_pure(19:23) = 1;    % αν βαλω μικροτερη τιμη, π.χ 0.6, πέφτει καλύτερα

fitnessPure = @(P) MF62_RMS_FY_Cost_weighted(P, SA_rad, Fz_N, Pinput, allFY, w);
fprintf('Starting GA for pure fit (dim=%d)...\n', Dpure);
[Pga_pure, ~] = ga(fitnessPure, Dpure, [], [], [], [], lb_pure, ub_pure, [], gaOpts);
problem = createOptimProblem('fmincon','objective',fitnessPure,'x0',Pga_pure,'options',fminOpts);
[Ppure_opt, finalPure] = fmincon(problem);
fprintf('Pure fit done. Weighted RMS = %.4f\n', finalPure);

%% === Combined fit (με SR) ===
% 5 params για G_yk: [Ck, By0, By1(dfz), Ek, Sk]
if doCombined
    Dk   = 5;  Dtot = Dpure + Dk;
    lb   = [-5*ones(1,Dpure),  0.5, 0.01, -1.0, -0.99, -0.2];
    ub   = [ 5*ones(1,Dpure),  2.5, 5.00,  1.0,  0.99,  0.2 ];
    lb(1) = 1.1;  ub(1) = 1.9;
    lb(19:23) = -2.5;  ub(19:23) = 2.5;

    x0 = [Ppure_opt(:)'  1.4 0.5 0.0 0.0 0.0];
    fitnessCombined = @(P) MF62_RMS_FY_Combined_Cost_weighted(P, SA_rad, SR, Fz_N, Pinput, allFY, w);
    fprintf('Starting GA for combined fit (dim=%d)...\n', Dtot);
    [Pga_comb, ~] = ga(fitnessCombined, Dtot, [], [], [], [], lb, ub, [], gaOpts);
    problem = createOptimProblem('fmincon','objective',fitnessCombined,'x0',Pga_comb,'options',fminOpts);
    [Popt_comb, finalComb] = fmincon(problem);
    fprintf('Combined fit done. Weighted RMS = %.4f\n', finalComb);
else
    Popt_comb = [Ppure_opt(:)' 1.4 0.5 0.0 0.0 0.0]; % dummy append
end

%% === Step 6: Validation plot (δείξε COMBINED ξεκάθαρα) ===
exampleFile = fullfile(dataFolder, 'L50-P8-IA0-V25.mat');
if exist(exampleFile,'file'), S = load(exampleFile); subTable = S.subTable;
else, files = dir(fullfile(dataFolder,'L*.mat'));
    if ~isempty(files), S = load(fullfile(dataFolder,files(1).name)); subTable = S.subTable; else, subTable = table(); end
end

Slip_Angle = linspace(min(allSA), max(allSA), 600)'; % deg
%Fy_pure_pred = MF62_FY_vector(Ppure_opt, deg2rad(Slip_Angle), 222.411, 8*6.89476*1000);

figure('Color','w'); hold on;
if ~isempty(subTable), scatter(subTable.SA, subTable.FY, 15, 'b','filled'); end
%plot(Slip_Angle, Fy_pure_pred, 'r', 'LineWidth', 2);

%leg = {'Measured FY','Pure fit'};
leg = {'Measured FY'};
if doCombined
    % πάρε μια ρεαλιστική SR για το plot
    if ismember('SR', subTable.Properties.VariableNames)
        sr_plot = median(subTable.SR, 'omitnan');
    else
        tmp = SR(abs(SR)>1e-6); if isempty(tmp), sr_plot = 0.1; else, sr_plot = median(tmp); end
    end
    Fy_comb_plot = MF62_FY_combined(Popt_comb, deg2rad(Slip_Angle), sr_plot, 222.411, 8*6.89476*1000);
    plot(Slip_Angle, Fy_comb_plot, 'k', 'LineWidth', 1.8);
    leg{end+1} = sprintf('Combined fit (SR = %.3g)', sr_plot);
end
legend(leg,'Location','best');
xlabel('Slip Angle (deg)'); ylabel('Lateral Force FY (N)');
title('FY vs Slip Angle: Real vs Predicted'); grid on;

%% =================== Functions ===================
function Fy = MF62_FY_vector(P, SA, Fz, Pinput)
    ex  = 1e-8; Fz0 = 1112.055; p0 = 82737.090; IA = 0;
    dfz = (Fz - Fz0) ./ Fz0; dp  = (Pinput - p0) ./ p0;
    C   = max(0.5, P(1));
    D   = Fz .* (P(2) - abs(P(3)) .* dfz) .* (1 + P(4)*dp + P(5)*dp.^2) .* (1 - P(6)*IA.^2);
    Kyg0 = Fz .* (P(7) + P(8)*dfz) .* (1 + P(9)*dp);
    SVyg = Fz .* (P(10) + P(11)*dfz) .* IA; SVy  = SVyg;
    Kya_scale = exp(P(12));   % θετικό και με μεγάλο εύρος
    Kya = Kya_scale*Fz0 .* (1 + P(13)*dp) .* (1 - P(14)*abs(IA)) .* ...
      sin(P(15)*atan((Fz./Fz0) ./ ((P(16) + P(17)*IA.^2) .* (1 + P(18)*dp))));
    Shy = (Kyg0 .* IA - SVyg) ./ (Kya + ex);
    E   = (P(19) + P(20)*dfz) .* (1 + P(21)*IA.^2 - (P(22) + P(23)*IA) .* sign(SA + Shy));
    E   = max(-2.5, min(2.5, E));        % χαλαρό clamp → δυνατότητα downturn
    By  = Kya ./ (C .* D + ex);
    ay  = SA + Shy;
    Fy  = D .* sin(C .* atan(By .* ay - E .* (By .* ay - atan(By .* ay)))) + SVy;
end

function rmsError = MF62_RMS_FY_Cost_weighted(P, SA, Fz, Pinput, realFY, w)
    Fy_pred = MF62_FY_vector(P, SA, Fz, Pinput);
    err = (realFY - Fy_pred);
    rmsError = sqrt( sum(w .* err.^2) / max(1,sum(w)) );
    if ~isfinite(rmsError), rmsError = 1e12; end
end

function Fy = MF62_FY_combined(P, SA, SR, Fz, Pinput)
    % P(1:23): pure, P(24:28): [Ck, By0, By1(dfz), Ek, Sk]
    pureP = P(1:23); Fy0 = MF62_FY_vector(pureP, SA, Fz, Pinput);
    Ck = max(0.5,  P(24)); By0 = max(1e-6, P(25)); By1 = P(26);
    Ek = max(-0.99, min(0.99, P(27))); Sk = P(28);
    Fz0 = 1112.055; dfz = (Fz - Fz0)./Fz0; Byk = max(1e-6, By0 + By1*dfz);
    kx  = SR + Sk; phi = Byk .* kx;
    Gyk = cos( Ck .* atan( phi - Ek .* (phi - atan(phi)) ) );
    Fy  = Fy0 .* Gyk;
end

function rmsError = MF62_RMS_FY_Combined_Cost_weighted(P, SA, SR, Fz, Pinput, realFY, w)
    Fy_pred = MF62_FY_combined(P, SA, SR, Fz, Pinput);
    err = (realFY - Fy_pred);
    rmsError = sqrt( sum(w .* err.^2) / max(1,sum(w)) );
    if ~isfinite(rmsError), rmsError = 1e12; end
end

function s = tern(cond,a,b), if cond, s=a; else, s=b; end
end


