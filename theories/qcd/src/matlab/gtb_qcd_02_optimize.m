%% GTB QCD Optimization (CVX + MOSEK)
% Reads kernels and constants, runs a Watsonian iteration loop, saves outputs.

clear; clc;

%% ------------------------------- Paths ----------------------------------
here = fileparts(mfilename('fullpath'));                 % ...\theories\qcd\src\matlab
% Go up 4 levels: matlab -> src -> qcd -> theories -> (repo root)
repo_root = fileparts(fileparts(fileparts(fileparts(here))));
qcd_root  = fullfile(repo_root, 'theories', 'qcd');

matrices_dir = fullfile(qcd_root, 'data', 'matrices');   % inputs from Mathematica
results_dir  = fullfile(qcd_root, 'data', 'results');    % outputs from MATLAB

assert(isfolder(matrices_dir), 'Expected input dir not found: %s', matrices_dir);
if ~isfolder(results_dir), mkdir(results_dir); end

%% ------------------------- Parameters setup -----------------------------
unit   = 139.57; % pion mass as energy unit
fpi    = 92/unit; % fpi in pion mass unit
alphas = 0.327631; % strong interaction coupling at ~2GeV
FS0asym = 1;                   
FP1asym = 16*pi*alphas*fpi^2;
FD0asym = 48*pi*alphas*fpi^2;% Brodsky-Lepage asymptotic form factors

% Discretization
M   = 50;         % number of interpolation points
l   = 10;         % number of partial waves per isospin
Na  = 49;         % number of coefficients for spectral density
Mrho   = M^2 + M*(M+1)/2;
Msigma = 2*M + 1; % number of amplitude params: Mrho + Msigma

% Interpolation points on the physical line
nu0 = -20;
k   = (1:M);
vnu = nu0 + (8 - 2*nu0) ./ (1 + cos((k - 0.5)*pi/M));
s0      = (2000/unit)^2; % s0~2GeV
n0      = sum(vnu <= s0);
vnuhigh = vnu(n0:M);  % #interpolation points above s0~2GeV

% Unphysical region (for imposing chiSB constraints)
Munphys   = 4;
vsunphys  = 1/2 + ((1:Munphys)-1)/2;

% Numerical control parameters
eS0 = 0.05; eD0 = 0.1; eP1 = 0.1;      % sum-rule errors
rS0 = 8; rP1 = 2.2; rD0 = 8;           % asymptotic multipliers
aS0b = 10; aP1b = 100; aD0b = 1000;    % bounds on spectral density coefficients
echi = 2e-3;                           % low-energy matching error
nSR  = 6;                              % number of sum-rule moments
nWatson = 5;                           % number of Watsonian iterations

%% ----------------------- Load all kernels/data -------------------------
reqFiles = { ...
  'partial_waves.json', ...
  'chisb.json', ...
  'form_factors.json', ...
  'pion_coupling.json', ...
  'spectral_density.json', ...
  'sum_rules.json' ...
};
for i = 1:numel(reqFiles)
  f = fullfile(matrices_dir, reqFiles{i});
  assert(isfile(f), 'Missing required file: %s', f);
end

K   = jsondecode(fileread(fullfile(matrices_dir,'partial_waves.json')));
B1 = K.B1;   B2 = K.B2;   A1 = K.A1;   A2 = K.A2;   Bhat1 = K.Bhat1;   Bhat2 = K.Bhat2;

C   = jsondecode(fileread(fullfile(matrices_dir,'chisb.json')));
fsigmachi = C.fsigmachi;  frhochi = C.frhochi;

F   = jsondecode(fileread(fullfile(matrices_dir,'form_factors.json')));
Lambda3 = F.Lambda3;      KFF     = F.KFF3;

AS  = jsondecode(fileread(fullfile(matrices_dir,'pion_coupling.json')));
Asymsigma = AS.Asymsigma; Asymrho = AS.Asymrho;

SP  = jsondecode(fileread(fullfile(matrices_dir,'spectral_density.json')));
atorhoS0 = SP.atorhoS0;   atorhoP1 = SP.atorhoP1;   atorhoD0 = SP.atorhoD0;

SR  = jsondecode(fileread(fullfile(matrices_dir,'sum_rules.json')));
intS0 = SR.intS0;         intP1    = SR.intP1;      intD0    = SR.intD0; 
qcdfesr = SR.qcdfesr;

%% ----------------------- Watsonian iteration loop ----------------------
ImhtS2D2F1_prev = [];
RehtS2D2F1_prev = [];
htS0P1D0_prev   = [];
FF3_prev        = [];

for ni = 0:nWatson
    if ni == 0
        watson = [];
    else
        % Watsonian functional for S0,P1,D0
        htS0P1D0fnl = abs(htS0P1D0_prev) .* FF3_prev ./ abs(FF3_prev);
        watson.ImhtS2D2F1 = ImhtS2D2F1_prev;
        watson.RehtS2D2F1 = RehtS2D2F1_prev;
        watson.ImhtS0P1D0 = imag(htS0P1D0fnl);
        watson.RehtS0P1D0 = real(htS0P1D0fnl);
    end

    % ---- Solve one CVX problem ----
    [sol, diag] = optimize_core( ...
        ni, ...                       % iteration index
        M,l,Na,Mrho,Msigma, ...
        Munphys, ...
        vnu,n0,vnuhigh, ...
        B1,B2,A1,A2, ...
        Bhat1,Bhat2, ...
        fsigmachi,frhochi, ...
        Lambda3,KFF, ...
        Asymsigma,Asymrho, ...
        atorhoS0,atorhoP1,atorhoD0, ...
        intS0,intP1,intD0,qcdfesr, ...
        eS0,eP1,eD0,FS0asym,FP1asym,FD0asym,rS0,rP1,rD0, ...
        aS0b,aP1b,aD0b,echi,nSR, ...
        vsunphys, ...
        watson ...
    );

    % --- Logging --- 
    fprintf('%d %f %s\n', ni, sol.lambda, sol.cvx_status);
    
    % --- Build result struct --- 
    result = struct( ...
    'lambda',        sol.lambda, ...
    'RehtS0',        sol.Reht6(1:n0), ...
    'ImhtS0',        sol.Imht6(1:n0), ...
    'RehtD0',        sol.Reht6(n0+1:2*n0), ...
    'ImhtD0',        sol.Imht6(n0+1:2*n0), ...
    'RehtS2',        sol.Reht6(2*n0+1:3*n0), ...
    'ImhtS2',        sol.Imht6(2*n0+1:3*n0), ...
    'RehtD2',        sol.Reht6(3*n0+1:4*n0), ...
    'ImhtD2',        sol.Imht6(3*n0+1:4*n0), ...
    'RehtP1',        sol.Reht6(4*n0+1:5*n0), ...
    'ImhtP1',        sol.Imht6(4*n0+1:5*n0), ...
    'RehtF1',        sol.Reht6(5*n0+1:6*n0), ...
    'ImhtF1',        sol.Imht6(5*n0+1:6*n0), ...
    'specS0',        sol.spec3(1:n0), ...
    'specP1',        sol.spec3(n0+1:2*n0), ...
    'specD0',        sol.spec3(2*n0+1:3*n0), ...
    'ReFFS0',        sol.ReFF(1:M), ...
    'ImFFS0',        sol.ImFF(1:M), ...
    'ReFFP1',        sol.ReFF(M+1:2*M), ...
    'ImFFP1',        sol.ImFF(M+1:2*M), ...
    'ReFFD0',        sol.ReFF(2*M+1:3*M), ...
    'ImFFD0',        sol.ImFF(2*M+1:3*M) ...
    );

    % --- Output ---
    jsonText = jsonencode(result);    % turn struct into JSON text
    fid = fopen(fullfile(results_dir,sprintf('gtb_qcd_W%02d.json', ni)),'w');
    fwrite(fid, jsonText, 'char');
    fclose(fid);
    
    % --- Prep previous-iteration info (for next Watson step) ---
    ImhtS2D2F1_prev = diag.ImhtS2D2F1;
    RehtS2D2F1_prev = diag.RehtS2D2F1;
    htS0P1D0_prev   = diag.htS0P1D0;
    FF3_prev        = diag.FF3;
end
