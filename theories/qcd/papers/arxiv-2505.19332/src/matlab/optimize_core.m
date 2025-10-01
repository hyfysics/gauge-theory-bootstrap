function [sol, diag] = optimize_core( ...
    ni, M,l,Na,Mrho,Msigma, ...
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
    watson  ...
)
%OPTIMIZE_CORE One CVX solve of the GTB program.
% Returns:
%   sol  : struct with fields (lambda, cvx_status, sigmarho, Reht6, Imht6,
%          spec3, ReFF, ImFF).
%   diag : diagnostics for Watson iteration (Im/Re partial waves, FF3, etc.)

    % -------- Indices for partial-wave blocks below s0 ----------
    indS0P1D0 = [1:n0, 2*l*M+(1:n0), M+(1:n0)];
    indS2D2F1 = [l*M+(1:n0), l*M+M+(1:n0), 2*l*M+M+(1:n0)];

    % --------- CVX ----------
    cvx_begin sdp quiet
    cvx_solver mosek

    variable rho(Mrho)
    variable sig(Msigma)
    variable ImFF(3*M)
    variable Bmat(3,3,3*n0) hermitian
    variable a(3*Na)

    % regularization
    norm(rho/Mrho,4) <= 1e2;

    % partial waves (physical region) & unitarity
    Imht = B2*rho + B1*sig;
    Reht = A2*rho + A1*sig;
    ht   = Reht + 1i*Imht;
    Imhh = Bhat2*rho + Bhat1*sig;
    norms([Reht,Imht],2,2) <= sqrt(2*Imhh);

    % split by channels up to s0
    htS0P1D0 = ht(indS0P1D0);
    ImhtS0P1D0 = imag(htS0P1D0);
    RehtS0P1D0 = real(htS0P1D0);
    ImhhS0P1D0 = Imhh(indS0P1D0);

    htS2D2F1 = ht(indS2D2F1);
    ImhtS2D2F1 = imag(htS2D2F1);
    RehtS2D2F1 = real(htS2D2F1);
    ImhhS2D2F1 = Imhh(indS2D2F1);

    % spectral densities from coefficients a
    aS0 = a(1:Na);
    aP1 = a(Na+1:2*Na);
    aD0 = a(2*Na+1:3*Na);
    spec = [atorhoS0*aS0; atorhoP1*aP1; atorhoD0*aD0];

    % spectral densities coefficient bounds (for imposing low energy behaviors)
    abs(aS0(2:Na)) <= aS0b; aS0(1) <= 10; aS0(1) >= 1;
    abs(aP1(2:Na)) <= aP1b; aP1(1) <= 10; aP1(1) >= 1;
    abs(aD0(2:Na)) <= aD0b; aD0(1) <= 10; aD0(1) >= 1;

    % Form factors (S0,P1,D0) up to all energies
    ReFF = KFF*ImFF + ones(3*M,1);
    FF   = ReFF + 1i*ImFF;

    % Form factors (S0,P1,D0) up to s0
    FF3    = [FF(1:n0); FF(M+1:M+n0); FF(2*M+1:2*M+n0)];
    ImFFh3 = [ImFF(1:n0); ImFF(M+1:M+n0); ImFF(2*M+1:2*M+n0)] ./ Lambda3;

    % PSD matrices per physical point (S0,P1,D0 up to s0)
    for i = 1:(3*n0)
        Bmat(:,:,i) == [1,                  htS0P1D0(i),                FF3(i); ...
                        conj(htS0P1D0(i)),  2*ImhhS0P1D0(i),            2*ImFFh3(i); ...
                        conj(FF3(i)),       2*ImFFh3(i),                spec(i)];
        Bmat(:,:,i) >= 0;
    end

    % QCD sum rules
    fesr = [intS0*aS0; intP1*aP1; intD0*aD0];
    w    = fesr - qcdfesr;
    wS0  = w(1:nSR);
    wP1  = w(nSR+1:2*nSR);
    wD0  = w(2*nSR+1:3*nSR);
    abs(wS0) <= eS0;
    abs(wP1) <= eP1;
    abs(wD0) <= eD0;

    % Form factor asymptotics (above s0)
    abs( (vnuhigh ./ log(vnuhigh))' .* FF(n0:M) )  <= rS0*FS0asym;
    abs( (vnuhigh)' .* FF(M+n0:2*M) )              <= rP1*FP1asym;
    abs( (vnuhigh)' .* FF(2*M+n0:3*M) )            <= rD0*FD0asym;

    % chiSB constraints (unphysical region)
    fchi    = fsigmachi*sig + frhochi*rho;
    fchiS0  = fchi(1:Munphys);
    fchiS2  = fchi(l*Munphys + (1:Munphys));
    fchiP1  = fchi(2*l*Munphys + (1:Munphys));

    R01tree = 3*(2*vsunphys - 1) ./ (vsunphys - 4);
    R21tree = 3*(2 - vsunphys) ./ (vsunphys - 4);

    norm([ fchiS0 - R01tree'.*fchiP1 ;  fchiS2 - R21tree'.*fchiP1 ]) <= echi;

    % Symmetric point amplitude & pion coupling
    Asym   = Asymsigma'*sig + Asymrho'*rho;
    lambda = (3*pi*Asym)/4;

    % Objective: feasibility (ni=0) or Watson iterations (ni>0)
    if ni == 0
        minimize(0)
    else
        v1 = watson.ImhtS2D2F1' * ImhtS2D2F1 + watson.RehtS2D2F1' * RehtS2D2F1 - sum(ImhhS2D2F1);
        v2 = watson.ImhtS0P1D0' * ImhtS0P1D0 + watson.RehtS0P1D0' * RehtS0P1D0 - sum(ImhhS0P1D0);
        maximize(v1 + v2)
    end

    cvx_end

    % ------------------- Pack outputs -------------
    % amplitude params
    sigmarho = [sig; rho];

    % 6 partial waves (Re/Im) below s0, ordered as S0,D0,S2,D2,P1,F1
    Imht6 = [ ...
        Imht(1:n0); ...
        Imht(M+(1:n0)); ...
        Imht(l*M+(1:n0)); ...
        Imht(l*M+M+(1:n0)); ...
        Imht(2*l*M+(1:n0)); ...
        Imht(2*l*M+M+(1:n0)) ...
    ];
    Reht6 = [ ...
        Reht(1:n0); ...
        Reht(M+(1:n0)); ...
        Reht(l*M+(1:n0)); ...
        Reht(l*M+M+(1:n0)); ...
        Reht(2*l*M+(1:n0)); ...
        Reht(2*l*M+M+(1:n0)) ...
    ];

    % 3 spectral densities (S0,P1,D0) up to s0
    spec3 = spec(1:3*n0);

    % 3 form factors (S0,P1,D0) at all energies
    ReFF = real(ReFF);
    ImFF = imag(FF);

    % ----- Return structs -----
    sol = struct( ...
        'lambda', lambda, ...
        'cvx_status', cvx_status, ...
        'sigmarho', sigmarho, ...
        'Reht6', Reht6, ...
        'Imht6', Imht6, ...
        'spec3', spec3, ...
        'ReFF', ReFF, ...
        'ImFF', ImFF ...
    );

    diag = struct( ...
        'ImhtS2D2F1', ImhtS2D2F1, ...
        'RehtS2D2F1', RehtS2D2F1, ...
        'htS0P1D0',   htS0P1D0, ...
        'FF3',        FF3 ...
    );
end
