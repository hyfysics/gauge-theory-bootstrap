%%%%%%%%%%%%%%%%%%%
% numerical setup %
%%%%%%%%%%%%%%%%%%%

mpi=139.57; fpi=92; % fpi and mpi

M = 50; % number of interpolation points on the physical line
l = 10; % number of partial waves per isospin
nB = 40; % number of BSpline basis functions
Mrho = M^2+M*(M+1)/2; Msigma = 2*M+1; % number of amplitude parameters
nu0 = -20; % nu0 maps to the center of the disk
vnu = nu0+(8-2*nu0)./(1+cos(([1:M]-1/2)*pi/M)); % list of interlation points s_j
s0 = (2000/mpi)^2; % s0~2GeV
n0 = sum(double(vnu<=s0)); % number of interpolation points below s0
vnuhigh = vnu(n0+1:M); % list of interlation points above s0

% for computing physical region partial waves
B1 = importdata('B1.mat'); A1 = importdata('A1.mat');
B2 = importdata('B2.mat'); A2 = importdata('A2.mat');
Bhat1 = importdata('Bhat1.mat'); Bhat2 = importdata('Bhat2.mat'); 

% for evaluating unphysical region partial waves: functional
fsigmaF = importdata('fsigmaF.mat'); frhoF = importdata('frhoF.mat');
% for evaluating unphysical region partial waves: chiSB constraints
Munphys=4; vsunphys=1/2+([1:Munphys]-1)/2;
fsigmachi = importdata('fsigmachi.mat'); frhochi = importdata('frhochi.mat');

% form factors
Lambda3 = importdata('Lambda3.mat'); % rescaling factors for ImFhat
KFF = importdata('KFF3.mat'); % dispersion relations kernel

% spectral densities
MBSpline = importdata('MBSpline3.mat'); % BSpline basis functions
intS0 = importdata('intS0.mat'); intP1 = importdata('intP1.mat'); intD0 = importdata('intD0.mat'); % sum rules numerical integration
qcdfesr = importdata('SVZSR.mat'); % sum rules from QCD

% numerical parameters
echi = 2*10^(-3); % low energy matching tolerance
eS0=1*10^(-7); eD0=5*10^(-6); eP1=6*10^(-6); % sum rule errors
FP1asym = 6.87076; FD0asym = 20.6123; %asymptotics of s*F(s) from QCD
rS0 = 5*10^(-2); rP1 = 2; rD0 = 6; %relaxation factors for bounding asymptotic s*F(s)

%take a point near Weinberg model
ni = 46; %ni=44,46,48 for blue, green, red
F0fix = ni/50*0.163108; %fix F0, maximize F1

%%%%%%%%%%%%%%%%%%%%%%%
% initial maxmization %
%%%%%%%%%%%%%%%%%%%%%%%

cvx_begin sdp
cvx_solver mosek
variable rho(Mrho) % double spectral density of amplitudes
variable sig(Msigma) % single spectral density of amplitudes
variable specC(3*nB) % coefficients for rhohat
variable ImFF(3*M) % imaginary part of form factors
variable Bmat(3,3,3*n0) hermitian % positive semidefinite matrices

% regularization
norm(rho/Mrho,4) <= 10^(2);

% partial waves unitarity as quadratic cones
Imht = B2*rho+B1*sig; Reht = A2*rho+A1*sig; ht = Reht+1i*Imht; % htilde
Imhh = Bhat2*rho+Bhat1*sig; % Imhhat
norms([Reht,Imht],2,2)<=sqrt(2*Imhh); % unitarity for all partial waves at all energies

%6 partial waves up to s0
indS0P1D0 = [1:n0,2*l*M+1:2*l*M+n0,M+1:M+n0]; % indices for S0, P1, D0
indS2D2F1 = [l*M+1:l*M+n0,l*M+M+1:l*M+M+n0,2*l*M+M+1:2*l*M+M+n0]; % indices for S2, D2, F1
htS0P1D0 = ht(indS0P1D0); ImhtS0P1D0 = imag(htS0P1D0); RehtS0P1D0 = real(htS0P1D0); ImhhS0P1D0 = Imhh(indS0P1D0); %S0,P1,D0
htS2D2F1 = ht(indS2D2F1); ImhtS2D2F1 = imag(htS2D2F1); RehtS2D2F1 = real(htS2D2F1); ImhhS2D2F1 = Imhh(indS2D2F1); %S2,D2,F1

% S0,P1,D0 rescaled spectral densities up to s0
spec = MBSpline*specC; specC >= 0;
specC(1)==1; specC(nB+1)==1; specC(2*nB+1)==1; % boundary condition rhohat(0)=1

% S0,P1,D0 form factors up to s0
ReFF = KFF*ImFF+1; FF = ReFF+1i*ImFF;
FF3 = [FF(1:n0);FF(M+1:M+n0);FF(2*M+1:2*M+n0)];
ImFFh3 = [ImFF(1:n0);ImFF(M+1:M+n0);ImFF(2*M+1:2*M+n0)]./Lambda3;

% positive semidefinite matrix: S0,P1,D0 up to s0
for i=[1:3*n0]
    Bmat(:,:,i) == [1,htS0P1D0(i),FF3(i); conj(htS0P1D0(i)),2*ImhhS0P1D0(i),2*ImFFh3(i); conj(FF3(i)),2*ImFFh3(i),spec(i)];
    Bmat(:,:,i)>=0;
end

% QCD sum rules
fesr=[intS0*spec(1:n0);intP1*spec(n0+1:2*n0);intD0*spec(2*n0+1:3*n0)]; 
w=fesr-qcdfesr; wS0=w(1:3); wP1=w(4:6); wD0=w(7:9);
norm(wS0) <= eS0; norm(wP1) <= eP1; norm(wD0) <= eD0; % sum rule errors
    
% form factor asymptotic bounds
abs(FF(n0+1:M)) <= rS0; abs(vnuhigh'.*FF(M+n0+1:2*M)) <= rP1*FP1asym; abs(vnuhigh'.*FF(2*M+n0+1:3*M)) <= rD0*FD0asym;

% chiSB constraints
fchi = fsigmachi*sig + frhochi*rho; % partial waves in unphysical region
fchiS0 = fchi(1:Munphys); fchiD0 = fchi(Munphys+1:2*Munphys);
fchiS2 = fchi(l*Munphys+1:l*Munphys+Munphys); fchiD2 = fchi(l*Munphys+Munphys+1:l*Munphys+2*Munphys);
fchiP1 = fchi(2*l*Munphys+1:2*l*Munphys+Munphys); fchiF1 = fchi(2*l*Munphys+Munphys+1:2*l*Munphys+2*Munphys);
R01tree = 3*(2*vsunphys-1)./(vsunphys-4); R21tree = 3*(2-vsunphys)./(vsunphys-4); % ratio from tree level
norm([fchiS0-R01tree'.*fchiP1;fchiS2-R21tree'.*fchiP1]) <= echi; % low energy matching
norm([fchiD0;fchiD2;fchiF1]) <= echi;

% initial functionals 
f3 = fsigmaF*sig + frhoF*rho; % evaluating partial waves at s*=3
f3S0 = f3(1); f3P1 = f3(2*l+1); f3D0 = f3(2); f3F1 = f3(2*l+2);
F0 = 2*(f3S0+5*f3D0); F1 = 2*(3*f3P1+7*f3F1); % two functionals projection
F0 == F0fix;  maximize(F1); %fix F0, maximize F1

cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Watsonian unitarization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Watsonian unitarization functionals
ImhtS2D2F1fnl = ImhtS2D2F1; RehtS2D2F1fnl = RehtS2D2F1;
ImhtS0P1D0fnl = 2*imag(FF3).*ImFFh3./(FF3.*conj(FF3));
RehtS0P1D0fnl = 2*real(FF3).*ImFFh3./(FF3.*conj(FF3));

cvx_begin sdp
cvx_solver mosek
variable rho(Mrho) % double spectral density of amplitudes
variable sig(Msigma) % single spectral density of amplitudes
variable specC(3*nB) % coefficients for rhohat
variable ImFF(3*M) % imaginary part of form factors
variable Bmat(3,3,3*n0) hermitian % positive semidefinite matrices

% regularization
norm(rho/Mrho,4) <= 10^(2);

% quadratic cone of unitarity for all partial waves at all energy
Imht = B2*rho+B1*sig; Reht = A2*rho+A1*sig; % Rehtilde, Imhtilde
ht = Reht+1i*Imht;
Imhh = Bhat2*rho+Bhat1*sig; % Imhhat
norms([Reht,Imht],2,2)<=sqrt(2*Imhh); % unitarity

%6 partial waves up to s0
indS0P1D0 = [1:n0,2*l*M+1:2*l*M+n0,M+1:M+n0]; % indices for S0, P1, D0
indS2D2F1 = [l*M+1:l*M+n0,l*M+M+1:l*M+M+n0,2*l*M+M+1:2*l*M+M+n0]; % indices for S2, D2, F1
htS0P1D0 = ht(indS0P1D0); ImhtS0P1D0 = imag(htS0P1D0); RehtS0P1D0 = real(htS0P1D0); ImhhS0P1D0 = Imhh(indS0P1D0); %S0,P1,D0
htS2D2F1 = ht(indS2D2F1); ImhtS2D2F1 = imag(htS2D2F1); RehtS2D2F1 = real(htS2D2F1); ImhhS2D2F1 = Imhh(indS2D2F1); %S2,D2,F1

% S0,P1,D0 rescaled spectral densities up to s0
spec = MBSpline*specC; specC >= 0;
specC(1)==1; specC(nB+1)==1; specC(2*nB+1)==1; % boundary condition rhohat(0)=1

% S0,P1,D0 form factors up to s0
ReFF = KFF*ImFF+1; FF = ReFF+1i*ImFF;
FF3 = [FF(1:n0);FF(M+1:M+n0);FF(2*M+1:2*M+n0)];
ImFFh3 = [ImFF(1:n0);ImFF(M+1:M+n0);ImFF(2*M+1:2*M+n0)]./Lambda3;

% positive semidefinite matrix: S0,P1,D0 up to s0
for i=[1:3*n0]
    Bmat(:,:,i) == [1,htS0P1D0(i),FF3(i); conj(htS0P1D0(i)),2*ImhhS0P1D0(i),2*ImFFh3(i); conj(FF3(i)),2*ImFFh3(i),spec(i)];
    Bmat(:,:,i)>=0;
end

% QCD sum rules
fesr=[intS0*spec(1:n0);intP1*spec(n0+1:2*n0);intD0*spec(2*n0+1:3*n0)]; 
w=fesr-qcdfesr; wS0=w(1:3); wP1=w(4:6); wD0=w(7:9);
norm(wS0) <= eS0; norm(wP1) <= eP1; norm(wD0) <= eD0; % sum rule errors
    
% form factor asymptotic bounds
abs(FF(n0+1:M)) <= rS0; abs(vnuhigh'.*FF(M+n0+1:2*M)) <= rP1*FP1asym; abs(vnuhigh'.*FF(2*M+n0+1:3*M)) <= rD0*FD0asym;

% chiSB constraints
fchi = fsigmachi*sig + frhochi*rho; % partial waves in unphysical region
fchiS0 = fchi(1:Munphys); fchiD0 = fchi(Munphys+1:2*Munphys);
fchiS2 = fchi(l*Munphys+1:l*Munphys+Munphys); fchiD2 = fchi(l*Munphys+Munphys+1:l*Munphys+2*Munphys);
fchiP1 = fchi(2*l*Munphys+1:2*l*Munphys+Munphys); fchiF1 = fchi(2*l*Munphys+Munphys+1:2*l*Munphys+2*Munphys);
R01tree = 3*(2*vsunphys-1)./(vsunphys-4); R21tree = 3*(2-vsunphys)./(vsunphys-4); % ratio from tree level
norm([fchiS0-R01tree'.*fchiP1;fchiS2-R21tree'.*fchiP1]) <= echi; % low energy matching
norm([fchiD0;fchiD2;fchiF1]) <= echi;

%Watsonian unitarization functional
v1 = ImhtS2D2F1fnl'*ImhtS2D2F1+RehtS2D2F1fnl'*RehtS2D2F1-sum(ImhhS2D2F1);
v2 = ImhtS0P1D0fnl'*ImhtS0P1D0+RehtS0P1D0fnl'*RehtS0P1D0-sum(ImhhS0P1D0);
maximize(v1+v2);

cvx_end

%%%%%%%%%%%%%%%%%%%%
% phase shift plot %
%%%%%%%%%%%%%%%%%%%%

deg=180/pi; %phase shift in unit of degree
%phase shift of the 6 partial waves
deltaS0=angle(htS0P1D0(1:n0))*deg; deltaD0=angle(htS0P1D0(2*n0+1:3*n0))*deg;
deltaS2=angle(htS2D2F1(1:n0))*deg; deltaD2=1/2*angle(1+1i*htS2D2F1(n0+1:2*n0))*deg;
deltaP1=angle(htS0P1D0(n0+1:2*n0)); deltaF1=1/2*angle(1+1i*htS2D2F1(2*n0+1:3*n0))*deg;

%plot S0,D0,S2,D2,P1,F1 phase shifts
figure
subplot(3,2,1)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaS0,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('S0 phase shift');
subplot(3,2,2)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaD0,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('D0 phase shift');
subplot(3,2,3)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaS2,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('S2 phase shift');
subplot(3,2,4)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaD2,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('D2 phase shift');
subplot(3,2,5)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaP1,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('P1 phase shift');
subplot(3,2,6)
plot(sqrt(vnu(1:n0))*mpi/1000,deltaF1,'go','MarkerSize',2);
xlabel('energy (GeV)');
ylabel('F1 phase shift');