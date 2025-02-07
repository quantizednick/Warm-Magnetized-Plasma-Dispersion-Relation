% Ex2 in Wave Propagation in Ionized Medium, Nikolaos Saltas, ECE NTUA
% Solving the Dispersion Relation in Warm Magnetized Plasma with Propagation at an
% angle theta to the magnetic field.

%refractive index n relative to temperature T for constant theta
clear; clf(figure(1)); clc;
%% Define Constant Plasma Parameters
e0 = 8.8541878128e-12; % vacuum permittivity (F/m)
c = 3e8; % speed of light in vacuum (m/s)
me = 9.1093837e-31; % electron mass (kg)
qe = 1.60217663e-19; % electron charge (C)
KB = 1; % Boltzmann constant Normalized as we incert it in T

% Define Values Related to Specific Problem
Bo = 5; % Magnetic flux density (T)
ne = 1e19; % Electron density (m^-3)
f = 100e9; % Frequency (Hz)
theta = 0; % Angle range from 0 to pi/2
T_ = linspace(0, 10, 200); % Temperature (KeV)
T = T_*1000;

% Calculate Other Useful Parameters
costh = cos(theta);
sinth = sin(theta);
omega = 2 * pi * f; % Frequency in rad/s
omegace = Bo * qe / me; % Electron cyclotron frequency
omegape = sqrt(ne * qe^2 / (me * e0)); % Plasma frequency
Vthef = @(T) sqrt(KB * T*1.6e-16 / me); % Thermal velocity

%% Solving for each theta
%Initialize arrays for storing solutions
nsq_1_real = zeros(size(theta));
nsq_2_real = zeros(size(theta));
nsq_1_imag = zeros(size(theta));
nsq_2_imag = zeros(size(theta));

% Loop over each theta value to compute n^2
for i = 1:length(T)
    Vthe = Vthef(T(i));
    
% Dielectric Tensor Elements considering k^2 = n^2*(omega/c)^2 
exx = @(nsq)  1 - omegape^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2) ./ ...
              (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2));
exy = @(nsq)  -1i * (omegace/omega) * omegape^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2) ./ ...
              (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2));
exz = @(nsq)  -1i * omegape^2 * nsq*(omega/c)^2 * Vthe^2 * costh * sinth ./ ...
              (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2));
eyx = @(nsq)  -exy(nsq); %Hermitian Property
eyy = @(nsq)  1 - omegape^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) ./ ...
              (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2));
eyz = @(nsq)  -1i * (omegace/omega) * omegape^2 * nsq*(omega/c)^2 * Vthe^2 * costh * sinth ./ ...
              (omega * (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2)));
ezx = @(nsq)  exz(nsq); %Hermitian Property
ezy = @(nsq)  -eyz(nsq); %Hermitian Property
ezz = @(nsq)  1 - omegape^2 * (omega^2 - omegace^2 - nsq*(omega/c)^2 * Vthe^2 * sinth^2) ./ ...
              (omega^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2) - omegace^2 * (omega^2 - nsq*(omega/c)^2 * Vthe^2 * costh^2));
%Dtensor = [exx exy exz; eyx eyy eyz; ezx ezy ezz];
%Using the Helmholtz Equation we determine the dispersion relation is a
%biquadratic equation. So we will get two solutions(modes) for n^2. For theta = 0
%we expect one to transition to L and the other to R, and for theta = pi/2,
%we expect one to transition to X and the other to O respectively

%An^4 + Bn^2 + C = 0
A = @(nsq) ezz(nsq) .* costh^2 + exx(nsq) .* sinth^2 + 2 .* sinth .* costh * exz(nsq);

B = @(nsq) -ezz(nsq) * exx(nsq) - ezz(nsq) * eyy(nsq) * costh^2 - exx(nsq) * eyy(nsq) * sinth^2 - exy(nsq)^2 * sinth^2 ...
           - eyz(nsq) * exz(nsq) * costh^2 + sinth * costh * exy(nsq) * exz(nsq) - sinth * costh * eyy(nsq) * exz(nsq) ...
           + exz(nsq)^2 + sinth * costh * exy(nsq) * eyz(nsq) - sinth * costh * eyy(nsq) * exz(nsq);

C = @(nsq) exx(nsq) * eyy(nsq) * ezz(nsq) + ezz(nsq) * exy(nsq)^2 + exx(nsq) * eyz(nsq) * exz(nsq)...
           + exy(nsq) * exz(nsq)^2 + exy(nsq) * exz(nsq) * eyz(nsq) - eyy(nsq) * exz(nsq)^2;

% Compute the discriminant and define equations to solve
discriminant = @(nsq) B(nsq)^2 - 4 * A(nsq) * C(nsq);
fun1 = @(nsq) (-B(nsq) + sqrt(discriminant(nsq))) / (2 * A(nsq)) - nsq;
fun2 = @(nsq) (-B(nsq) - sqrt(discriminant(nsq))) / (2 * A(nsq)) - nsq;
n0 = 1;
options = optimoptions('fsolve', 'Display', 'off');
Sol1 = fsolve(fun1, n0, options);
Sol2 = fsolve(fun2, n0, options);
% Store real and imaginary parts separately
nsq_1_real(i) = real(Sol1);
nsq_2_real(i) = real(Sol2);
nsq_1_imag(i) = imag(Sol1);
nsq_2_imag(i) = imag(Sol2);
end

%% Plot the real parts of the two modes
figure(1);
plot(T, sqrt(nsq_1_real), 'b', 'LineWidth', 2); hold on;
plot(T, sqrt(nsq_2_real), 'r', 'LineWidth', 2);
xlabel('T(eV)');
ylabel('Re(n^2)');
title('Real Part of n relative to Temperature');
legend('Mode 1', 'Mode 2');
grid on;

%Plot the imaginary parts of the two modes
%figure(2);
%plot(T, nsq_1_imag, 'b', 'LineWidth', 2); hold on;
%plot(T, nsq_2_imag, 'r', 'LineWidth', 2);
%xlabel('T(KeV)');
%ylabel('Im(n^2)');
%title('Imaginary Part of n^2 vs. Temperature');
%legend('Mode 1', 'Mode 2');
%grid on;
