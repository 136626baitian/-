% Rayleigh extinction & backscatter profile @1030nm
clear; clc;

% Constants
kB = 1.380649e-23;        % Boltzmann constant (J/K)
n_air = 1.0003;           % Refractive index of air
rho = 0.03;               % Depolarization factor
N0 = 2.547e25;            % Molecular number density at STP (m^-3)
lambda = 1030e-9;         % Wavelength in meters

% Altitude range (0 - 50 km)
h = linspace(0, 50000, 500); % Altitude in meters
T0 = 288.15;    % Sea level temperature (K)
P0 = 101325;    % Sea level pressure (Pa)
L = -0.0065;    % Temperature lapse rate (K/m)
R = 287.05;     % Gas constant (J/kg/K)
g = 9.80665;    % Gravitational acceleration (m/s^2)

% Temperature and pressure profile (troposphere approximation)
T = T0 + L * h;
P = P0 * (T ./ T0).^(-g / (L * R));

% Molecular number density profile
N = P ./ (kB * T);  % molecules/m^3

% Rayleigh scattering cross section (per molecule)
sigma_rayleigh = (24 * pi^3) / (lambda^4 * N0^2) * ...
                 ((n_air^2 - 1)/(n_air^2 + 2))^2 * ...
                 (6 + 3 * rho) / (6 - 7 * rho);

% Extinction coefficient (1/m)
alpha_rayleigh = N * sigma_rayleigh;

% Backscatter coefficient (1/m/sr)
beta_rayleigh = (3 / (4 * pi)) * alpha_rayleigh;

% Plot
figure;
plot(alpha_rayleigh, h/1000, 'b-', 'LineWidth', 1.5); hold on;
plot(beta_rayleigh, h/1000, 'r--', 'LineWidth', 1.5);
xlabel('Coefficient (1/m or 1/m/sr)');
ylabel('Altitude (km)');
title('Rayleigh Extinction and Backscatter at 1030 nm');
legend('\alpha (extinction)', '\beta (backscatter)');
grid on;
