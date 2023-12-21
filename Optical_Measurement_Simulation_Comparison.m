% Script to plot the modeled (simulated MZI transfer function) versus the measured data

% Author: Emmanuel Santacruz
% Date:         2023-12-10
% Last Updated: 2023-12-10

folder = 'C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\OPTICAL MEASUREMENTS\CHIP 4';

% File Naming Convention
% Unbalanced_(number) -> number is the path length difference between the two branches
% (f\b) -> two iterations of the same measurement

% Test Names:
% esantacr_MZI_Unbalanced_12_f.mat
% esantacr_MZI_Unbalanced_290_f.mat
% esantacr_MZI_Unbalanced_590_f.mat

% Modelling Using Lumerical neff approximations
lambda0 = 1550e-9;           % Central Wavelength [m]
% Coefficients given by Taylor Expansion from Lumerical MODE
n1 = 2.476350025193613;     
n2 = -1.063976189658311;
n3 = 0.022030567074899;
neff = @(wavelength) n1 + n2.*((wavelength-lambda0)*10^(0))+n3.*((wavelength-lambda0)*10^(0)).^2;
% Complex Propagation
alpha = 1.000001011;    % Propagation Loss Coeff [m^-1] per meter
beta = @(wavelength) (2*pi*neff(wavelength)./(wavelength.*(10^(0))) - 1i*alpha/2*ones(1,length(wavelength)));
% MZI Transfer Function
T_MZI= @(deltaL,wavelength) (0.25*abs(1+exp(-1i*beta(wavelength)*deltaL)).^2);
T_MZI_dBm = @(deltaL,wavelength) 10.*log10(T_MZI(deltaL,wavelength));


% First we open the GC reference file to have on hand

file = "esantacr_Reference_GC_f.mat";
filename = fullfile(folder,file);

if ~exist(filename, "file")
    message = sprintf("%s does not exist",filename);
    uiwait(warndlg(message));
end

if contains(filename,"_f.mat","IgnoreCase",false)
    power_channel = 2;
elseif contains(filename,"_b.mat","IgnoreCase",false)
    power_channel = 3;
end

data = load(filename);
wavelength_space_measurement = data.scandata.wavelength;
gc_power_spectrum_measurement = data.scandata.power(:,power_channel);

% Open the GC reference simulated on lumerical interconnect

file = strrep(file,"_f.mat","_f_interconnect.mat");
filename = fullfile(folder,file);

if ~exist(filename, "file")
    message = sprintf("%s does not exist",filename);
    uiwait(warndlg(message));
end

data = load(filename);
wavelength_space_lum = 10^(-9).*data.lum.x0; % Convert to [m]
gc_power_spectrum_lum = flip(data.lum.y0);

% Open MZI file and extract data

file = "esantacr_MZI_Unbalanced_590_f.mat";
filename = fullfile(folder, file);

if ~exist(filename, "file")
    message = sprintf("%s does not exist",filename);
    uiwait(warndlg(message));
end

% Path Length Difference
% I have found that the FSR of the optical measurements is much smaller than
% what it is predicted to be, from the following one can choose what it
% should be, or an approximation to the measurements taken.

%deltaL = str2num(regexprep(file,{'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '}))*10^(-6);
deltaL = 1050*10^(-6); % Correct for 590_f.mat
%deltaL = 535*10^(-6); % Correct for 290_f.mat

if ~exist(filename,"file")
    message = sprintf("%s does not exist", filename);
    uiwait(warndlg(message));
end

if contains(filename,"_f.mat","IgnoreCase",false)
    power_channel = 2;
elseif contains(filename,"_b.mat","IgnoreCase",false)
    power_channel = 3;
end

data = load(filename);
wavelength_space_measurement = data.scandata.wavelength;
power_space_measuremet = data.scandata.power(:,power_channel);

% Determine the min index and max index for which the power is less the power threshold

[max_power,max_ind] = max(power_space_measuremet);     % Find the maximum power and its index
N_measurement = length(wavelength_space_measurement);
window_size = ceil(N_measurement/20);                         % Points to the left and points to the right
ind_1 = 1;                                  % Lower index
ind_2 = 1;                                  % Upper index

% Determine the correct vlaues for the indecies to chop-off 
if max_ind > window_size
    ind_1 = max_ind - window_size;
else 
    ind_1 = 1;
end

if max_ind + window_size > N_measurement
    ind_2 = N_measurement;
else
    ind_2 = max_ind + window_size;
end

% Wavelengths to Plot
min_wavelength = wavelength_space_measurement(ind_1);
max_wavelength = wavelength_space_measurement(ind_2);

% Chop-off Noise nad Remove GC Response
gc_power_spectrum_measurement = gc_power_spectrum_measurement(ind_1:ind_2);
wavelength_space_measurement = wavelength_space_measurement(ind_1:ind_2);
power_space_measuremet = power_space_measuremet(ind_1:ind_2);
power_spectrum_measurement = power_space_measuremet - gc_power_spectrum_measurement;
% Smooth out data
power_spectrum_sm_measurement = smooth(wavelength_space_measurement,power_spectrum_measurement,0.08,'loess');

% Open Interconnect Simulation and extract data
% file name esantacr_MZI_Unbalanced_590_f_interconnect.mat
file = strrep(file,"_f.mat","_f_interconnect.mat");
filename = fullfile(folder,file);

if ~exist(filename,"file")
    message = sprintf("%s does not exist", filename);
    uiwait(warndlg(message));
end

% Extract Data

data = load(filename);
wavelength_space_lum = 10^(-9).*flip(data.lum.x0); % Convert to [m]
power_space_lum = flip(data.lum.y0);

% Find indecies for these wavelength
A = abs(min_wavelength.*ones(1,length(wavelength_space_lum))-wavelength_space_lum);
B = abs(max_wavelength.*ones(1,length(wavelength_space_lum))-wavelength_space_lum);
[~, lum_ind_1] = min(A);
[~, lum_ind_2] = min(B);

gc_power_spectrum_lum = gc_power_spectrum_lum(lum_ind_1:lum_ind_2);
wavelength_space_lum = wavelength_space_lum(lum_ind_1:lum_ind_2);
power_space_lum = power_space_lum(lum_ind_1:lum_ind_2);
power_spectrum_lum = power_space_lum - gc_power_spectrum_lum;
% Smoouth out data
power_spectrum_sm_lum = smooth(wavelength_space_lum,power_spectrum_lum,0.08,'loess');

% Plot Signal
figure(1); clf;
title(sprintf("Original Optical Measurement of %s", file));
xlabel("Wavelength [m]");
ylabel("Power in [dBm]");
hold on;
plot(wavelength_space_measurement,power_space_measuremet);
hold off;

% Plot Calculated
figure(2); clf;
title(sprintf("Calculated Transfer Function of MZI given deltaL = %i um", deltaL));
xlabel("Wavelength [m]");
ylabel("Power in [dB]");
hold on;
plot(wavelength_space_measurement,T_MZI_dBm(deltaL,wavelength_space_measurement));
hold off;

% Plot Comparison between the calculated and measured transferfucntion
figure(3); clf;
title("Comparison Between Optical Measurement and Calculated Transfer Function");
xlabel("Wavelength [m]");
ylabel("Power");
hold on;
plot(wavelength_space_measurement,power_spectrum_measurement);
plot(wavelength_space_measurement,T_MZI_dBm(deltaL,wavelength_space_measurement));
plot(wavelength_space_lum,power_spectrum_lum);
legend('Smoothed Optical Measurement','Calculated Transfer Function','Lumertical Interconnect Sim','Location','SE');
hold off;