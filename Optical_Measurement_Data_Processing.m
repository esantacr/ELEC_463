% Script to process data from optical measurements

% Author: Emmanuel Santacruz
% Date:         2023-12-10
% Last Updated: 2023-12-10

folder = 'C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\OPTICAL MEASUREMENTS\CHIP 4';

% File Naming Convention
% esantacr_MZI_(convention)_(number)_(f\b).mat
% EOPhi_(number) -> number is the length of the polling arm
% Balanced_(number) -> number is the length of the straight MZI Branch
% Poling_(number) -> number is the expected value of r33.
% Unbalanced_(number) -> number is the path length difference between the two branches
% (f\b) -> two iterations of the same measurement

% Test Names:
% esantacr_MZI_Balanced_25_f.mat
% esantacr_MZI_EOPhi_35_f.mat
% esantacr_MZI_Poling_8377_f.mat
% esantacr_MZI_Unbalanced_290_f.mat

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
wavelength_space = data.scandata.wavelength;
gc_power_spectrum = data.scandata.power(:,power_channel);


% Open file and extract data
file = "esantacr_MZI_Poling_8377_f.mat"; 
filename = fullfile(folder, file);

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
wavelength_space = data.scandata.wavelength;
power_space = data.scandata.power(:,power_channel);

% Determine the min index and max index for which the power is less the power threshold

[max_power,max_ind] = max(power_space);     % Find the maximum power and its index
N = length(wavelength_space);
window_size = N/10;                         % Points to the left and points to the right
ind_1 = 0;                                  % Lower index
ind_2 = 0;                                  % Upper index

% Determine the correct vlaues for the indecies to chop-off 
if max_ind > window_size
    ind_1 = max_ind - window_size;
else 
    ind_1 = 1;
end

if max_ind + window_size > N
    ind_2 = N;
else
    ind_2 = max_ind + window_size;
end

% Plot Signal
figure(1); clf;
title(sprintf("Original Optical Measurement of %s", file));
xlabel("Wavelength [m]");
ylabel("Power in [dBm]");
hold on;
plot(wavelength_space,power_space);
hold off;

% Chop-off Noise
gc_power_spectrum = gc_power_spectrum(ind_1:ind_2);
wavelength_space = wavelength_space(ind_1:ind_2);
power_space = power_space(ind_1:ind_2);
power_spectrum = power_space - gc_power_spectrum;
% Smooth out data
power_spectrum_sm = smooth(wavelength_space,power_spectrum,0.08,'loess');

% Plot Figures
figure(2); clf;
title(sprintf("Adjusted Optical Measurement of %s",file));
xlabel("Wavelength [m]");
ylabel("Power in [dBm]");
hold on;
plot(wavelength_space,power_spectrum);
hold off;

% Compare Smoothed Version of Data
figure(3); clf;
title(sprintf("Comparison of Post-Processed Data"));
xlabel("Wavelength [m]");
ylabel("Power is [dBm]");
hold on;
plot(wavelength_space,power_spectrum,'b.');
plot(wavelength_space,power_spectrum_sm,'r-');
hold off;