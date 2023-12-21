% Script to perfrom the analysis of Active Data

% Author: Emmanuel Santacruz
% Student Number: 95127841
% Date:         2023-12-11
% Last Updated: 2023-12-11

% LOCAL FOLDER WHERE DATA IS KEPT
folder = "C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\ACTIVE DATA\3rd try";

% Many Voltages were tested here is a list of them
Voltages_Tested = [-30 -20 -14 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 9 10 20 30 40 100];
num_tests = length(Voltages_Tested);            % Number of tests performed

% Variables of Interest
Max_Wavelength = zeros(1,num_tests);
Max_Power = zeros(1,num_tests);
Min_Power = zeros(1,num_tests);
delta_power = zeros(1,num_tests);

for i = 1:num_tests
    % NAMING CONVENTION
    file = sprintf("%iV.mat",Voltages_Tested(i));
    filename = fullfile(folder,file);
    % CHECK IF FILE EXISTS
    if ~exist(filename,"file")
        message = sprintf("%s does not exist", filename);
        uiwait(warndlg(message));
    end
    % LOAD DATA
    data = load(filename);
    wavelength_spectrum = data.wavelength;
    % SAVE GC POWER SPECTRUM
    if i == 1
        gc_power_spectrum = data.power(:,1);
        gc_wavelength_spectrum = data.wavelength;
    end
    % SIZE OF POWER ARRAY
    [~,columns] = size(data.power);
    % CHOOSE CORRECT CHANNEL
    if columns == 1
        power_spectrum = data.power(:,1);
    else
        power_spectrum = data.power(:,2); 
    end
    % FIND MAXIMUMS AND MINIMUMS
    [max_power,max_power_ind] = max(power_spectrum);
    [min_power,min_mower_ind] = min(power_spectrum);
    % STORE VALUES
    Max_Power(i) = max_power;
    Min_Power(i) = min_power;
    Max_Wavelength(i) = wavelength_spectrum(max_power_ind);
    delta_power(i) = min_power - max_power;
    % FOR FUTURE ANALYSIS
    if Voltages_Tested(i) == -4
        power_spectrum_neg_4 = power_spectrum;
    elseif Voltages_Tested(i) == 0
        power_spectrum_0 = power_spectrum;
    elseif Voltages_Tested(i) == 4
        power_spectrum_4 = power_spectrum;
    end
end

% PLOT MAX POWER AND WAVELENGTH VS VOTAGE
figure(1); clf;
title("Maximum Power and Wavelength Measured as a function of Modulated Voltage");
xlabel("Modulation Voltage [V]");
yyaxis right;
ylabel("Maximum Power [dBm]");
hold on;
plot(Voltages_Tested,Max_Power);
yyaxis left;
ylabel("Maximum Wavelength [m]");
plot(Voltages_Tested,Max_Wavelength);
hold off;

% PLOT ATTENUATION VS VOLTAGE
figure(2); clf;
title("Change in Power Spectrum Amplitude as a function of Modulated Voltage");
xlabel("Modulation Voltage [V]");
ylabel("del-Power [dBm]");
hold on;
plot(Voltages_Tested,delta_power);
hold off;

% PLOT MIN AND MAX POWER VS VOLTAGE
figure(3); clf;
title("Maximum Power and Minimum Power as function of Modulated Voltage");
xlabel("Modulation Voltage [V]");
yyaxis right;
ylabel("Maximum Power [dBm]");
hold on;
plot(Voltages_Tested,Max_Power);
yyaxis left;
ylabel("Minimum Power [dBm]");
plot(Voltages_Tested,Min_Power);
hold off;

% EXTRACT THE DATA FOR -4V, 0V and 4V TO CREATE A COMPARISON
[~,ref_ind] = max(gc_power_spectrum);           % GC Max Wavelength
window_size = ceil(length(gc_wavelength_spectrum)/20);   % Size of Window To Be Plotted
num_points = length(gc_wavelength_spectrum);
% FIND INDECES
ind_1 = 1;
ind_2 = 1;
% LOW INDEX
if ref_ind - window_size <= 1
    ind_1 = 1;
else
    ind_1 = ref_ind - window_size;
end
% HIGH INDEX
if ref_ind + window_size >= num_points
    ind_2 = num_points;
else
    ind_2 = ref_ind + window_size;
end
% CHOP OFF UNNECESSARY PARTS
wavelength_spectrum = gc_wavelength_spectrum(ind_1:ind_2);
gc_power_spectrum = gc_power_spectrum(ind_1:ind_2);
power_spectrum_neg_4 = power_spectrum_neg_4(ind_1:ind_2);
power_spectrum_0 = power_spectrum_0(ind_1:ind_2);
power_spectrum_4 = power_spectrum_4(ind_1:ind_2);
% SMOOTH SPECTRUM
gc_power_spectrum = smooth(wavelength_spectrum,gc_power_spectrum,0.08,'loess');
power_spectrum_neg_4 = smooth(wavelength_spectrum,power_spectrum_neg_4,0.08,'loess');
power_spectrum_0 = smooth(wavelength_spectrum,power_spectrum_0,0.08,'loess');
power_spectrum_4 = smooth(wavelength_spectrum,power_spectrum_4,0.08,'loess');
% REMOVE GC RESPONSE
power_spectrum_neg_4 = power_spectrum_neg_4 - gc_power_spectrum;
power_spectrum_0 = power_spectrum_0 - gc_power_spectrum;
power_spectrum_4 = power_spectrum_4 - gc_power_spectrum;

% FIND SHIFT IN WAVELENGTH
[~, max_power_index_neg_4] = max(power_spectrum_neg_4);
[~, max_power_index_0] = max(power_spectrum_0);
wavelength_shift = abs(wavelength_spectrum(max_power_index_0) - wavelength_spectrum(max_power_index_neg_4));

% PLOT TRANSFER FUNCTION SIDE BY SIDE
dim = [0.15 0 0.3 0.3];
shift_str = "Shift ~ 530 pm/4V or 132.5 pm/V  ";
fsr_str = "FSR ~ 1.4 nm  VpiL ~ 5.28 V";
str = {strcat(shift_str,fsr_str)};
figure(4); clf;
title("Plot of Modulated MZI at voltages: (-4,0,4)V");
xlabel("Wavelength [m]");
ylabel("Power [dBm]");
hold on;
plot(wavelength_spectrum,power_spectrum_neg_4);
plot(wavelength_spectrum,power_spectrum_0);
plot(wavelength_spectrum,power_spectrum_4);
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off;
legend('-4V','0V','4V','Location',"NE");
