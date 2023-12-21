% Use a a linearized model for the change in n of an electro optic mod.
% Note r33 is assumed to be the 40 pm/V calculated in the lab

% Voltages Tested
Voltages_Tested = [-100 -90 -80 -70 -60 -4 -3 -2 -1 0 1 2 3 4 40 50 60 70 80 90 100];
% Number of Voltaegs Tested
N = length(Voltages_Tested);
% Create important variables
neff_simulated = zeros(1,N);
% Path to folder
folder = 'C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\EO_NEFF_SIMULATIONS';

for i = 1:N
    filename = fullfile(folder,sprintf("%iV.mat",Voltages_Tested(i)));
    if ~exist(filename,'file')
        message = sprintf('%s does not exist', filename);
        uiwait(warndlg(message));
    end
    data = load(filename);
    neff = real(data.neff);
    neff_simulated(i) = neff;
end

% Plot the results
figure(1); clf;
title("Modulation of neff with EO active polymer");
xlabel("Modulation Voltage [V]");
ylabel("neff");
hold on;
plot(Voltages_Tested,neff_simulated);
hold off;
