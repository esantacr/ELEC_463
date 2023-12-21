% Script to perform the plotting for the POLLING Data

% Author: Emmanuel Santacruz
% Student Number: 95127851
% Date:         2023-12-11
% Last Updated: 2023-12-11

% LOCAL FOLDER WHERE DATA IS KEPT
folder = "C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\POLLING DATA";

% THERE ARE TWO FILES TO CHOOSE FROM:
% WEDNESDAYS SECTION: elec463poling.csv
% THURSDAYS SECTION: defbuffer1.csv

file = "esantacr_poling_processing.csv";
filename = fullfile(folder,file);

% CHECK IF FILE EXISTS

if ~exist(filename,"file")
    message=sprintf("%s does not exist",filename);
    uiwait(warndlg(message));
end

data = readmatrix(filename);
time = data(:,3);
current_1 = data(:,10);
temperature_1 = data(:,12);
current_2 = data(:,6);
temperature_2 = data(:,8);
current_3 = data(:,14);
voltage_3 = data(:,15);
current_4 = data(:,16);
voltage_4 = data(:,17);


% Plot Polling Curve - Wednesday's Section
figure(1); clf;
title("Poling Curve - Wednesday's Section");
xlabel("time [s]");
yyaxis left;
ylabel("Temperature [celcius]");
hold on;
plot(time,temperature_1);
yyaxis right;
ylabel("Current [A]");
plot(time,current_1);
hold off 

% Plot Polling Curve - Thursday's Section
figure(2); clf;
title("Poling Curve - Thursday's Section");
xlabel("time [s]");
yyaxis left;
ylabel("Temperature [celcius]");
hold on;
plot(time,temperature_2);
yyaxis right;
ylabel("Current [A]");
plot(time,current_2);
hold off

% Plot I-V Curve - Wednesday's Section
figure(3); clf;
title("Current Voltage Relationship - Thursdays's Section");
xlabel("Current");
ylabel("Voltage [V]")
hold on;
plot(current_3,voltage_3);
hold off;

% Plot I-V Curve - Thursdays's Section
figure(4); clf;
title("Current Voltage Relationship - Wednesday's Section");
xlabel("Current");
ylabel("Voltage [V]")
hold on;
plot(current_4,voltage_4);
hold off;


