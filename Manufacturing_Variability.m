% Script to perfrom the analysis of Active Data

% Author: Emmanuel Santacruz
% Student Number: 95127841
% Date:         2023-12-11
% Last Updated: 2023-12-11

% LOCAL FOLDER WHERE DATA IS KEPT
folder = "C:\Users\emman\OneDrive - The University Of British Columbia\UNIVERSITY\FOURTH YEAR\ELEC 463\PROJECT\MANUFACUTING VARIABILITY";

% Corner Analysis Performed
Width_Analysis  = [480 500 520];
Height_Analysis = [210 220 230];

W = length(Width_Analysis);
H = length(Height_Analysis);

Legend = cell(H*W,1);

for i = 1:H*W
    k = floor((i+2)/3);
    l = mod(i,3);
    if l == 0, l = 3; end
    Legend{i} = sprintf("(%i,%i)",Width_Analysis(k),Height_Analysis(l));
end
% Variables of interest
neff_bunch = zeros(10,W*H);
ng_bunch = zeros(10,W*H);
n_coeff = zeros(3,W*H);

% Important Constants
c = 299792458;  % [m/s]
lambda0 = 1.55; % [micro meters]
count = 1;

% Matlab equations
neff_eq = @(nx, lambda) (nx(1) + nx(2).*(lambda-lambda0) + nx(3).*(lambda-lambda0).^2);
alpha = 0.1;    % Propagation Loss Coeff [m^-1] per meter
beta = @(nx,lambda) ((2*pi*neff_eq(nx,lambda))./(lambda)) - 1i*alpha/2.*ones(1,length(lambda));
T_MZI = @(nx,lambda,deltaL) (0.25*abs(1+exp(-1i*beta(nx,lambda)*deltaL)).^2);
T_MZI_dBm = @(nx,wavelength,deltaL) 10.*log10(T_MZI(nx,wavelength,deltaL));

for i = 1:W
    %EXTRACT DATA
    width = Width_Analysis(1,i);
    for j = 1:H
        % NAMING CONVENTION
        height = Height_Analysis(1,j);
        file = sprintf("%i_%i.mat",width,height);
        filename = fullfile(folder,file);
        if ~exist(filename,"file")
            message = sprintf("%s does not exist", filename);
            uiwait(warndlg(message));
        end
        % LOAD DATA
        data = load(filename);
        frequency = data.f;
        lambdas = c./frequency;
        lambdas = lambdas * 1e6; % Convert to microns
        vg = data.vg;
        ng = c./vg;
        neff = data.neff;
        % SAVE DATA
        neff_bunch(:,count) = real(neff);
        ng_bunch(:,count) = real(ng);
        % FIND CURVE FIT
        X=[2.4 0 0]; 
        format long
        X = lsqcurvefit (neff_eq, X, lambdas, neff);
        n_coeff(:,count) = real(X);
        count = count+1;
    end
end

% Plot the results
% neff Variability
figure(1); clf;
title("Variability in Effective Index");
xlabel("Wavelength [micro meters]");
ylabel("Effective Index neff");
hold on;
for i = 1:W*H
    plot(lambdas,neff_bunch(:,i))
end
legend(Legend);
hold off;
% ng Variability
figure(2); clf;
title("Variability in Group Index");
xlabel("Wavelength [micro meters]");
ylabel("Group Index ng");
hold on;
for i = 1:W*H
    plot(lambdas,ng_bunch(:,i))
end
legend(Legend);
hold off;
lambdas2 = linspace(1.4,1.6,1000);
% Transfer Function Variability
dL = 4.50;    % Given in micro meters
figure(3); clf;
title("MZI Transfer Function Variability, dL = 4.5 micron, FSR = 125nm");
xlabel("Wavelength [micro m]");
ylabel("Power [dBm]");
hold on;
for i = 1:W*H
    plot(lambdas2,T_MZI_dBm(n_coeff(:,i),lambdas2,dL))
end
legend(Legend)
hold off;
