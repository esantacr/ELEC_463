% User provides a matrix of neff values vs. wavelength

% Matlab curve fits to an expression.

 
load('Lumerical_Simulation.mat');

 

neff = real(neff)  % take the real part of the effective index.

 

c=299792458;  % speed of light, m/s

lambdas = c ./ f;  % f is the matrix of frequency points, 

                   % where the effective index is recorded.

lambdas = lambdas * 1e6  % convert to microns.

lambda0 = 1.55;   % replace with desired centre wavelength

 

figure; plot (lambdas, neff,'o','MarkerSize',10); hold on;

 

% use Matlab anonymous function for the effective index expression:

neff_eq = @(nx, lambda) (nx(1) + nx(2).*(lambda-lambda0) + nx(3).*(lambda-lambda0).^2); 

 

% initial guess.

X=[2.4 0 0]; 

 

plot ( lambdas, neff_eq(X, lambdas), 'r')

 

% curve fit to find expression for neff.

format long

X = lsqcurvefit (neff_eq, X, lambdas, neff)

 

r=corrcoef(neff,neff_eq(X, lambdas));

r2=r(1,2).^2;

disp (['Goodness of fit, r^2 value: ' num2str(r2) ])

 

lambdas2=linspace(min(lambdas), max(lambdas), 100);

 

plot ( lambdas2, neff_eq(X, lambdas2), 'k')

title(' Neff vs Wavelength - Non-Modulated Simulation')

xlabel ('Wavelength [nm]');

ylabel ('Effective Index');

 

legend ('Data','Initial Guess','Curve Fit')