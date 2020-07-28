function [wavenumberstruct] = generateWavenumberArray(dispersioncalibration,carscalibration,varargin)
%GENERATEWAVENUMBERARRAY 
%   Inputs:
%   dispersioncalibration,carscalibration: Nx1 vector (mean spectrum)
%   wavenumberrange: (2x1) vector in the form [2270 2340] (default value).
%
%   This function generates the wavenumberarray obtained by the dispersion
%   calibration with the Xe-PenRay lamp. It offers two features:
%   1. only a dispersioncalibration:
%   this is useful when measurements are
%   performed with the SPEX spectrometer because the spectral resolution is
%   sufficiently low to accomodate both lines used for calibration (473.415
%   nm and 479.261 nm) AND the entire CARS signal without moving the
%   grating.
%
%   2. dispersioncalibration + carscalibration:
%   When using the McPhearson spectrometer, the spectral resolution is a
%   little better. Unfortunately this means, that the two spectral lines
%   for calibration are very close to the edges of the camera and already a
%   little vignetted. To prevent vignetting of the CARS signal, which is
%   very close to the 'left' calibration line, for actual measurements, the
%   grating has to be moved a little to center the CARS signal. Thus, two
%   calibrations are necessary, one for the dispersion and one for the
%   position of the CARS signal.
%
%   If you want to go for option 1, just pass [] to the second argument of
%   this function.
%   
%   Optional key-value pair arguments:
%   Wavenumberrange: Not a required input parameter and defaults
%   to [2270 2340], which accomodates the first two vibrational bands of
%   the CARS signal. Change this only if you know what you're doing....
%   Might be useful if you only expect cold signals, where the relative
%   population of the hot band is negligible.
%
%   Showplots: shows explanatory plots if set to 1
%
%   Probewavelength: defaults to 532 nm. Is added as a parameter for adding
%   the possibility to use different wavelength. Highly untested! Use with
%   caution. This affects the conversion from wavelength to raman shift in
%   wavenumbers.
%
%
%   The output of the function is a struct containing the wavenumberarray
%   (used for generateCARSLibrary) and the region of interest, defined by the first and last
%   selected wavenumber in units pixel. The latter is used in
%   preprocessCARSSpectra to crop the spectrum to the desired region of
%   interest.
%

p = inputParser;
addParameter(p,'Wavenumberrange',[2270 2340]);
addParameter(p,'Probewavelength',532);
addParameter(p,'Showplots',0,@isnumeric);    % show explanatory plots
addParameter(p,'ManualLineInput',-1,@isnumeric);
parse(p,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of the dispersion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Line positions
% if the linepositions are not manually given, assume that the dispersion
% calibration contains two lines at the following locations
if p.Results.ManualLineInput ~= -1
    % check that input is valid
    if size(p.Results.ManualLineInput,2)~=2
        error('ManualLineInput expects a Nx2 matrix, where the first row is the line position in px and the second row the corresponding wavelength as retrieved from the spectral database of your calibration lamp.');
    end
    % atm, just take the first and the last line and compute linear
    % dispersion. Later this could be improved by using a polynomial fit
    % over multiple lines.
	line1 = p.Results.ManualLineInput(2,1);
    line2 = p.Results.ManualLineInput(2,end);
    peak1 = p.Results.ManualLineInput(1,1);
    peak2 = p.Results.ManualLineInput(1,end);
    dispersion = (line2-line1)/(peak2-peak1);
else
    % 1. find pixel positions of the two peaks at 473.415 nm and 479.261 nm.
    % look for the maximum in the left and the right half of the spectrum
    [~,peak1] = max(dispersioncalibration(1:floor(length(dispersioncalibration)/2)));
    [~,peak2] = max(dispersioncalibration(floor(length(dispersioncalibration)/2):end));
    peak2 = peak2+floor(length(dispersioncalibration)/2)-1; %shift by half the width, because we started searching at the beginning of the right half
    line1 = 473.415; % left line wavelength in nm
    line2 = 479.261; % right line wavelength in nm
    % compute the dispersion in nanometer per px
    dispersion = (line2-line1)/(peak2-peak1);
end

% generate a linearly spaced vector that uses the prior computed dispersion
% as gradient and has the correct values at peak1 and peak2
t = -dispersion*(peak1-1)+line1;
x = 0:1:size(dispersioncalibration,1)-1;
wavelengtharray_dispersion = dispersion*x+t;
% wavelengtharray_dispersion=linspace(wavelength_leftmostpixel,wavelength_rightmostpixel,length(dispersioncalibration));

% check if carscalibration was supplied or not
% if not, just use the wavelengtharray based on the dispersion alone. If
% yes, shift the wavelengthvector accordingly
if isempty(carscalibration)
    wavelengtharray_cars=wavelengtharray_dispersion;
    carscalibration = dispersioncalibration;
else
    % find the peak in the carscalibration. It should only be one!
    [~,peak_cars]=max(carscalibration);
    % difference in peak location
    shift_px = peak_cars-peak1;     % unit px
    shift_nm = shift_px*dispersion; % unit nm
    % shift the array to the left
    wavelengtharray_cars = wavelengtharray_dispersion-shift_nm;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of wavenumberrange %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalwavenumberarray = 1 ./ (wavelengtharray_cars *0.000000001)/100 - 1/(p.Results.Probewavelength*0.000000001)/100; % in cm^-1

% crop to region of interest
rightend = find(totalwavenumberarray<p.Results.Wavenumberrange(1),1,'first');
leftend = find(totalwavenumberarray<p.Results.Wavenumberrange(2),1,'first');
roi = leftend+1:rightend-1;
wavenumberarray = totalwavenumberarray(roi);

%%%%%%%%%%
% OUTPUT %
%%%%%%%%%%
% the input variables are also stored in the struct for later traceability
% of what was used for this evaluation.
wavenumberstruct.wavelengtharray_dispersion=wavelengtharray_dispersion;
wavenumberstruct.wavelengtharray_cars=wavelengtharray_cars;
wavenumberstruct.dispersioncalibration=dispersioncalibration;
wavenumberstruct.carscalibration=carscalibration;
wavenumberstruct.wavenumberrange=p.Results.Wavenumberrange;
wavenumberstruct.wavenumberarray = wavenumberarray;
wavenumberstruct.roi = roi;

% show the effect of the calibration and/or shifting
if p.Results.Showplots
figure;
hold all;
plot(dispersioncalibration,'DisplayName','Dispersion');
plot(carscalibration,'DisplayName','CARS');
title('Calibration spectra for spetrometer in Pixels');
ylabel('Signal in a.u.');
xlabel('Spectral direction in px');
legend('show')
box on;grid on;

figure;
hold all;
plot(wavelengtharray_dispersion, dispersioncalibration,'DisplayName','Dispersion');
plot(wavelengtharray_dispersion, carscalibration,'DisplayName','CARS');
title('Calibration spectra for spetrometer in nm, unshifted');
ylabel('Signal in a.u.');
xlabel('Spectral direction in nm');
legend('show')
box on;grid on;

figure;
hold all;
plot(wavelengtharray_dispersion, dispersioncalibration,'DisplayName','Dispersion');
plot(wavelengtharray_cars, carscalibration,'DisplayName','CARS');
title('Calibration spectra for spetrometer in nm, shifted. Left maxima should overlap.');
ylabel('Signal in a.u.');
xlabel('Spectral direction in nm');
legend('show')
box on;grid on;

figure;
hold all;
plot(totalwavenumberarray,carscalibration);
plot([p.Results.Wavenumberrange(1),p.Results.Wavenumberrange(1)],[-1000 65000],'k--');
plot([p.Results.Wavenumberrange(2),p.Results.Wavenumberrange(2)],[-1000 65000],'k--');
title('Spectral calibration over Raman-Shift');
ylabel('Signal in a.u.');
xlabel('Raman-Shift in cm^-^1');
ylim([min(carscalibration) max(carscalibration)])
legend('CARS Calibration','ROI')
box on;grid on;
end

end

