function options = marsft_options(lib,varargin)
% This function generates an option struct used by marsfit.

p = inputParser;

%% Fitting parameters for spectra
% Syntax: 'Key',[Starting value, lower bound, upper bound] for parameters
% that should be fitted. 'Key',fixed value for a parameter that should be
% kept fixed.
addParameter(p,'T',[1000,lib.Ts(1),lib.Ts(end)]); % Temperature in K
addParameter(p,'Mult',[floor(length(lib.linwidmults)/2), 1, length(lib.linwidmults)]); % Linewidth multiplier in units index!
addParameter(p,'xN2',[0.8,0.3,1.0]); % N2 Mole fraction
addParameter(p,'linewidth',[1.2,lib.preconvolution*1.05,2]); % instrumental linewidth
addParameter(p,'IntExpansion',[1,0.9,1.1]); % Intensity expansion. 
addParameter(p,'WavenumberShift',[0,-3,3]); % Wavenumbershift
addParameter(p,'WavenumberExpansion',[0,-5,5]); % Wavenumber expansion (tries to compensate for inaccuracies in spectral calibration)

%% Constants for spectral fitting
addParameter(p,'BufferGasSusc',8.5); % mole fraction averaged suscepbtibility of the buffer gas

%% METHOD
% 'mean','single-detailed','single-quick'.
% Mean just fits the mean spectrum
% single-detailed fits all spectra in the stack with the given set of parameters
% single-quick fits the mean spectrum and fixes all parameters except temperature for the single shot fits.
addParameter(p,'Method','single-quick');   

%% VARIANCE WEIGHTED FITTING
% supply the variance data. if it is set to 1, no weighting by the variance
% is performed
addParameter(p,'VarianceData',-1); % supply a custom variance vector (e.g. from stokes noise or a camera model or a combination)
% this selects the variance model that is used for weighting the least squares solution
% Note: if variancedata is supplied, this parameter is overriden!
addParameter(p,'VarianceModel','none'); % can be 'none','shot-nosie','std'. 

%% parameters for the genetic algorithm
addParameter(p,'PopulationSize',-1);        % Set a fixed population size
addParameter(p,'NormPopulationSize',25);    % Set the population size as a multiple of the number of free parameters
addParameter(p,'MinPopulationSize',50);    % Set the minimum population size
addParameter(p,'FunctionTolerance',1e-12);  % Function tolerance for GA
addParameter(p,'UseVectorized',true);     % Compute all children in one generation in a vectorized way
addParameter(p,'MaxStallGenerations',-1);   % Fixed number of MaxStallGenerations.
addParameter(p,'NormMaxStallGenerations',3);     % Set the MaxStallGenerations as a multiple of the number of free parameters
addParameter(p,'Display','summary');            % Display: 'off','iter'

%% additional parameters for runtime control
addParameter(p,'Plot','off');       % show plots.

parse(p,varargin{:});

%% store all arguments in options struct
for fn = fieldnames(p.Results)'
    options.(fn{:})=p.Results.(fn{:});
end

%% SANITY CHECKS
% if a parameter is set to be fixed, set the lower and upper bounds to the
% parameter (which causes GA to treat it as a constant)
if numel(options.T) == 1; options.T=repmat(options.T,1,3);end
if numel(options.Mult) == 1; options.Mult=repmat(options.Mult,1,3);end
if numel(options.xN2) == 1; options.xN2=repmat(options.xN2,1,3);end
if numel(options.linewidth) == 1; options.linewidth=repmat(options.linewidth,1,3);end
if numel(options.IntExpansion) == 1; options.IntExpansion=repmat(options.IntExpansion,1,3);end
if numel(options.WavenumberShift) == 1; options.WavenumberShift=repmat(options.WavenumberShift,1,3);end
if numel(options.WavenumberExpansion) == 1; options.WavenumberExpansion=repmat(options.WavenumberExpansion,1,3);end

% return some errors in case the fitparameters violate constraints from the
% library
if any(options.linewidth<=lib.preconvolution);error('Desired minimum instrumental linewidth is below preconvolution value %g in library!',lib.preconvolution);end
% add error for temperature limits outside library
if options.T(3)>lib.Ts(end);error('Upper temperature limit is out of range defined in library.');end
if options.T(2)<lib.Ts(1);error('Lower temperature limit is out of range defined in library.');end
% also for line width multiplier
if any(options.Mult<1);error('Linewidth multiplier is fitted by index in library. This cannot be smaller than 1.');end
if options.Mult(3)>length(lib.linwidmults);error('Upper limit of line width multiplier index exceeds precomputed range in library.');end
% check that the starting value is always within bounds
if options.T(1)<options.T(2) || options.T(1)>options.T(3); error('Temperature starting value out of bounds');end
if options.Mult(1)<options.Mult(2) || options.Mult(1)>options.Mult(3); error('Linewidth multiplier starting value out of bounds');end
if options.xN2(1)<options.xN2(2) || options.xN2(1)>options.xN2(3); error('N2 mole fraction starting value out of bounds');end
if options.linewidth(1)<options.linewidth(2) || options.linewidth(1)>options.linewidth(3); error('Linewidth starting value out of bounds');end
if options.IntExpansion(1)<options.IntExpansion(2) || options.IntExpansion(1)>options.IntExpansion(3); error('Intensity expansion starting value out of bounds');end
if options.WavenumberShift(1)<options.WavenumberShift(2) || options.WavenumberShift(1)>options.WavenumberShift(3); error('Wavenumber shift starting value out of bounds');end
if options.WavenumberExpansion(1)<options.WavenumberExpansion(2) || options.WavenumberExpansion(1)>options.WavenumberExpansion(3); error('Wavenumber expansion starting value out of bounds');end

%% POPULATE FITPAR MATRIX
% create the fitpar matrix based on the provided options
fitpar = [options.T;...             % temperature
    options.linewidth;...           % instrumental linewidth
    options.xN2;...                 % N2 mole fraction
    options.Mult;...                % linewidth multiplier
    options.IntExpansion;...        % intensity expansion
    options.WavenumberShift;...     % wavenumber shift
    options.WavenumberExpansion;];  % wavenumber expansion

options.fitpar = fitpar;


%% SET PARAMETERS THAT MIGHT DEPEND ON THE NUMBER OF FREE PARAMETERS
% get the number of free parameters
options.numfree = sum((fitpar(:,1)~=fitpar(:,2)) & (fitpar(:,1) ~= fitpar(:,3)));

% if no fixed population size is set, use the normalized population size
% and multiply by the number of free parameters
if options.PopulationSize == -1;options.PopulationSize = options.NormPopulationSize*options.numfree;end
% same for maxstallgenerations
if options.MaxStallGenerations == -1;options.MaxStallGenerations = options.NormMaxStallGenerations*options.numfree;end

if options.PopulationSize < options.MinPopulationSize; options.PopulationSize = options.MinPopulationSize;end;
end