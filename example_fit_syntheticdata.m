clear
close all
clc

% this is an example on how to use marsft for fitting
% run example_generate_library.m with the default library first or
% uncomment the following section

% lib=marsft_genlibrary();
% 
% % save the library for future use
% save('library.mat','-struct','lib')

% load a precomputed library
lib=load('library.mat');

%% simulate a spectrum
T = 1900;
Mult = 3;
xN2 = 0.8;
LineWidth = 1.2;

[sim,h] = marsft_sim('T',T,'LineWidthMultiplier',Mult,'linewidth',LineWidth,'xN2',xN2);
h.plotSpec(sim,1)

% make it look like an experiment, ie. create experimental grid and add
% some noise. create multiple shots
numshots = 100;
experiment.wavenumberarray = linspace(2270,2350,256); % typicall range for n2 ro-vib cars
experiment.signal = interp1(sim.wavenumberarray,sim.spectra.CARS,experiment.wavenumberarray); % interpolate on the experimental grid
experiment.signal = experiment.signal.*(1+(rand(numshots,length(experiment.signal))-0.5)/5); % add some multiplicative noise (simulates 10% stokes noise)
experiment.signal = experiment.signal+((rand(numshots,length(experiment.signal))-0.5)/2); % add some additive noise
plot(experiment.wavenumberarray,experiment.signal(:,:),'r-','LineWidth',1,'Color',[0.6 0.6 0.6]) % plot it
plot(experiment.wavenumberarray,experiment.signal(1,:),'r-','LineWidth',2) % plot it


%% Now fit it the spectra using the single-quick mode
% you can turn live plotting on or off
plotting='off';
options = marsft_options(lib,'VarianceModel','none','Plot',plotting,'Method','single-quick');
f1 = marsft(lib,experiment.wavenumberarray,experiment.signal,options);
%% Now fit it again, this time with one of the variance models
options = marsft_options(lib,'VarianceModel','std','Plot',plotting,'Method','single-quick');
f2 = marsft(lib,experiment.wavenumberarray,experiment.signal,options);
%% Now fit it again, but set some degrees of freedom fixed
% This is a little more representative, because some parameters only depend
% on the experimental setup or at least change only at a daily basis
% Fix the Wavenumber Shift and Expansion by supplying single values. These
% will be automatically used as initial condition as well as the lower and
% upper bound. This effectively fixes the value and allows the population
% size to be reduced, see marsft_options for more instructions.
% If you want to limit the range, please supply an array like this:
% [initial condition, lower bound, upper bound]. See linewidth:
options = marsft_options(lib,'VarianceModel','std','Plot',plotting,'Method','single-quick',...
    'WavenumberShift',0,'WavenumberExpansion',0,...
    'linewidth',[1,0.81,2]); % note, that linewidth cannot be <= preconvolution line width. if selected too small, this will throw an error.
% fit the spectra using the new options
% note how the result is equally good (or bad), notice lower runtime for
% mean fit due to reduced number of degrees of freedom
f3 = marsft(lib,experiment.wavenumberarray,experiment.signal,options);
%% Compare the results
figure;
hold all;
plot([f1.T]);
plot([f2.T]);
plot([f3.T]);
xlabel('Sample')
ylabel('Temperature in K')
box on
grid on
title(sprintf('T1: %.2f +- %.2f K, T2: %.2f +- %.2f K, T3: %.2f +- %.2f K',mean([f1.T]),std([f1.T]),mean([f2.T]),std([f2.T]),mean([f3.T]),std([f3.T])))
