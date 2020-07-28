% this is an example on how to fit multiple spectra in parallel

clear
close all
clc

% load the library
lib=load('library_1bar.mat');
%% some options
options = marsft_options(lib,'Plot','off','Method','mean','Display','summary');

%% distribute the library to the workers
% unfortunately, workers cannot share memory. so beware, your library will
% be copied in memory to every worker.
% this automatically starts a parallel pool with your default preferences.
% if you want to control it, start in manually prior to execution of this
% line.
Dlib=parallel.pool.Constant(lib);
%% simulate a spectrum
sim = marsft_sim('T',2000,'LineWidthMultiplier',3,'linewidth',1,'xN2',0.8);
waven = linspace(2270,2340,256);
signal = interp1(sim.wavenumberarray,sim.spectra.CARS,waven);
%% Compare runtime of sequential vs. parallel execution
% set random number seed to default for reproducibility
rng('default')
totalseq=tic;
for ii=1:12
    marsft(lib,waven,signal,options);
end
fprintf('Runtime for sequential fitting: %.2f s.\n',toc(totalseq));
% rng(s)

% parallel
% attention: the runtime values for the summary output represent the actual runtime of the spectrum per core. compare the
% total runtime!
totalpar=tic;
parfor ii=1:12
    %use Dlib.Value to tell the workers to use the already distributed
    %library instead of distributing it at runtime
    marsft(Dlib.Value,waven,signal,options);
end
fprintf('Runtime for parallel fitting: %.2f s.\n',toc(totalpar));
% rng(s)