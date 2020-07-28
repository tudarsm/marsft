%% EXAMPLE USAGE
% This script shows the example usage of the MARSFT code for evaluating and
% simulating CARS spectra

% clean up
clear
close all
clc

tic;
fprintf('Simulating CARS spectrum...')
% use output h to get a handle on some subfunctions. currently only
% plotting, might be extended in the future
% use key-value type arguments to adjust the simulation parameters. check
% out the input parser of marsft_sim to see all options.
[s,h]=marsft_sim('T',2000,'LineWidthMultiplier',2);
fprintf('done in %.2f s.\n',toc);
h.plotSpec(s,1)

% you can also access library spectra. note, that this finds the nearest
% temperature and linewidth multiplier value available! i.e. if your
% library goes up to 2500 K, requesting a 3000 K spectrum will return the
% 2500 K spectrum. use with caution.
% generate the library first using example_generate_library.m
lib=load('library.mat');
fprintf('Retrieve spectrum from library...')
tic;
s2=marsft_sim('T',2000,'library',lib);
fprintf('done in %.2f s.\n',toc);
h.plotSpec(s2,2)
