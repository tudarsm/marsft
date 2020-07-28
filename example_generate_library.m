% this is an example on how to generate different libraries

clear
close all
clc

%% use this to create a default library
lib=marsft_genlibrary(); 

%% you can pass options to adapt it to your situation
% in this case, we want a different pressure and a custom temperature range
% also, we have a minimum instrumental of 1 cm-1 FWHM
lib=marsft_genlibrary('P',10,'MinInstrumental',1,'Ts',280:320); 

%% in any case, you have to save the library manually
save('library.mat','-struct','lib')