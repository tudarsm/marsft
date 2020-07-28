function library = marsft_genlibrary(varargin)
% This function generates a library of the real and imaginary parts of the
% spectra. These depend only on pressure, temperature and the linewidth
% multiplier. Pressure is treated as a constant parameter but could in
% principle be tabulated as another dimension (e.g. for application with
% varying pressure such as IC engines).


p = inputParser;

% standard parameters for marsft_sim
addParameter(p,'type','theosusc');      % use theoretical susceptibility output of marsft_sim
addParameter(p,'P',1);                  % Pressure in atm
addParameter(p,'Model','I');            % Which model? So far isolated lines ('I')
addParameter(p,'phi',0);                % Polarization angle in degress
addParameter(p,'theta',0);              % Angle between pump1 and probe in degrees
addParameter(p,'ROI',[2200 2390]);      % Boundaries of plot/fit range
addParameter(p,'HWFactors','JK');       % Selected Hermann-Wallis-Factors

% additional parameters for library
addParameter(p,'Ts',280:2500);              % Range of temperatures for tabulation. Make sure this is within the original range as tabulated by marsft_tabulation
addParameter(p,'linwidmults',2:0.04:4);      % Array containing all linewidth multipliers
addParameter(p,'MinInstrumental',0.8);      % Minimum allowed instrumental function
addParameter(p,'Wavenumberresolution',-1);  % Resolution of the resulting wavenumberarray. Use -1 to use the default (MinInstrumental / 5)
parse(p,varargin{:});

% store arguments in struct for later use
for fn = fieldnames(p.Results)'
    s.(fn{:})=p.Results.(fn{:});
end

% set the default wavenumberresolution if not specified otherwise
if s.Wavenumberresolution == -1
    s.Wavenumberresolution = s.MinInstrumental / 10;
end

% generate the resulting wavenumberarray
wavenumberarray = s.ROI(1):s.Wavenumberresolution:s.ROI(2);

% preallocate the arrays for the real and imaginary component of the
% theoretical susceptibility as single for memory efficiency
chireal = zeros(length(s.Ts),length(s.linwidmults),length(wavenumberarray),'single');
chires2 = zeros(length(s.Ts),length(s.linwidmults),length(wavenumberarray),'single');
% for every temperature and linewidth multiplier, generate an entry in the
% library.
fprintf('Creating a database with %d temperature and %d linewidth multiplier entries...\nGo grab a coffee, that is %.2f MB of spectra...\n',length(s.Ts),length(s.linwidmults),length(chireal(:))*2*4/1024/1024);
strcr = '';
indicator={'Z','M','A','R','S','F','T',' ','R','U','L','E'};
startgen=tic;
toctot=NaN(10,1);
for tidx = 1:length(s.Ts)
    ticset=tic;
    if tidx == 1
        eta = NaN;
        tocset = NaN;
    else
        meantime = mean(toctot(~isnan(toctot)));
        eta = (length(s.Ts)-tidx+1)*meantime/60;
    end
    strout = sprintf('Processing set %d of %d. Last set took %.2f s, estimated time left: %.2f min. %s...',tidx,length(s.Ts),tocset,eta,indicator{mod(tidx,length(indicator))+1});
    fprintf([strcr strout]);
    strcr = repmat('\b',1,length(strout));
    parfor lmidx = 1:length(s.linwidmults)
        % simulate the real and imaginary part of the spectrum
        scurr = marsft_sim('T',s.Ts(tidx),'linewidth',s.MinInstrumental,'LineWidthMultiplier',s.linwidmults(lmidx),'type',s.type,'Model',s.Model,'P',s.P,'phi',s.phi,'theta',s.theta,'ROI',s.ROI,'HWFactors',s.HWFactors);
        % generate the gaussian for preconvolution
        convrange = -2.5*s.MinInstrumental:scurr.res_wavenumber:2.5*s.MinInstrumental; % range for convolution kernels
        kernel = gaussian(convrange,0,s.MinInstrumental);
        % convolve the real (chirl) and resonant (chirl^2 + chiim^2) part with the kernel, interpolate
        % on output array and store in library
        chireal(tidx,lmidx,:)=interp1(scurr.wavenumberarray,conv(scurr.chi_real,kernel,'same'),wavenumberarray);
        chires2(tidx,lmidx,:)=interp1(scurr.wavenumberarray,conv(scurr.chi_res,kernel,'same'),wavenumberarray);
    end
    tocset=toc(ticset);
    toctot(mod(tidx,10)+1)=tocset;
end

% rotate the database such that the wavenumber direction is in the first
% dimenstion. this ensures that the spectral datapoints are consecutive in
% memory which makes accessing it much(!) faster.
chireal=shiftdim(chireal,2);
chires2=shiftdim(chires2,2);

% organize data in struct
library.wavenumberarray = wavenumberarray;
library.chireal = chireal;
library.chires2 = chires2;
% store the temperature and linewidth variables in the library
library.Ts = s.Ts;
library.linwidmults = s.linwidmults;
% store some of the conditions in the library,too
library.preconvolution = s.MinInstrumental;
library.res = s.Wavenumberresolution;
library.P = s.P;
% also store polarization and molecular information in library
library.tnres = sind(s.theta)*sind(s.phi)/3+cosd(s.theta)*cosd(s.phi);
% run simulation again just to get the molacular parameters. scurr is unavailable
% outside parfor loop, so just call it again. settings do not affet the
% molecular parameters, so just use defaults.
scurr = marsft_sim();
library.mol = scurr.mol;
fprintf('\nLibrary generated in %.2f min.\n',toc(startgen)/60);
end