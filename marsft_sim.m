function [s,h] = marsft_sim(varargin)

% h is an optional output to subfunction handles
h.plotSpec = @plotSpec;

p = inputParser;
addParameter(p,'type','cars');          % cars, theosusc or library
addParameter(p,'T',2000);               % Temperature
addParameter(p,'P',1);                  % Pressure in atm!
addParameter(p,'linewidth',sqrt(2)*0.75);         % instrument function
addParameter(p,'xN2',0.8);              % N2 mole fraction
addParameter(p,'chinr_buffergas',8.5);            % Buffer gas non resonant susceptibility
addParameter(p,'LineWidthMultiplier',1);% Linewidth multiplier
addParameter(p,'Model','I');            % Which model? So far isolated lines ('I')
addParameter(p,'HWFactors','JK');       % Which model for the Hermann-Wallis-Factors? JK, LBY, TB
addParameter(p,'library',-1);           % library for accessing spectra
addParameter(p,'phi',0);                % Polarization angle in degress
addParameter(p,'theta',0);              % Angle between pump1 and probe in degrees
addParameter(p,'ROI',[2200 2390]);      % Boundaries of plot/fit range
addParameter(p,'ROIExp',50);            % Expand the ROI internally to avoid wing effects by this amount
parse(p,varargin{:});

% store arguments in struct for later use
for fn = fieldnames(p.Results)'
    s.(fn{:})=p.Results.(fn{:});
end

% if type library, return a spectrum from the library
if isstruct(p.Results.library)
    library=p.Results.library;
    % get the total number density
    totnum = sconst('NA')*273.15/p.Results.T*(p.Results.P/sconst('molarvolume'));
    % get closest spectrum
    [~,idxT] = min(abs(library.Ts-p.Results.T));
    [~,idxLW] = min(abs(library.linwidmults-p.Results.LineWidthMultiplier));
    chireal = library.chireal(:,idxT,idxLW);
    chires2 = library.chires2(:,idxT,idxLW);
    
    % return actual temperature and linewidth multiplier in struct
    s.T = library.Ts(idxT);
    s.LineWidthMultiplier = library.linwidmults(idxLW);
    
    % so far, only N2 implemented. everything that is not N2 contributes
    % with the same non resonant background.
    % calculate non-resonant contribution from N2, account for polarization
    chi_nr_r = library.tnres*library.mol.N2.CHINR*totnum*p.Results.xN2;
    
    % add the buffer gas nr contribution
    xBuffer = 1-p.Results.xN2;
    % attention: the unit of the constant is in cm3/(erg amagat), convert
    % to susceptibility per molecule using Avogadro and molar volume!
    chi_nr = chi_nr_r + library.tnres*totnum*p.Results.chinr_buffergas/(sconst('NA')/sconst('molarvolume'))*xBuffer;
    
    % convolve the resulting spectrum with the additional desired gaussian
    addgausswidth = sqrt(p.Results.linewidth^2-library.preconvolution^2);
    
    % create the kernel
    convrange = -10:library.res:10;
    kernel = gaussian(convrange,0,addgausswidth);
    
    
    % compute the spectrum
    % this is basically equation (6) in the paper
    s.spectra.CARS = sqrt(conv(chires2*p.Results.xN2^2 + 2*chi_nr*chireal*p.Results.xN2 + chi_nr^2,kernel,'same'));
    s.wavenumberarray = library.wavenumberarray;
    return
end

% Js and Vs to include. Change only if you know what you're doing.
s.J = 0:100;
s.V = 0:6;

% molecular params
s = getmol(s);
% term energies
s = calculateTermEnergies(s);
% calculate line positions
s = calculateLinePositions(s);
% and now the things that are temperature dependent
% calculate line widths in units cm-1/bar
s = calculateLineWidths(s);
% boltzmann fraction
s = calculateBoltzmannFraction(s);
% calculate the peak suscepbtibility chiamp
s=chiamplitudes(s);
% limit the transitions to the desired range
s=selectTransitions(s);
% create the theoretical suscebitility spectrum
s = createTheoSusc(s);
% if desired, create the convolved spectrum
if strcmp(s.type,'cars')
    s = createCARSspectrum(s);
end

end


function s = getmol(s)
% spectroscopic constants for a given species. See section I-C-4 in carsdoc.txt for some
% explanation what those parameters are. The parameters are stored here
% directly to speed things up and to have everything in one file
% Naming of parameters are kept the same way as it is in CARSFT for easier
% orientation.
% If other (diatomic) species are to be included, just add them here.
N2.WE       =  0.235854024E+04;
N2.WX       =  0.14305770E+02;
N2.WY       = -0.50668000E-02;
N2.WZ       = -0.10950000E-03;
N2.BE       =  0.19982600E+01;
N2.ALPHAE   =  0.17303500E-01;
N2.DE       =  0.57740000E-05;
N2.BETAE    =  0.15500000E-07;
N2.GAME     = -0.31536099E-04;
N2.DELTE    =  0;
N2.H0       =  0.30000000E-11;
N2.HE       =  0.18000000E-11;
N2.RE       =  0.20743101E+01;
N2.GAM      =  0.47640000E+01;
N2.DGAMDR   =  0.74000000E+01;
N2.GNE      =  0.60000000E+01;
N2.GNO      =  0.30000000E+01;
N2.DFSPHI   =  0;
N2.AG       =  0.11660000E+01;
N2.CHINR    =  0.85000000E+01;
N2.AC1      =  0.24158000E-11;
N2.CMASS    =  0.28013000E+02;

% store the molecular parameters in the struct
s.mol.N2=N2;

% store a cell array of included species
% for later versions, this can be used to include more species
% this cell array is iterated in every species dependent subfunction.
s.species={'N2'};
end

function s = calculateTermEnergies(s)
% Compute the term values for all species.

% the following works for all diatomic molecules with provided molecular
% data
% iterate through all species in s.species
for specimen = s.species
    % CHINR is in units 10^18 cm^3/(erg*amagat), convert CHINR per
    % molecule at standard conditions (101325 Pa, 273.15 K)
    s.mol.(specimen{:}).CHINR = s.mol.(specimen{:}).CHINR / (sconst('NA')/sconst('molarvolume'));
    
    
    % Anharmonic oscillator, see following book, which should be available
    % on google books. The equations are right below figure 47.
    %     @book{herzberg2013molecular,
    %   title={Molecular spectra and molecular structure},
    %   author={Herzberg, Gerhard},
    %   volume={1},
    %   year={2013},
    %   publisher={Read Books Ltd}
    % }
    
    % kept upper case W for omega for consistency, this is equation
    % (III,77) in herzberg
    s.mol.(specimen{:}).W0 = s.mol.(specimen{:}).WE-s.mol.(specimen{:}).WX+3/4*s.mol.(specimen{:}).WY;
    s.mol.(specimen{:}).W0X0 = s.mol.(specimen{:}).WX-3/2*s.mol.(specimen{:}).WY;
    s.mol.(specimen{:}).W0Y0 = s.mol.(specimen{:}).WY;
    
    % set the nuclear spin as degeneracy for even and odd Js
    % ATTENTION: J of course starts counting at zero, matlab at one
    % so 1,3,5 ... corresponds to the EVEN Js (0,2,4,...)
    s.mol.(specimen{:}).GJ(1:2:length(s.J))=s.mol.(specimen{:}).GNE; % EVEN
    s.mol.(specimen{:}).GJ(2:2:length(s.J))=s.mol.(specimen{:}).GNO; % ODD
    
    % this is equation (III,76)
    s.mol.(specimen{:}).G0V = s.mol.(specimen{:}).W0.*(s.V) - s.mol.(specimen{:}).W0X0*(s.V).^2 + s.mol.(specimen{:}).W0Y0*(s.V).^3;
    
    % Vibrating rotator
    % See Herzberg, see equations (III,124-126)
    % to make it a little easier to read, some shortcuts to (v+1/2) and
    % J(J+1). Also, this is vectorized, make sure that v and J are in
    % different dimensions.
    % again, keep the uppercase lettering for easier readability.
    VPH = (s.V)+0.5;
    JJ = ((s.J).*(s.J+1))';
    Bv = s.mol.(specimen{:}).BE+VPH.*(-s.mol.(specimen{:}).ALPHAE+VPH.*s.mol.(specimen{:}).GAME);
    Dv = s.mol.(specimen{:}).DE+VPH.*(s.mol.(specimen{:}).BETAE+VPH.*s.mol.(specimen{:}).DELTE);
    Hv = s.mol.(specimen{:}).H0+VPH.*s.mol.(specimen{:}).HE;
    s.mol.(specimen{:}).FJ = Bv.*JJ-Dv.*JJ.^2+Hv.*JJ.^3;
    % Term values for the rotating vibrator (III-121), higher orders of vph
    % available in CARS.MOL, use all. Note: WE corresponds to omega_e*x_e
    % and so forth.
    s.mol.(specimen{:}).T = s.mol.(specimen{:}).FJ + ...
        s.mol.(specimen{:}).WE*VPH - ...
        s.mol.(specimen{:}).WX*VPH.^2+...
        s.mol.(specimen{:}).WY*VPH.^3+...
        s.mol.(specimen{:}).WZ*VPH.^4;
    
end
end

function s = calculateLinePositions(s)
% for every specimen, calculate the line position based on the CARS
% selection rules:
% DV = 1 (only ro-vibrational implemented)
% DJ = -2,0,2 (O,Q,S branch)
% Pure Rotational not yet implemented! This can be easily done if another
% case is introduced (let's say, 'ROT'). However, this needs some adaptions
% to the later calculations of the amplitudes of Chi(3). In particular, the
% Placzek-Teller-coefficients and Hermann-Wallis-Factors have to be
% included.

for specimen = s.species
    % Q branch: DJ = 0
    s.transitions.(specimen{:}).Q = diff(s.mol.(specimen{:}).T,1,2);
    % O branch: DJ = -2
    % initiallize with NaNs to have same size arrays with the index still
    % representing ground state J
    s.transitions.(specimen{:}).O = NaN(size(s.transitions.(specimen{:}).Q));
    % subtract a shifted array to account for DJ=-2.
    s.transitions.(specimen{:}).O(3:end,:) = s.mol.(specimen{:}).T(1:end-2,2:end)-s.mol.(specimen{:}).T(3:end,1:end-1);
    % S branch: DJ = +2
    % same as for O, shift in the other direction
    s.transitions.(specimen{:}).S = NaN(size(s.transitions.(specimen{:}).Q));
    s.transitions.(specimen{:}).S(1:end-2,:) = s.mol.(specimen{:}).T(3:end,2:end)-s.mol.(specimen{:}).T(1:end-2,1:end-1);
    
end
end

function s = calculateLineWidths(s)
% this function returns the linewidths in units cm-1/bar, i.e. normalized
% to the pressure to make it independent.
for specimen = s.species
    switch specimen{:}
        case 'N2'
            % [1] From Koszykowski, Rahn, Palmer, Coltrin: Theoretical and
            % Experimental Studies of High-Resolution Inverse Raman Spectra
            % of N2 at 1-10 atm, Phys Chem 91,1,1987
            % AND
            % [2] Rahn, Palmer, Koszykowski, Greenhalgh: Comparison of rotationally
            % inelastic collision models for q-branch raman spectra of N2,
            % Chemical Physical Letters 133,6,1987. Primary source probably
            % L. A. Rahn and R. E. Palmer, "Studies of nitrogen self-broadening at high temperature with inverse Raman spectroscopy," J. Opt. Soc. Am. B 3, 1164-1169 (1986)
            % AND
            % [3] Sitz, Greg O., and R. L. Farrow. "Pump probe measurements of state to state rotational energy transfer rates in N2 (v= 1)." The Journal of chemical physics 93.11 (1990): 7883-7893.
            
            % some handles for better readability
            % data from [2]
            alpha = 0.0231;
            tcorr = sqrt(295./s.T).*(1-exp(-0.1487))./(1-exp(-0.1487*s.T/295)); % temperature correction, [2] eq (1)
            beta = 1.67;
            delta = 1.21; % should be 1.26!
            Ei = s.mol.N2.FJ(:,1)*sconst('h')*sconst('c');
            Ej = s.mol.N2.FJ(:,1)'*sconst('h')*sconst('c');
            dEij = Ej-Ei; % this is \Delta E_{ij}
            kT = sconst('kB')*s.T; % this is kT
            
            % this is without the pressure, this will be multiplied later
            % on so it is not included in the linewidth calculation yet
            % the transponings are necessary to account for the switching
            % between ij and ji.
            gammaji = tril((alpha*tcorr*((1+1.5.*Ei/(kT*delta))./(1+1.5.*Ei/kT)).^2.*exp(-beta*dEij/kT))',-1);
            % this is microscopic reversibility [1], eq (4.2):
            gammaij = triu((((2*s.J+1)./(2*s.J+1)'.*gammaji.*exp(dEij/kT)'))',1);
            
            % combine the upper and lower triangular matrices to the entire
            % S matrix
            S = gammaji + gammaij;
            
            % account for selection rules DJ = +-2
            % this basically means, that Ji and Jj are either both even or
            % uneven. this introduces a checkerboard like pattern. this is
            % multiplied onto the calculated S matrix to make the
            % forbidden entries zero.
            checkerboard = ones(s.J(end)+1);
            checkerboard(2:2:(s.J(end)+1)^2)=0;
            S = S .* checkerboard;
            
            % this is equation [2] (4.3). take the upper triangular matrix of
            % gammaij (because j<i) without the diagonal, multiply by 2,
            % sum.
            % transpose it to have it as a column vector
            s.linewidths.N2 = sum(2*S,1)';
            
            
    end
end
end

function s = calculateBoltzmannFraction(s)
% calculate the differential boltzmann fraction used for the amplitudes of
% chi(3)
% max. vibrational state and vibrational partition function qv:
% a shortcut to hc/kT for easier readability.
hc_KT = sconst('hc_k')./s.T; % this is hc/kT

for specimen = s.species
    
    
    
    % populations
    % the following looks a little complicated due to the vectorization
    % it just is (2J+1)*GJ*exp(-FJ) for every temperature in the range
    % (accounting for the different degeneracies for even and odd Js)
    % to do this vectorized, FJ is reshaped from NJxNV to NJ*NVx1 and the
    % result is reshaped back. See equation (III,161)
    % This is the numerator in the Boltzmann equation for all rotational
    % energy levels
    RS = (2.*s.J+1)'.*s.mol.(specimen{:}).GJ'.*exp(reshape(reshape(-s.mol.(specimen{:}).FJ,[],1).*hc_KT,length(s.J),length(s.V),[]));
    % This is the population distribution, accounting for the population of
    % the vibrational states
    VS(1,:,:) = exp(-s.mol.(specimen{:}).G0V'.*hc_KT); % add a singleton dimension for later repmat
    % this is now the combination of vibrational and rotational states
    RSVS = repmat(VS,length(s.J),1,1).*RS;
    fB = RSVS./sum(sum(RSVS,1),2); % get the boltzmann distribution by normalizing to the sum of all states (i.e. the partition function)
    
    % now, for every transition, get the differential population
    % distribution
    % Note: diff uses does i-(i+1), we want the opposite, so invert it.
    delta.Q = -diff(fB,1,2);
    delta.O = NaN(size(delta.Q));
    delta.O(3:end,:,:) = -(fB(1:end-2,2:end,:)-fB(3:end,1:end-1,:));
    delta.S = NaN(size(delta.Q));
    delta.S(1:end-2,:,:) = -(fB(3:end,2:end,:)-fB(1:end-2,1:end-1,:));
    
    % Extract the highest J level with significant population difference for the later
    % estimation of the desired grid size. This is useful because higher
    % J's tend to have lower linewidth. But if they are not populated, then
    % it doesn't make sense to resolve them.
    % Threshold is 0.001*max, i.e. 0.1% contribution. This does not mean that
    % the lines with lower values are disregarded, only that their
    % resolution is lower. limit to max s.J+1 (due to matlab starting arrays at index 1).
    s.JMax=min([round(1.5*find((delta.Q(:,1)./max(delta.Q(:,1)))>0.001,1,'last')) max(s.J)+1]);
    
    % rotational partition function, used later for rotational diffusion
    QR = s.mol.(specimen{:}).GJ.*(2*s.J+1)*exp(-s.mol.(specimen{:}).FJ(:,1)*hc_KT); % rotational partition function for every temperature
    
    % store for later use
    s.delta.(specimen{:}).Q = delta.Q;
    s.delta.(specimen{:}).O = delta.O;
    s.delta.(specimen{:}).S = delta.S;
    s.fB_rot = RS/QR; % for rotational diffusion. this is a normalized population distribution without accounting for the population of the vibrational band
end

end

function s = chiamplitudes(s)
% this function computes the amplitudes of the resonant susceptibility
% note: at this point, this is without the number density N, which will be
% treated later
for specimen = s.species
    % constants
    zeta   = s.mol.(specimen{:}).AC1;
    vp1 = s.V(2:end);   % this is v+1, i.e. ground state +1
    xi = s.mol.(specimen{:}).AG;
    
    % Placzek-Teller coefficients
    bjjO = s.J.*(s.J-1)./((2*s.J-1).*(2*s.J+1));
    bjjQ = s.J.*(s.J+1)./((2*s.J-1).*(2*s.J+3));
    bjjS = (s.J+1).*(s.J+2)./((2*s.J+3).*(2*s.J+1));
    
    % Herrmann-Wallis factors
    eta = 2*s.mol.(specimen{:}).BE/s.mol.(specimen{:}).WE;    % eta in Marrocco, gamma in Klemperer
    % for O and S-branch
    eta_os = 4*s.mol.(specimen{:}).BE*s.mol.(specimen{:}).GAM/(s.mol.(specimen{:}).WE*s.mol.(specimen{:}).DGAMDR*s.mol.(specimen{:}).RE); % this is from the left central part in Buckingham, p47. This corresponds to 4*BE/WE*(alpha||-alpha_|_)_e/(alpha||-alpha_|_)'_e
    FS = (1-eta_os.*(2*s.J+3)).^2;    % buckingham, p47
    FO = (1+eta_os*(2*s.J-1)).^2;    % buckingham, p47
    
    % ADD A SELECTOR FOR THE HW MODEL HERE!
    switch s.HWFactors
        case 'JK'
            FQ = 1-3*eta^2*s.J.*(s.J+1)/2; % klemperer, below equation (14)
        case 'LBY'
            a1 = -2.7; % see Marrocco paper, Table1
            FQ = (1-3*eta^2*(a1+1)*s.J.*(s.J+1)/4).^2;
        case 'TB'
            a1 = -2.7;  % see Marrocco paper, Table1
            p2p1_iso = 0.31;    % see Marrocco paper, Table1
            p2p1_ani = 0.57;    % see Marrocco paper, Table1
            FQ_iso = 1-(3*(a1+1)/2-4*p2p1_iso)*eta^2*s.J.*(s.J+1);
            FQ_ani = 1-(3*(a1+1)/2-4*p2p1_ani)*eta^2*s.J.*(s.J+1);
        otherwise
            error('Hermann-Wallis-Factors should either be JK, LBY or TB')
    end
    
    % O-branch
    Opar = 2/15*bjjO.*FO*xi^2;
    Operp = 3/4*Opar;
    
    % Q-branch
    switch s.HWFactors
        case 'TB'
            Qpar = (1+4/45*xi^2*bjjQ).*FQ_iso;
            Qperp = 1/15*xi^2*bjjQ.*FQ_ani;
        otherwise
            Qpar = (1+4/45*xi^2*bjjQ).*FQ;
            Qperp = 1/15*xi^2*bjjQ.*FQ;
    end
    % S-branch
    Spar = 2/15*bjjS.*FS*xi^2;
    Sperp = 3/4*Spar;
    
    % compute the amplitudes
    s.chiamp.(specimen{:}).O = 1e18 * s.delta.(specimen{:}).O .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Opar'+sind(s.theta)*sind(s.phi)*Operp');
    s.chiamp.(specimen{:}).Q = 1e18 * s.delta.(specimen{:}).Q .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Qpar'+sind(s.theta)*sind(s.phi)*Qperp');
    s.chiamp.(specimen{:}).S = 1e18 * s.delta.(specimen{:}).S .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Spar'+sind(s.theta)*sind(s.phi)*Sperp');
end
end

function s = selectTransitions(s)
for specimen = s.species
    % transitions outside range have been blanked with NaN
    % serialize all transitions with ground state J and V. Also store a
    % flag for the transition (-2 for O, 0 for Q, +2 for S)
    J_extended = repmat(s.J',size(s.delta.(specimen{:}).Q,2),1);
    V_extended = reshape(repmat(s.V(1:end-1)',1,size(s.delta.(specimen{:}).Q,1))',[],1);
    O = -2*ones(size(J_extended));
    Q = 0*ones(size(J_extended));
    S = 2*ones(size(J_extended));
    transitions = [s.transitions.(specimen{:}).O(:),V_extended,J_extended,O;...
        s.transitions.(specimen{:}).Q(:),V_extended,J_extended,Q;...
        s.transitions.(specimen{:}).S(:),V_extended,J_extended,S;];
    
    % blank all transitions that are far outside the plotrange (and then
    % some...)
    transitions(transitions(:,1)<s.ROI(1)-s.ROIExp,:)=NaN;
    transitions(transitions(:,1)>s.ROI(2)+s.ROIExp,:)=NaN;
    
    % remove all lines containing nan
    transitions(any(isnan(transitions),2),:)=[];
    
    % store in struct
    s.transitions_inc.(specimen{:})=transitions;
    
    % extract the linewidths for these transitions
    s.linewidths_inc.(specimen{:}) = s.linewidths.(specimen{:})(transitions(:,3)+1);
    
    % extract the chiamps (based on V, J and branch)
    % this is easiest in a for loop
    
    % preallocate
    chiamp = zeros(size(s.linewidths.(specimen{:})));
    for ii =1:size(transitions,1)
        % which branch?
        if transitions(ii,4)==-2
            branch = 'O';
        elseif transitions(ii,4)==0
            branch = 'Q';
        elseif transitions(ii,4)==2
            branch = 'S';
        end
        % extract the chiamps
        chiamp(ii,:) = s.chiamp.N2.(branch)(transitions(ii,3)+1,transitions(ii,2)+1);
    end
    
    % store in output struct
    s.chiamp_inc.(specimen{:})=chiamp;
end
end

function s = createTheoSusc(s)
% This creates the spectrum of theoretical susceptibilities
% At the moment, this is only implemented for N2
% generate the wavenumberarray based on the lowest and highest transition
% value. grid size is 0.1 * minimum linewidth

min_wavenumber = s.ROI(1)-s.ROIExp;
max_wavenumber = s.ROI(2)+s.ROIExp;
% choose the resolution to have 10 points in the smallest linewidth for a
% significantly populated highest J. this has proven to be quite robust (in terms of comparison
% with CARSFT)
s.res_wavenumber = round(s.LineWidthMultiplier*0.1*min(s.linewidths.N2(1:s.JMax)*s.P),1,'significant');
s.wavenumberarray = min_wavenumber:s.res_wavenumber:max_wavenumber;
% Note: DO NOT explicitly include the line centers in the wavenumberarray because unevenly
% spaced grids generate trouble when convolving.

% total number density
s.N = sconst('NA')*273.15/s.T*(s.P/sconst('molarvolume'));

% get indices of transitions with significant contribution to the signal
% Kidx = find(s.chiamp_inc.N2>0.01*max(s.chiamp_inc.N2));
Kidx = find(s.chiamp_inc.N2>0.001*max(s.chiamp_inc.N2)); % use all transitions with more than 0.1% contribution

% array of final linewidth (including multiplier for foreign gas broadening
% and pressure), units: cm-1
gamma = s.linewidths_inc.N2*s.LineWidthMultiplier*s.P; % multiply with pressure

switch s.Model
    case 'I'
        % in this case, the spectrum consists of the lorentzians of all
        % individual transitions
        s.chi = complex(zeros(1,length(s.wavenumberarray),'single'));
        for K = Kidx'
            % solve each transition and sum it
            s.chi = s.chi + s.N*s.chiamp_inc.N2(K)./(s.transitions_inc.N2(K,1)-s.wavenumberarray-1i*gamma(K)/2);
        end
        
    case 'V'
        % this is basically the same as the isolated lines model, with the
        % difference that the individual lines are convolved with a
        % gaussian with Doppler width before adding them to the spectrum
        s.chi = complex(zeros(1,length(s.wavenumberarray),'single'));
        for K = Kidx'
            % generate a lorentzian for every transition
            currtrans = s.N*s.chiamp_inc.N2(K)./(s.transitions_inc.N2(K,1)-s.wavenumberarray-1i*gamma(K)/2);
            % calculate the doppler width for the current transition
            doppler_width = 4.3014e-7*s.transitions_inc.N2(K,1)*sqrt(s.T/s.mol.N2.CMASS);
            doppler_kernel = gaussian(-10*doppler_width:s.res_wavenumber:10*doppler_width,0,doppler_width);
            % add the convolved spectrum to the complex spectrum
            s.chi = s.chi + conv(currtrans,doppler_kernel,'same');
        end
        
    case 'R'
        % We need the denominator from equation (4) in Hall and Greenhalgh.
        % the numerator and the term before the fraction is just the
        % ordinary lorentzian equation for chi(3). basically, the
        % transitions are weighted with the factor 1/(1+i<[f/tau]/Dv>)
        f = s.fB_rot(:,1:end-1);    % this is f in equation (4)
        tau_j = 1./(s.LineWidthMultiplier*s.P*s.linewidths.N2); % this is tau_J in equation(4)
        omega_vvp1 = diff(s.mol.N2.T,1,2); % this is omega_v,v+1(omega_r) (i.e. all line positions of the q branches for each v)
        wavenumberarray(1,1,:)=s.wavenumberarray; % add two singleton dimensions
        Dv = wavenumberarray-repmat(omega_vvp1,1,1,length(s.wavenumberarray))-1i./tau_j;
        
        denominator = squeeze(1+1i*sum(f./tau_j./Dv,1)); % expectation value in J direction
        denominator = denominator./denominator.^2; % this is still somewhat unclear
        
        % preallocate the array
        s.chi = complex(zeros(1,length(s.wavenumberarray),'single'));
        % iterate through the transitions
        for K = Kidx'
            % rotational diffusion is only available for q-branch
            if s.transitions_inc.N2(K,4)==0
                % this is equation (4) in ref Hall and Geenhalgh
                v = s.transitions_inc.N2(K,2)+1; % index of the groundstate v
                
                % standard lorentzian for this transition
                currtrans = s.N*s.chiamp_inc.N2(K)./(s.transitions_inc.N2(K,1)-s.wavenumberarray-1i*gamma(K)/2);
                currtrans = currtrans.*denominator(v,:);
                s.chi = s.chi + currtrans;
            else
                % use isolated lines for O and S branch
                s.chi = s.chi + s.N*s.chiamp_inc.N2(K)./(s.transitions_inc.N2(K,1)-s.wavenumberarray-1i*gamma(K)/2);
            end
        end
    otherwise
        error('Model %s not implemented.',s.Model);
end

% non-resonant contribution to chi
chinr = computeNonResonantContribution(s);
% store the real part and the resonant squared part in the result for
% library generation
s.chi_real = real(s.chi);
s.chi_res = real(s.chi).^2 + imag(s.chi).^2;
s.chi = s.xN2*s.chi+chinr;
s.spectra.theosusc = abs(s.chi).^2; % this is the squared susceptibility!

end


function chinr = computeNonResonantContribution(s)
% Non-Resonant contribution is polarization dependend
polfactor = sind(s.theta)*sind(s.phi)/3+cosd(s.theta)*cosd(s.phi);

% so far, only N2 implemented. everything that is not N2 contributes
% with the same non resonant background.
% PANR is the mole fraction of the other species
chinr = polfactor*s.mol.N2.CHINR*s.N*s.xN2; % this is the non-resonant contribution of N2

% buffer gas non resonant susceptibility
buffergas = 1-s.xN2;
% again, this is given in 10^18 cm^3/(erg*amagat), convert to contribution
% per molecule
chinr_buffergas = s.chinr_buffergas/(sconst('NA')/sconst('molarvolume'));
chinr = chinr + polfactor*s.N*chinr_buffergas*buffergas;

end

function s = createCARSspectrum(s)

% convolve the square susceptibility with the convolved instrument
% functions

convrange = -10:s.res_wavenumber:10; % range for convolution kernels
kernel = gaussian(convrange,0,s.linewidth);

% convolve the square of the spectrum and take the sqrt again
s.spectra.CARS = sqrt(conv(s.spectra.theosusc,kernel,'same'));

end

function c = sconst(input)
%return spectroscopic constant

switch input
    case 'h'
        c=6.62607004e-34; % plancks constant in m^2 kg / s (i.e. Js)
    case 'c'
        c = 29979245800; % speed of light in cm/s
    case 'kB'
        c = 1.38064852e-23; % boltzmann in m^2 kg / (s^2 K)
    case 'NA'
        c = 6.02214086e23; % avogadro's constant (mol^-1)
    case 'R'
        c = 8.314;          % ideal gas constant, (J/(mol K))
    case 'molarvolume'
        c = sconst('R')*273.15/101325*1e6;      % ideal gas molar volume at 101325 Pa and 273.15 K in L/kmol
    case 'hc_k'
        c = sconst('c')*sconst('h')/sconst('kB');     % hc/k (c in cm/s)
    otherwise
        error('Constant %s not yet defined',input)
end
end

function plotSpec(s,figid)
% this is a function to conveniently plot a result struct of marsft
% if figid is provided, use it. otherwise create a new one.
if nargin == 2
    figure(figid);clf;
else
    figure;
end
hold all;
if isfield(s.spectra,'theosusc')
yyaxis right
plot(s.wavenumberarray,s.spectra.theosusc,'-.','Color',[0.5 0.5 0.5])
ylabel('Theoretical Susc.')
yyaxis left
ha=gca;
ha.YAxis(2).Color='k';
ha.YAxis(1).Color='k';
end
plot(s.wavenumberarray,s.spectra.CARS,'k-','LineWidth',2);

box on
grid on
xlabel('Raman shift in cm-1')
ylabel('CARS Signal')
title(sprintf('T: %d, P: %.2f, XN2: %.2f, Instr.: %.2f, Linewidth Mult.: %.2f, HW: %s',s.T,s.P,s.xN2,s.linewidth,s.LineWidthMultiplier,s.HWFactors));

end