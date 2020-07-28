function [f,fh] = marsft(lib,wavenumberarray,data,options)
%MARSFT This function tries to fit a cars spectrum to the given
%experimental data
% start a tic for the total runtime
fh.plotSpec = @plotSpec;
global singleshottime
singleshottime=tic;
global strcr
strcr='';

% if not supplied, use default options
if nargin == 3
    options = marsft_options(lib);
end

% convert temperature and linewidth multiplier to indexes for fitting
Tstart = find(lib.Ts==options.fitpar(1,1));
Tmin = find(lib.Ts==options.fitpar(1,2));
Tmax = find(lib.Ts==options.fitpar(1,3));
options.fitpar(1,:)=[Tstart, Tmin, Tmax];


% prepare the experimental data
% CARSFIT only normalizes data to peak == 1. After the fitting process, it
% is renormalized with the peak of the theoretical susceptibility. For the
% fitting quality, this does not really make a difference. The problem with
% normalizing to the peak intensity is that this is very much prone to
% noise. Normalizing to the sum should cancel out most of the noise and
% allows for very tight boundaries in the vertical expansion fitparameter.
normfact = sum(data,2);
data = data./normfact;

% normalize the variance with the data
% input units are sqrt(counts), just as data
% is custom variance data supplied?
if options.VarianceData ~= -1
    % yes, use it and scale with data:
    options.variance = repmat(options.VarianceData,length(normfact),1)./normfact;
    options.VarianceModel = 'custom';
else
    % no, use a standard variant:
    switch options.VarianceModel
        case 'none'
            % use no estimate of the variance...
            options.variance = 1;
        case 'std'
            % for a dataset with constant conditions, such as in a laminar
            % flame, the temporal variance is a good measure of
            % the noise
            options.variance = (std(data).^2.*ones(size(data)));
        case 'shot-noise'
            % if shot noise is the dominant noise source, the square root of the data itself
            % is a good estimate of variance
            shotnoise = data;
%             shotnoise(shotnoise>0) = sqrt(shotnoise(shotnoise>0));
%             shotnoise(shotnoise<0) = -sqrt(-shotnoise(shotnoise<0));
            options.variance = sqrt(abs(shotnoise)); % clearly, shot-noise cannot be negative. in case of detector noise, take the absolute value in order not to bias the output
    end
end
% in any case, set the global mean to 1 so the magnitude of the residual value is unaffected
% options.variance = options.variance./mean(options.variance(:));

% fit the spectrum. depending on the method, there are different options
% available
switch options.Method
    case 'mean'
        options.variance = 1;
        if strcmp(options.Display,'summary')
            fprintf('MARSFT: Fitting mean spectrum...')
            tmean=tic;
        end
        f = fitSpectrum(lib,wavenumberarray,mean(data,1),options);
        if strcmp(options.Display,'summary')
            fprintf('done in %.2f s. T = %d K.\n',toc(tmean),f.T)
        end
        
        if strcmp(options.Plot,'on')
            figure(99);clf;plotSpec(f);
        end
    case 'single-detailed'
        if strcmp(options.VarianceModel,'none')
            % if model 'none' is selected, create a vector with ones for the
            % variances of the single shots
            singleshotvariance=ones(size(data,1),1);
        else
            singleshotvariance = options.variance;
        end
        
        for ii = 1:size(data,1)
            options.variance = singleshotvariance(ii,:);
            if ii == 1
                f(1) = fitSpectrum(lib,wavenumberarray,data(ii,:),options);
                f(size(data,1))=f(1);
                % empty the last entry....this is a stupid workaround because it
                % seems to be impossible to preallocate a struct array with empty
                % values and just fieldnames...
                f(end)=f(end-1);
            else
                f(ii) = fitSpectrum(lib,wavenumberarray,data(ii,:),options);
            end
            if strcmp(options.Display,'summary')
                printProgress(f,ii)
            end
            
            if strcmp(options.Plot,'on')
                subplot(1,2,1);cla;plotSpec(f(ii));
                subplot(1,2,2);plot([f.T]);title(sprintf('Mean: %.2f K, RMS: %.2f K',mean([f.T]),std([f.T])));box on;grid on;
                drawnow
            end
        end
    case 'single-quick'
        % fit the mean first
        meanvariance = 1;
        if strcmp(options.VarianceModel,'none')
            % if model 'none' is selected, create a vector with ones for the
            % variances of the single shots
            singleshotvariance=ones(size(data,1),1);
        else
            singleshotvariance = options.variance;
        end
        % use the mean variance for weighting!
        options.variance=meanvariance;
        if strcmp(options.Display,'summary')
            fprintf('MARSFT: Fitting mean spectrum...')
            tmean=tic;
        end
        f_mean = fitSpectrum(lib,wavenumberarray,mean(data,1),options);
        if strcmp(options.Display,'summary')
            fprintf('done in %.2f s. T = %d K.\n',toc(tmean),f_mean.T)
        end
        % and now, fix all parameters except temperature
        options.fitpar(:,1)=f_mean.fitresult;
        % fix everything but temperature
        options.fitpar(2:end,2)=options.fitpar(2:end,1);
        options.fitpar(2:end,3)=options.fitpar(2:end,1);
        options.PopulationSize = options.NormPopulationSize*1;
        options.MaxStallGenerations = options.NormMaxStallGenerations*1;
        % preallocate the output array of structs
        f(size(data,1))=f_mean;
        % empty the last entry....this is a stupid workaround because it
        % seems to be impossible to preallocate a struct array with empty
        % values and just fieldnames...
        f(end)=f(end-1);
        if strcmp(options.Plot,'on')
            figure(99);clf;
            subplot(2,2,1);plotSpec(f_mean);
        end
        % note: the mean fit is disregarded!
        for ii = 1:size(data,1)
            if singleshotvariance~=1
                options.variance = singleshotvariance(ii,:);
            end
            f(ii) = fitSpectrum(lib,wavenumberarray,data(ii,:),options);
            if strcmp(options.Display,'summary')
                printProgress(f,ii)
            end
            if strcmp(options.Plot,'on')
                subplot(2,2,2);cla;plotSpec(f(ii));
                subplot(2,2,3);plot([f.T]);title(sprintf('Mean: %.2f K, RMS: %.2f K',mean([f.T]),std([f.T])));box on;grid on;
                subplot(2,2,4);plot([f.err]);title(sprintf('Goodness of fit: %g',f(ii).err));box on;grid on;
                drawnow
            end
        end
    otherwise
        error('Unknown method...');
end
end

function f = fitSpectrum(lib,wavenumberarray,data,options)
runtime=tic;
A=zeros(length(options.fitpar));
b=ones(length(options.fitpar),1);
% turn off ga display for type summary
if strcmp(options.Display,'summary')
    options.Display='off';
end
ga_options = optimoptions('ga','Display',options.Display,'FunctionTolerance',options.FunctionTolerance,'PopulationSize',options.PopulationSize,'UseVectorized',options.UseVectorized,'MaxStallGenerations',options.MaxStallGenerations,'PlotFcn',[]);
lb = options.fitpar(:,2);
ub = options.fitpar(:,3);
[x,~,~,~,population,scores] = ga(@(c) calcError(c,lib,wavenumberarray,data,options.variance,options.BufferGasSusc),size(options.fitpar,1),A,b,[],[],lb,ub,[],[1,4],ga_options);
[err,bestfit]=calcError(x,lib,wavenumberarray,data,options.variance,options.BufferGasSusc);
x(1)=lib.Ts(x(1));

% store in struct
f.fitresult = x;
f.T = x(1);
f.bestfit = bestfit;
f.data = data;
f.res = bestfit-data;
f.wavenumberarray = wavenumberarray;
f.options = options; %store the options struct for later reproducability
f.variance = options.variance;
f.err = err;
f.r2 = 1-err/sum((data-mean(data)).^2);
f.runtime=toc(runtime);
end

function [err,sim,data,res] = calcError(coeffs,library,wavenumber,data,variance,buffergassusc)

% get the set of chireal and chires based on the temperature and linewidth
% pairs. linear indexing is not faster as for loop, so stick with that for
% easier readability.

chireal = zeros(size(coeffs,1),size(library.chireal,1),'single');
chires2 = zeros(size(coeffs,1),size(library.chireal,1),'single');
for ii = 1:size(coeffs,1)
    chireal(ii,:) = library.chireal(:,coeffs(ii,1),coeffs(ii,4));
    chires2(ii,:) = library.chires2(:,coeffs(ii,1),coeffs(ii,4));
end

% LIBRARY
totnum = 6.022e23*273.15./library.Ts(coeffs(:,1))'*(library.P/22414);

% so far, only N2 implemented. everything that is not N2 contributes
% with the same non resonant background.
% PANR is the mole fraction of the other species
chinr = library.tnres*library.mol.N2.CHINR*totnum.*coeffs(:,3);

buffergas = 1-coeffs(:,3);
chinr_buffergas = buffergassusc/(6.022e23/22414);
chinr = chinr + library.tnres*totnum*chinr_buffergas.*buffergas;

% convolve the resulting spectrum with the additional desired gaussian
addgausswidth = sqrt(coeffs(:,2).^2-library.preconvolution^2);

% create the kernel
convrange = -10:library.res:10;
kernel = gaussian(convrange,0,addgausswidth);
wavenumberarray = library.wavenumberarray;
spectra = chires2.*coeffs(:,3).^2 + 2*chinr.*chireal.*coeffs(:,3) + chinr.^2;
sim=getConvAndInterp(coeffs,wavenumber,wavenumberarray,spectra,kernel);

% normalize simulation to the area
sim = sim./sum(sim,2).*coeffs(:,5);

% variance weighted least squares error

if variance == 1
    res = (data-sim);
    err = sum(res.^2./variance,2)./length(data);
else
    res = (data-sim);
    % this actually is the reduced chi square stats that works if variance is in correct units  
%     err = (1-sum((res.^2./variance),2).^2./(length(data)-1).^2).^2;
    % this is more robust
    err = sum((res.^2./variance),2).^2./length(data);
end


end

function sim = getConvAndInterp(coeffs,wavenumber,wavenumberarray,spectra,kernel)

sim = zeros(size(coeffs,1),length(wavenumber));
% converting to 16 bit integers for speed, saves about 30% time for the
% convolution
spectra=uint16(spectra./max(spectra,[],2)*65536);
kernel=uint16(kernel./max(kernel,[],2)*65536);
for ii = 1:size(coeffs,1)
    % take the sqrt of the convolved spectra
    spectrum = sqrt(conv(spectra(ii,:),kernel(ii,:),'same'));
    %     interpolate the result on the experimental grid
    destarray = linspace(wavenumberarray(1)-coeffs(ii,7)+coeffs(ii,6),wavenumberarray(end)+coeffs(ii,6),length(wavenumberarray));
    sim(ii,:) = interp1(destarray,spectrum,wavenumber);
end

end

function plotSpec(f)
cla;hold off;
yyaxis left
plot(f.wavenumberarray,f.data,'k-')
hold on;
plot(f.wavenumberarray,f.bestfit,'r-')
plot(f.wavenumberarray,f.res)
yyaxis right
plot(f.wavenumberarray,f.variance,'-','Color',[0.75 0.75 0.75])
box on
grid on
title(sprintf('T: %.2f, Linewidth: %.2f, Mult: %.2f, N2: %.2f, R^2: %g', f.fitresult(1),f.fitresult(2),f.fitresult(4),f.fitresult(3),f.r2));
end

function printProgress(f,ii)
global strcr
global singleshottime
strout = sprintf('Finished set %d of %d, T = %d, T Mean = %.2f K, T RMS = %.2f K...fitted in %.2f s.',ii,length(f),f(ii).T,mean([f(1:ii).T]),std([f(1:ii).T]),toc(singleshottime));
singleshottime=tic;
fprintf([strcr strout]);
strcr = repmat('\b',1,length(strout));
if ii==length(f)
    fprintf('\nFINISHED\n');
end
end

function state = gaplot1drange(options,state,flag)
% gaplot1drange Plots the mean and the range of the population.
% normalize the fitparameters to the lower and upper bounds to make them
% fit in one plot
generation = state.Generation;
score = (state.Population-options.InitialPopulationRange(1,:))./(options.InitialPopulationRange(2,:)-options.InitialPopulationRange(1,:));
smean = mean(score);
Y = smean;
L = smean - min(score);
U = max(score) - smean;
numpar = length(smean);
switch flag
    
    case 'init'
        hold all;
        for ii=1:numpar
            %             subplot(numpar,1,ii);
            if ismember(ii,[6,7])
                plotRange(ii) = errorbar(generation,Y(ii),L(ii),U(ii),'LineWidth',2);
            else
                plotRange(ii) = plot(generation,Y(ii),'LineWidth',1);
            end
            
            set(plotRange(ii),'Tag',sprintf('gaplot1drange%02d',ii));
        end
        title('Range of normalized Populations, Mean','interp','none')
        xlabel('Generation','interp','none')
        box on
        grid on
        legend('Temperature','LineWidth','XN2','LinMult','Int.Exp','Wave.Shift','Wave.Exp','location','northoutside','orientation','horizontal')
    case 'iter'
        for ii=1:numpar
            hold all;
            %             subplot(numpar,1,ii);
            plotRange = findobj(get(gca,'Children'),'Tag',sprintf('gaplot1drange%02d',ii));
            
            try
                newX = [get(plotRange,'Xdata') generation];
                newY = [get(plotRange,'Ydata') Y(ii)];
                newL = [get(plotRange,'Ldata') L(ii)];
                newU = [get(plotRange,'Udata') U(ii)];
                set(plotRange,'Xdata',newX,'Ydata',newY,'Ldata',newL,'Udata',newU);
            catch
                newX = [get(plotRange,'Xdata') generation];
                newY = [get(plotRange,'Ydata') Y(ii)];
                set(plotRange,'Xdata',newX,'Ydata',newY);
            end
        end
        xlim([1 generation+1]);
end
end

function state = gaplotbi(options,state,flag)
% gaplot1drange Plots the mean and the range of the population.
% normalize the fitparameters to the lower and upper bounds to make them
% fit in one plot
generation = state.Generation;
score = (state.Population(1,:)-options.InitialPopulationRange(1,:))./(options.InitialPopulationRange(2,:)-options.InitialPopulationRange(1,:));
Y = score;
numpar=length(score);
switch flag
    
    case 'init'
        hold all;
        for ii=1:numpar
            %             subplot(numpar,1,ii);
            plotRange(ii) = plot(generation,Y(ii),'LineWidth',1);
            
            set(plotRange(ii),'Tag',sprintf('gaplotbi%02d',ii));
        end
        title('Range of normalized Populations, Mean','interp','none')
        xlabel('Generation','interp','none')
        box on
        grid on
        legend('Temperature','LineWidth','XN2','LinMult','Int.Exp','Wave.Shift','Wave.Exp','location','northoutside','orientation','horizontal')
    case 'iter'
        for ii=1:numpar
            hold all;
            %             subplot(numpar,1,ii);
            plotRange = findobj(get(gca,'Children'),'Tag',sprintf('gaplotbi%02d',ii));
            
            newX = [get(plotRange,'Xdata') generation];
            newY = [get(plotRange,'Ydata') Y(ii)];
            set(plotRange,'Xdata',newX,'Ydata',newY);
            
        end
        xlim([1 generation+1]);
end
end
