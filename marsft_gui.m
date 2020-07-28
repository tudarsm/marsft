% marsft_gui
% features: simulate cars spectra in (almost) real time and compare
% manually to measured data. this is useful to estimate the linewidth
% multiplier range that should be included in the database

% make the handles struct global
clear global h
global h
% create the figure
h.figure = figure(1337);

% clear the figure
clf;
% set the size
h.figure.Units = 'normalized';
h.figure.Position(3) = 0.6;
h.figure.Position(4) = 0.6;
h.figure.Position(1) = 0.2;
h.figure.Position(2) = 1-1.5*h.figure.Position(4);
% create the axes for plotting and style it
h.ax = axes();
h.ax.Position(1)=0.1;
h.ax.Position(2)=0.1;
h.ax.Position(3)=0.5;
h.ax.OuterPosition(4)=1;

% add buttons to load experimental data
h.load = uicontrol(h.figure,'String','Load exp. data','Style','togglebutton','Value',0);   % Temperature
h.load.Units='normalized';
h.load.Position(1)=0.65;
h.load.Position(2)=h.ax.Position(2);
h.load.Position(3)=0.1;
h.load.Position(4)=0.075;
h.load.Callback=@loadData;
% add buttons to load experimental data
h.clear = uicontrol(h.figure,'String','Clear exp. data','Style','togglebutton','Value',0);   % Temperature
h.clear.Units='normalized';
h.clear.Position(1)=h.load.Position(1)+0.1;
h.clear.Position(2)=h.ax.Position(2);
h.clear.Position(3)=0.1;
h.clear.Position(4)=0.075;
h.clear.Callback=@clearData;


% create the controls for marsft
h.s.T = uicontrol(h.figure,'String','Temperature','Style','slider','Min',0,'Max',3000,'Value',2000,'SliderStep',[1/3000,100/3000]);   % Temperature
h.s.P = uicontrol(h.figure,'String','Pressure','Style','slider','Min',0,'Max',50,'Value',1,'SliderStep',[1/50,10/50]);   % Pressure
h.s.Inst = uicontrol(h.figure,'String','Instrumental','Style','slider','Min',0.5,'Max',10,'Value',1,'SliderStep',[1/100,1/10]);   % Linewidth
h.s.LinWidMult = uicontrol(h.figure,'String','Linewidth Multiplier','Style','slider','Min',0.1,'Max',10,'Value',1,'SliderStep',[1/1000,1/100]);   % linewidth multiplier
h.s.xN2 = uicontrol(h.figure,'String','N2 Mole fraction','Style','slider','Min',0,'Max',1,'Value',1,'SliderStep',[0.01,0.1]); % N2 mole fraction
% model not implemented yet
h.s.phi = uicontrol(h.figure,'String','Phi','Style','slider','Min',0,'Max',359,'Value',0,'SliderStep',[1/360,1/36]); % phi
h.s.theta = uicontrol(h.figure,'String','Theta','Style','slider','Min',0,'Max',359,'Value',0,'SliderStep',[1/360,1/36]); % theta
% roi not implemented yet
h.s.waveshift = uicontrol(h.figure,'String','Waven. Shift','Style','slider','Min',-10,'Max',10,'Value',0,'SliderStep',[1/1000,1/10]); % horizontal shift
h.s.waveexp = uicontrol(h.figure,'String','Waven. Expansion','Style','slider','Min',-5,'Max',5,'Value',0,'SliderStep',[1/1000,1/100]); % wavenumber expansion
% buffer gas susceptibility
h.s.buffer = uicontrol(h.figure,'String','Buffer gas susc.','Style','slider','Min',0,'Max',25,'Value',8.5,'SliderStep',[25/1000,25/100]); %
% model
h.s.model = uicontrol(h.figure,'String','Model','Style','slider','Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]); %

% set the positions of all controls
cheight = 0.05;
cwidth = 0.15;
cwidth2 = 0.075;
cdist = 0.01;
starty = h.ax.Position(2)+h.ax.Position(4)-cheight;
fn = fieldnames(h.s);
for ii = 1:length(fn)
    h.s.(fn{ii}).Units='normalized';
    h.s.(fn{ii}).Position(1)=0.65; % x coordinate
    h.s.(fn{ii}).Position(2)=starty-(ii-1)*(cheight+cdist); % y coordintae
    h.s.(fn{ii}).Position(3)=cwidth; % width
    h.s.(fn{ii}).Position(4)=cheight; % height
    % add edittexts field next to it with the same width and height which
    % carries the value of the slider
    h.sv.(fn{ii}) = uicontrol(h.figure,'Style','edit','Units','normalized');
    h.sv.(fn{ii}).Position(1)=0.66+cwidth;
    h.sv.(fn{ii}).Position(2)=h.s.(fn{ii}).Position(2);
    h.sv.(fn{ii}).Position(3)=cwidth2; % width
    h.sv.(fn{ii}).Position(4)=cheight; % height
    % set the string
    h.sv.(fn{ii}).String = sprintf('%d',h.s.(fn{ii}).Value);
    % initialize the model string to a valid value
    if strcmp(h.s.(fn{ii}).String,'Model')
        switch h.s.(fn{ii}).Value
            case 0
                modelstr = 'I';
            case 1
                modelstr = 'V';
            case 2
                modelstr = 'R';
        end
        h.sv.(fn{ii}).String = modelstr;
    end
    % add labels next to it
    h.st.(fn{ii}) = uicontrol(h.figure,'Style','text','Units','normalized');
    h.st.(fn{ii}).Position(1)=h.sv.(fn{ii}).Position(1)+cwidth2+0.01;
    h.st.(fn{ii}).Position(2)=h.s.(fn{ii}).Position(2);
    h.st.(fn{ii}).Position(3)=cwidth2; % width
    h.st.(fn{ii}).Position(4)=cheight; % height
    % set the string
    h.st.(fn{ii}).String = sprintf('%s',h.s.(fn{ii}).String);
    
    % set the callbacks
    h.s.(fn{ii}).Callback = @cCallbackSlider;
    h.sv.(fn{ii}).Callback = @cCallbackEdit;
end
plotCars();

% callback for sliders
function cCallbackSlider(hObject,callbackdata)
global h
% if a slider was moved, update the edittext
% get a handle on the edittext
fn = fieldnames(h.s); % all fieldnames
id = find(structfun(@(x) x==hObject,h.s)); % get the id for the fieldnames array
% update the value
h.sv.(fn{id}).String = sprintf('%d',h.s.(fn{id}).Value);
% if it is the model slider, update the model identifier
if strcmp(h.s.(fn{id}).String,'Model')
    h.s.(fn{id}).Value=round(h.s.(fn{id}).Value);
    switch h.s.(fn{id}).Value
        case 0
            modelstr = 'I';
        case 1
            modelstr = 'V';
        case 2
            modelstr = 'R';
    end
    h.sv.(fn{id}).String = modelstr;
end
plotCars();
end

% callback for sliders
function cCallbackEdit(hObject,callbackdata)
global h
% if a slider was moved, update the edittext
% get a handle on the edittext
fn = fieldnames(h.s); % all fieldnames
id = find(structfun(@(x) x==hObject,h.sv)); % get the id for the fieldnames array
% update the value
h.s.(fn{id}).Value = str2double(h.sv.(fn{id}).String);
plotCars();
end

function plotCars()
% get the values from the slider struct
global h
T = h.s.T.Value;
P = h.s.P.Value;
XN2 = h.s.xN2.Value;
LineWidthMult = h.s.LinWidMult.Value;
CHINR = h.s.buffer.Value;
linwid = h.s.Inst.Value;
WMIN = 2270;
WMAX = 2340;
PHI = h.s.phi.Value;
THETA = h.s.theta.Value;
Model = h.sv.model.String;
wavenshift = h.s.waveshift.Value;
wavenexp = h.s.waveexp.Value;

s=marsft_sim('chinr',CHINR,'ROI',[WMIN WMAX],'T',T,'P',P,'Model',Model,'xN2',XN2,'LineWidthMultiplier',LineWidthMult,'linewidth',linwid,'phi',PHI,'theta',THETA);
% check if experimental data is available
if ~isfield(h,'dat')
    hold off
    h.p = plot(h.ax,s.wavenumberarray,s.spectra.CARS./sum(s.spectra.CARS),'k-');
else
    % plot the experimental spectrum
    hold off;
    expdat = h.dat(2,:)./sum(h.dat(2,:));
    h.pexp = plot(h.ax,h.dat(1,:),expdat,'r-');
    hold on;
    % interpolate the simulation on the wavenumber array
    destarray = s.wavenumberarray+wavenshift;
    destarray = linspace(destarray(1)-wavenexp,destarray(end)+wavenexp,length(destarray));
    marsft_interp = interp1(destarray,s.spectra.CARS,h.dat(1,:));
    marsft_interp = marsft_interp./sum(marsft_interp);
    h.p = plot(h.ax,h.dat(1,:),marsft_interp,'k-');
    h.pres = plot(h.ax,h.dat(1,:),expdat-marsft_interp,'Color',[0 0.5 1]);
end
%style the axes
h.ax.Box='on';
h.ax.XGrid='on';
h.ax.YGrid='on';
h.ax.Title.String='Simulated CARS Spectrum';
h.ax.XLabel.String='Wavenumber in cm^{-1}';
h.ax.XLabel.String='Intensity in a.u.';
% h.ax.XLim=[WMIN WMAX];
% h.ax.YLim(1)=0;
% h.ax.YLim(2)=max(s.spectra.CARS./sum(s.spectra.CARS));

end

function loadData(hObject,callbackdata)
global h
dat = load(uigetfile());
h.dat = dat.dat;
hObject.Value=0;
plotCars();
end

function clearData(hObject,callbackdata)
global h
if isfield(h,'dat')
    h=rmfield(h,'dat');
    plotCars();
end
hObject.Value=0;
end