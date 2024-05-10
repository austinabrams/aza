%% Automated Zeta Potential Analysis .m file
% Keithley Controller and Zeta Tracker
% Last updated October 3, 2022

clear all; close all;

yourmail = {''};  % Insert Email address(es) to send final FIG & JPEG
folder = ''; % Insert Local address to save FIG and DATA
chan = ''; %  % Insert Name of Experiment to go on Plot Title
d = 39.5e-6;    % 100e-9     % Insert Diameter/height channel
L = 7.83e-3;        % Insert Length (m) -> use calipers in lab
type = 'Cap';       % Enter Cap (capillary) or Sol (Solomon)
if strcmp(type,'Sol')
    A = 500e-9*5e-6;         % X-sect Area
    L = 5e-3;
else
    A = pi*(d/2)^2;     % X-sect Area
end
solv = 1;        % Type 1 (H2O) or 2 (MeOH) for solvent
sig = 1.75;        % Conductivity S/m 
pH = '7.2';             % pH measured w/ freshly calibrated meter
buf_mM = '10';
buf = '';       % Buffer name Tris / CAPS / H2SO4
elec_mM = '0';
elec = '';      % Electrolyte name NaCl / KCl
dateFilt = '';          % Date solution was filtered
L_ratio = round(L/8e-3);

% Range: L) 0-1000 V, I=0-1 mA; R) 0-1100 V, I=0-105 uA.
GPIB = 24;          % Open NI-VISA and Lookup primary address (after GPIB0:)
boardInd = 0;       % Lookup board index w/ NI-VISA (# before colon, 1 or 0) 
process = 2;        % Process (1: constant voltage polarity; 2: alternating voltage polarity with or without autoSwitch; 3: combination of 1 and 2, 4: maintain net EOF+EP transport distance per cycle)
tfinal = 5*24*3600;   % Total seconds of Zeta monitoring
% I_range = 200e-9;   % Enter predefined range, or comment out to calculate based on conductivity and channel area

% Enter Process 2 settings below:
I_autorange = 0;    % Enter 1 to use autorange, 0 to use I_range
splitter = 0;       % E1 to use splitter card
channels = '1';  % enter channels 'a,b,c' to alternate measurement
bufStor = 0;    % 0: disables storage, 1:enables
d1_smooth = 1*L_ratio;  % Control data smoothing (n>1: range = d1_smooth*3 pts, 0 or 1: range=3 pts)
speed = 1;      % Integration time in power line cycles (1: 1PLC (~10/s), 2: 0.1PLC (~20/s), 3: 0.01PLC
d1_thres = 0.33;                   % If D1now = d1_thres * aveD1_n1n2, endpt reached
minCurrChange = 0.04;
% P1, P2 for autoswitch endpoint detection
% 0.01PLC and 0.1PLC both produce 20 pts/s

% Process 1: Run constant voltage
if process==1    
	vout1 = +10*1e3*L;  %+30*1e3*L
    E_ratio = vout1/L/30e3;
    tCycle1 = tfinal;
    autoSwitch = 0;
    minSlope = NaN;    t_PlatHold = NaN;
    I_autorange = 1;    % Enter 1 to use autorange, 0 to use I_range
    I_range = 0;
    liveplot = 1;

% Process 2: Flip voltage polarity every tcycle
elseif process==2
	vout1 = +30*1e3*L;  %+30*1e3*L
    E_ratio = vout1/L/30e3;
    I_autorange = 1;    % Enter 1 to use autorange, 0 to use I_range
    autoSwitch = 1;
    tCycle1 = 90; % L_ratio
    tCycle = 180*L_ratio; 	% Max duration of cycle (auto switch if plateau found)
    t_PlatHold = 5*L_ratio; % time after endpoint is reached to hold V
    tCycle2 = tCycle;    
    vout2 = -vout1;
%     vout2 = -60*1e3*L;
    add0Vcyc = 0;       % Apply 0V for 5s before Vswitch
    E = vout1/L;             % E-field
    Ipred = sig*E*A;        % Predicted current reading (A)
    I_range = Ipred*1.1;   % Enter predefined range
    minSlope = 6e-7*Ipred; % Don't find endpoint / autoSwitch if slope now < minSlope
    d1_smooth = 1;  % Control data smoothing (n>1: range = d1_smooth*3 pts, 0 or 1: range=3 pts)
    liveplot = 1;

% Process 3: Add Volt-hold intervals when t>=hold(i) && found==1
elseif process==3
	vout1 = +62.5*1e3*L;  %+30*1e3*L
    E_ratio = vout1/L/30e3;
    I_autorange = 1;
    autoSwitch = 1;
    tCycle1 = 15*L_ratio;
    tCycle = 60*L_ratio/E_ratio; 	% Max duration of cycle (auto switch if plateau found)
    t_PlatHold = 5*L_ratio/E_ratio;
    tCycle2 = tCycle;
    vout2 = -vout1;
    vout_hold = vout1; % normal, set to vout1; 0V- set to 0
    E = vout1/L;             % E-field
    Ipred = sig*E*A;        % Predicted current reading (A)
    I_range = Ipred*1.1;   % Enter predefined range
    minSlope = 100; % Don't find endpoint / autoSwitch if slope now < minSlope
    d1_smooth = 1;  % Control data smoothing (n>1: range = d1_smooth*3 pts, 0 or 1: range=3 pts)
    liveplot = 0;
    add0Vcyc = 0;       % Apply 0V for 5s before Vswitch

    t_hold = [3 2 7200; 3 2 120; 8 2 24*60; 58 2 5*24*60]; % mins bias, mins aza, until cutoff mins
    t_i = t_hold(1,2); % start with 2 min aza, then 3 min bias
    for i = 1:3000
        if t_i < t_hold(1,3)
            n=1;
            hold_V(i) = t_i;
            hold_t(i) = t_hold(n,1);
            t_i = t_i + t_hold(n,1) + t_hold(n,2);
        elseif t_i < t_hold(2,3)
            n=2;
            hold_V(i) = t_i;
            hold_t(i) = t_hold(n,1);
            t_i = t_i + t_hold(n,1) + t_hold(n,2);
        elseif t_i < t_hold(3,3)
            n=3;
            hold_V(i) = t_i;
            hold_t(i) = t_hold(n,1);
            t_i = t_i + t_hold(n,1) + t_hold(n,2);
        else
            break;
        end
    end
    hold_V = hold_V*60;         % convert hold start mins to sec's
    hold_t = hold_t*60;         % convert hold mins to sec's
    a = 1;                      % index of hold_V cycles

% Process 4: maintain net EOF+EP transport distance per cycle
elseif process == 4    
  	vout1 = +30*1e3*L;  %+30*1e3*L
    E_ratio = vout1/L/30e3;
    I_autorange = 1;    % Enter 1 to use autorange, 0 to use I_range
    autoSwitch = 1;
    net_CV_cycle = 3; %[zeta(EOF)-zeta(EP)]*E*eps/mu*(t_traverse/L0)
    tCycle1 = 15*net_CV_cycle; % L_ratio
    zeta_EP = +15; % mV, sign matters
    tCycle = 180; 	% Max duration of cycle (auto switch if plateau found)
%     t_PlatHold = 5*L_ratio; % time after endpoint is reached to hold V
    tCycle2 = tCycle;    
    vout2 = -vout1;
%     vout2 = -60*1e3*L;
    add0Vcyc = 0;       % Apply 0V for 5s before Vswitch
    E = vout1/L;             % E-field
    Ipred = sig*E*A;        % Predicted current reading (A)
    I_range = Ipred*1.1;   % Enter predefined range
    minSlope = 6e-7*Ipred; % Don't find endpoint / autoSwitch if slope now < minSlope
    d1_smooth = 1;  % Control data smoothing (n>1: range = d1_smooth*3 pts, 0 or 1: range=3 pts)
    liveplot = 0;
    
end

% DO NOT CHANGE PARAMETERS BEYOND THIS POINT

% Specify Keithley (1=6517B, 2=6517A, 3=2410)
% Voltage ranges. 1:0-1000. 2:0-100.3:0-1100
vout = vout1;
vout_last = vout1;
if autoSwitch == 0
    t_PlatHold = NaN;
end

E = vout1/L;             % E-field
%E2 = vout2/L;
Ipred = sig*E*A;        % Predicted current reading (A)

if solv == 1
    eps = 78.4*8.85e-12;	% Permittivity (water) * vaccuum (MeOH = 32.7, h2o = 78.4)
    mu = 8.9e-4;          % Viscosity (water) (MeOH = 5.85e-4, h2o = 8.9e-4)
    t1_skip = 2*L_ratio/E_ratio;      	% Enter sec's -> start reference slope
    t2_skip = 4*L_ratio/E_ratio;     	% Enter sec's -> start comparison slope (endpoint? == 33% run slope)

elseif solv == 2
    eps = 32.7*8.85e-12;	% Permittivity (water) * vaccuum (MeOH = 32.7, h2o = 78.4)
    mu = 5.85e-4;          % Viscosity (water) (MeOH = 5.85e-4, h2o = 8.9e-4)
    t1_skip = 3;      	% Enter sec's -> start reference slope
    t2_skip = 5;     	% Enter sec's -> start comparison slope (endpoint? == 33% run slope)

end

name = strcat(chan,',',32,...
...%  	num2str(d*1e6),'µm,',32,...
...% 	num2str(E/1e3),'kV/m,'...
num2str(E/1e3),'kV');
...% ,32,num2str(tCycle),'s,');
...% num2str(E2/1e3),'kV @',32,num2str(tCycle),'s'...
...%     res,32,dateFilt,32,...
...% 	buf_mM,'mM',32,buf,32,elec_mM,'mM',32,elec,',',32,...
...% 	'pH',32,pH,',',32,...
...%     '\sigmaEA=',num2str(round(1e9*Ipred)),32,'nA')


% Open GPIB Connections to Keithley
obj6517 = instrfind('Type', 'gpib', 'BoardIndex', boardInd, 'PrimaryAddress', GPIB, 'Tag', '');
if isempty(obj6517)
    obj6517 = gpib('NI', boardInd, GPIB);
else
    fclose(obj6517);
    obj6517 = obj6517(1);
end
fopen(obj6517);

% Reset System
% fprintf(obj6517, '*RST');             	% resets to SCPI parameters
fprintf(obj6517,':SYST:RNUM:RES');
fprintf(obj6517,':SYST:ZCOR ON');
fprintf(obj6517,':SYST:ZCH OFF');

% Setting data storage into the buffer
fprintf(obj6517,':DISP:ENAB ON');
fprintf(obj6517,':SENS:FUNC "CURR"');          % Measure 'VOLTage[:DC]','CURRent[:DC]','RESistance','CHARge'

fprintf(obj6517,':SYST:TSC ON');             	% Enable external temperature readings
fprintf(obj6517,':UNIT:TEMP C');             	% Select temperature units

fprintf(obj6517,':SYST:TST:TYPE REL');        	% RELative or RTClock for timestamp
fprintf(obj6517,':SYST:TST:REL:RES');        	% RESet realtive-time stamp to 0 sec. Max = 99,999.999999 sec
fprintf(obj6517,':TRAC:TST:FORM ABS');     	% Timestamp format ABSolute or DELTa
fprintf(obj6517,':TRAC:ELEM TST,RNUM,ETEM,VSO'); 	% Read READing,TSTamp,RNUMber,CHANnel,ETEMperature,VSOurce

fprintf(obj6517,':FORM:DATA ASC');
fprintf(obj6517,':FORM:ELEM READ,TST,RNUM,ETEM,VSO'); 	% Read READing,TSTamp,RNUMber,CHANnel,ETEMperature,VSOurce

if bufStor
    fprintf(obj6517,':TRAC:FEED:CONT ALW');        % Specify buffer control (NEVer, NEXT, ALWays, PRETrigger)
    fprintf(obj6517,':TRAC:POIN MAX');         	% Buffer size (1 to 50,000 (6517B only), MAX (max of 6517A), DEF=100)
    fprintf(obj6517,':TRAC:CLE');                  % Clears the buffer
else
    fprintf(obj6517,':TRAC:FEED:CONT NEV');
end

% fprintf(obj6517,':TRIG:COUN INF');          	% Set measure count to INFinite
% fprintf(obj6517,':TRIG:SOUR IMM');            	% Select control source (HOLD, IMMediate, TIMer, MANual, BUS, TLINk, EXTernal)

if I_autorange == 1
    fprintf(obj6517,':SENS:CURR:RANG:AUTO ON' ); 	% Set range based on present input signal
else
    fprintf(obj6517,[':SENS:CURR:RANG ', num2str(I_range)]); 	% Set range based on present input signal
end
fprintf(obj6517,[':SOUR:VOLT:RANG ', num2str(vout)]);
%% splitter function
if splitter
    fprintf(obj6517,[':ROUT:SCAN:INT (@', channels, ')']);
    fprintf(obj6517,':ROUT:SCAN:STIM 0');  	% settling time = x.x sec
    fprintf(obj6517,':ROUT:SCAN:SMET CURR');   % to scan CURR/VOLT
    fprintf(obj6517,':ROUT:SCAN:VSL 0');      	% 1: Sets V-source limit to 200V
    fprintf(obj6517,':ROUT:SCAN:LSEL INT'); % initiate scan
end

%% Apply voltage, reset time, and restart plotting
if speed == 1
    fprintf(obj6517,':SENS:CURR:NPLC 1');          % Integration rate line cycles (0.01-10)
    n_1s = 13;
elseif speed == 2
    fprintf(obj6517,':SENS:CURR:NPLC 0.1');
    n_1s = 26;
elseif speed == 3
    fprintf(obj6517,':SENS:CURR:NPLC 0.01');
    n_1s = 39;
end
Time = NaN(tfinal * n_1s, 1);   % Timestamp from Keithley
Time_k = NaN(tfinal * n_1s, 1);   % Timestamp from Matlab
Volt = NaN(tfinal * n_1s, 1);
CurrP = NaN(tfinal * n_1s, 1);
CurrN = NaN(tfinal * n_1s, 1);
Curr0V = NaN(tfinal * n_1s, 1);
dIdtP = NaN(tfinal * n_1s, 1);
dIdtN = NaN(tfinal * n_1s, 1);
dIdt0V = NaN(tfinal * n_1s, 1);
IBeg = NaN(tfinal * n_1s, 1);
IEnd = NaN(tfinal * n_1s, 1);
ZetaP = NaN(tfinal * n_1s, 1);
ZetaN = NaN(tfinal * n_1s, 1);
RNum = NaN(tfinal * n_1s, 1);   % Keithley reading number
ETemp = NaN(tfinal * n_1s, 1);

D1pre = NaN;
ZTnow = NaN;
tnow = 0;
tStartCycle = tnow;
n = 1;
c_maxt = 0;                 % add x * 10,000 s to time (each time Keithley timestamp resets to 0)
RNum_last = NaN;
ave_runSlope = NaN;

E_mult = 1;
d1_dt = 1;        % 1st deriv smoothing: d1(n) = I(n - d1_dt) to I(n + d1_dt)
compr = t2_skip - t1_skip;          % Compare to pt x sec's behind tnow for endpoint determination
% Update these each cycle
n1 = n + round(t1_skip * n_1s);
n2 = n + round(t2_skip * n_1s);
if d1_smooth > 0
    d1_dn = round(d1_smooth);
else
    d1_dn = 1;
end
nStartCycle = n;
tNextCycle = tCycle1;
found = 1;
ratio = NaN;
compr = floor(compr*n_1s);

% CREAte plot starting line 282
axI = subplot(3,1,1); % top left subplot
axI.FontSize = 10;
axI.FontName = 'Calibri';
hIP = plot(axI,Time(1:n),CurrP(1:n),'r-');
hold on;
hIN = plot(axI,Time(1:n),CurrN(1:n),'k-');
hI0V = plot(axI,Time(1:n),Curr0V(1:n),'g-');
hIEnd = plot(axI,Time(1:n),IEnd(1:n),'b*');
ylabel(axI,'|Curr| (nA)','FontName','Calibri','FontSize',14,'FontWeight','Bold');
if process == 2
    legend(axI,strcat('+',num2str(vout1),'V'),strcat(num2str(vout2),'V'))
elseif process == 3
    legend(axI,strcat('+',num2str(vout1),'V'),strcat(num2str(vout2),'V'))
end

axd1I = subplot(3,1,2);
axd1I.FontSize = 10;
axd1I.FontName = 'Calibri';
hdIdtP = plot(axd1I,Time(1:n),dIdtP(1:n),'r-');
hold on;
hdIdtN = plot(axd1I,Time(1:n),dIdtN(1:n),'k-');
hdIdt0V = plot(axd1I,Time(1:n),dIdt0V(1:n),'g-');
ylabel(axd1I,'dIdt (nA/s)','FontName','Calibri','FontSize',14,'FontWeight','Bold');
% ylim([-0.5 0.5])
if process == 2
    legend(axI,strcat('+',num2str(vout1),'V'),strcat(num2str(vout2),'V'))
elseif process == 3
    legend(axI,strcat('+',num2str(vout1),'V'),strcat(num2str(vout2),'V'))
end

axZ = subplot(3,1,3);
axZ.FontSize = 10;
axZ.FontName = 'Calibri';
hZetaP = plot(axZ,Time(1:n),ZetaP(1:n),'ro');
hold on;
hZetaN = plot(axZ,Time(1:n),ZetaN(1:n),'ko');
xlabel(axZ,'Time (s)','FontName','Calibri','FontName','Calibri','FontSize',14,'FontWeight','Bold')
ylabel(axZ,'Zeta (mV)','FontName','Calibri','FontSize',14,'FontWeight','Bold')
legend(strcat('+',num2str(vout),'V'),strcat('-',num2str(vout),'V'))

sgt = sgtitle(name);
sgt.FontSize = 12;
sgt.FontWeight = 'Normal';
sgt.FontName = 'Calibri';

drawnow;
% CREAte plot starting line 282

% To Keithley: apply voltage, measure current
fprintf(obj6517,[':SOUR:VOLT ', num2str(vout)]);
fprintf(obj6517,'OUTP ON');            % Start Voltage
timer=tic;                          % Start Matlab timer
fprintf(obj6517,':INIT:CONT ON');              % Start reading current
pause(1.3);
while tnow <= tfinal
    fprintf(obj6517,':FETC?');    % :FETCh? / :DATA:FRESh? / TRACe:DATA? -> get last / new reading / all readings
    tnow = toc(timer);
    tempI = fscanf(obj6517);       % ASCii string in order: READing,TSTamp,CHANnel,VSOurce
    tempI = str2num(tempI);
    if n > 1 && isnan(Time_k(n-1)) == 0 && tempI(1,2) < Time_k(n-1)     % Passed 10,000 s Keith timestamp reset point
        c_maxt = c_maxt + 10000;
    end
    if n > 1 && isnan(RNum(n-1)) == 0
        RNum_last = RNum(n-1);
    end
    while isnan(RNum_last) == 0 && tempI(1,3) == RNum_last
        fprintf(obj6517,':FETC?');    % :FETCh? / :DATA:FRESh? / TRACe:DATA? -> get last / new reading / all readings
        tnow = toc(timer);
        tempI = fscanf(obj6517);       % ASCii string in order: READing,TSTamp,RNUMber,CHANnel,VSOurce
        tempI = str2num(tempI);
        if tempI(1,2) < Time_k(n-1)     % Passed 10,000 s Keith timestamp reset point
            c_maxt = c_maxt + 10000;
        end
    end
    Time_k(n) = tempI(1,2) + c_maxt;
    Inow = 1e9*abs(tempI(1,1));
    RNum(n) = tempI(1,3);
    ETemp(n) = tempI(1,4);
    
    % 1st Derivative
    if n >= 2*d1_dn + 1
        if Volt(n-d1_dn) > 0
            D1pre = (Inow - CurrP(n - 2*d1_dn)) ...   % D1pre is d1_dn behind D1now
            / (tnow - Time(n - 2*d1_dn));
            dIdtP(n-d1_dn) = D1pre;
        elseif Volt(n-d1_dn) < 0
            D1pre = (Inow - CurrN(n - 2*d1_dn)) ...   % D1pre is d1_dn behind D1now
            / (tnow - Time(n - 2*d1_dn));
            dIdtN(n-d1_dn) = D1pre;
        else
            D1pre = (Inow - Curr0V(n - 2*d1_dn)) ...   % D1pre is d1_dn behind D1now
            / (tnow - Time(n - 2*d1_dn));
            dIdt0V(n-d1_dn) = D1pre;
        end
    end

    if found == 0 && n >= n2 + d1_dn
    % Calc quotient of current slope over slope x sec's behind
        if vout > 0
            P1_dIdt = mean(dIdtP(n1 + d1_dn : n - d1_dn - compr));
            ratio = abs(D1pre)/abs(P1_dIdt);
        else
            P1_dIdt = mean(dIdtN(n1 + d1_dn : n - d1_dn - compr));
            ratio = abs(D1pre)/abs(P1_dIdt);
        end
        % Endpt found at n-d1_dn if: D1pre <= d1_thres * d1I_dt(n1)
        % MODIFIED LINE BELOW ON 9/29/22, check for 4% current change
        if abs(Inow-IStartCycle)/IStartCycle >= minCurrChange && ratio <= abs(d1_thres)   % && abs(P1_dIdt) >= minSlope 
            found = 1;
            tEnd = Time(n-d1_dn);
            t_EOF = tEnd - tStartCycle;
            signSlope = round(P1_dIdt/abs(P1_dIdt)); % 1 if positive slope, -1 if neg.
            ZTnow = 1e3 * signSlope * L^2 * mu / eps / ...
                vout / t_EOF;
            if process == 4
                t_PlatHold = abs(net_CV_cycle*(L^2 * mu / ((ZTnow - zeta_EP)/1000) / eps / vout)) - t_EOF; % desired net CVs * L/net_vel - time_traverse = t_remaining... L/net_vel = sec/CV
                if t_EOF + t_PlatHold > tCycle % if cycle time > tCycle
                    t_PlatHold = tCycle - t_EOF; % every cycle has a max time of tCycle, otherwise infinite cycle times in stagnant regime
                elseif t_PlatHold < 5
                    t_PlatHold = 5;
                end
            end
            if autoSwitch
                tNextCycle = tEnd + t_PlatHold;
            end
%             if process == 3 && a <= length(hold_V) && tEnd >= hold_V(a)
%                 vout = vout_hold;
%                 tNextCycle = tnow + hold_t(a);
%                 fprintf(obj6517,[':SOUR:VOLT ' num2str(vout)]);
%                 a=a+1;
%             end
            if vout > 0
                ZetaP(n-d1_dn) = ZTnow;
                IEndnow = CurrP(n-d1_dn);
            else
                ZetaN(n-d1_dn) = ZTnow;
                IEndnow = CurrN(n-d1_dn);
            end
            IEnd(n-d1_dn) = IEndnow;
            % plot last cycle
            set(hIP,'xdata',Time(1:n),'ydata',CurrP(1:n));
            set(hIN,'xdata',Time(1:n),'ydata',CurrN(1:n));
            set(hI0V,'xdata',Time(1:n),'ydata',Curr0V(1:n));
            set(hdIdtP,'xdata',Time(1:n),'ydata',dIdtP(1:n));
            set(hdIdtN,'xdata',Time(1:n),'ydata',dIdtN(1:n));
            set(hdIdt0V,'xdata',Time(1:n),'ydata',dIdt0V(1:n));
            set(hIEnd,'xdata',Time(1:n),'ydata',IEnd(1:n));
            set(hZetaP,'xdata',Time(1:n),'ydata',ZetaP(1:n));
            set(hZetaN,'xdata',Time(1:n),'ydata',ZetaN(1:n));
            drawnow;
        end
    end
    % Switch Voltage at Cycle end
    if tnow >= tNextCycle
        if process == 1
            break;
        end
        if add0Vcyc
            fprintf(obj6517,':SOUR:VOLT 0');
            pause(5);
            tnow = toc(timer);
        end
        tThisCycle = tnow - tStartCycle;
        n_1s = (n - nStartCycle) / (tThisCycle);
        n1 = n + round(t1_skip * n_1s);
        n2 = n + round(t2_skip * n_1s);
        d1_dn = round(d1_dt * n_1s);
        
        nStartCycle = n;
        tStartCycle = tnow;
        IStartCycle = Inow;
        IBeg(n) = IStartCycle;
        ratio = NaN;

%       Reset X/Y limits
        if tnow > 600
            x_lim = [tnow-8.1*tThisCycle,tnow+2.1*tThisCycle];
            axI.XLim = x_lim;
            axd1I.XLim = x_lim;
        end
        % axZ.XLim = axI.XLim;

        % Switch voltage
        found = 0;
        if splitter
            fprintf(obj6517,':ROUT:SCAN:LSEL NONE'); % terminate scan
        end
        % added on 6/16/22 to have constant cycle times
        if process == 3 && a <= length(hold_V) && tnow >= hold_V(a) % && autoswitch == 0
            if vout_hold == 0    
                vout_last = vout;
            else
                vout_last = vout_hold;
            end
            found = 1;
            tNextCycle = tnow + hold_t(a);
            vout = vout_hold;
            a=a+1; % end addition 6/20/22
        elseif vout_last == vout1
            vout = vout2;
            vout_last = vout2;
            tNextCycle = tnow + tCycle2;
        else
            vout = vout1;
            vout_last = vout1;
            tNextCycle = tnow + tCycle;
        end
        fprintf(obj6517,[':SOUR:VOLT ' num2str(vout)]);
        pause(0.5);
        if splitter
            fprintf(obj6517,[':ROUT:SCAN:INT (@', channels, ')']);
        end
    else
        IBeg(n) = NaN;
    end

    % update table and fig
    Time(n) = tnow;
    Volt(n) = vout;
    if vout < 0
        CurrN(n) = Inow;
    elseif vout > 0
        CurrP(n) = Inow;
    else
        Curr0V(n) = Inow;
    end
    if process==1
        set(hIP,'xdata',Time(1:n),'ydata',CurrP(1:n));
        set(hdIdtP,'xdata',Time(1:n),'ydata',dIdtP(1:n));
        drawnow;
    end
    display(strcat('Time: ',32,num2str(tnow),32,...
        's, Next Cycle:',32,num2str(tNextCycle-tnow),32,...
        's, I:',32,num2str(Inow),32,'nA, dIdt:',32,num2str(D1pre),32,...
        'nA/s, Zeta:',32,num2str(ZTnow),'mV, Res:',32,num2str(n_1s),32,'pt/sec'));
    n=n+1;
    if liveplot
        set(hIP,'xdata',Time(1:n),'ydata',CurrP(1:n));
        set(hIN,'xdata',Time(1:n),'ydata',CurrN(1:n));
        set(hI0V,'xdata',Time(1:n),'ydata',Curr0V(1:n));
        set(hdIdtP,'xdata',Time(1:n),'ydata',dIdtP(1:n));
        set(hdIdtN,'xdata',Time(1:n),'ydata',dIdtN(1:n));
        set(hdIdt0V,'xdata',Time(1:n),'ydata',dIdt0V(1:n));
        set(hIEnd,'xdata',Time(1:n),'ydata',IEnd(1:n));
        set(hZetaP,'xdata',Time(1:n),'ydata',ZetaP(1:n));
        set(hZetaN,'xdata',Time(1:n),'ydata',ZetaN(1:n));
        drawnow;
    end
end

%% Copy code starting here and paste in CMD line (if stopped early)
% if splitter
%     fprintf(obj6517,':ROUT:SCAN:LSET NONE');
% end

fprintf(obj6517,':OUTP OFF');                              % Turn off voltage
fprintf(obj6517, '*RST');                                  % Reset GPIB to default
fclose(obj6517);

n = n-1;
% Data reorganization
Time = Time(1:n);
Time_k = Time_k(1:n);
Volt = Volt(1:n);
CurrP = CurrP(1:n);
CurrN = CurrN(1:n);
Curr0V = Curr0V(1:n);
dIdtP = dIdtP(1:n);
dIdtN = dIdtN(1:n);
dIdt0V = dIdt0V(1:n);
IBeg = IBeg(1:n);
IEnd = IEnd(1:n);
ZetaP = ZetaP(1:n);
ZetaN = ZetaN(1:n);
ETemp = ETemp(1:n);
RNum = RNum(1:n);

set(hIP,'xdata',Time,'ydata',CurrP);
set(hIN,'xdata',Time,'ydata',CurrN);
set(hI0V,'xdata',Time,'ydata',Curr0V);
set(hdIdtP,'xdata',Time,'ydata',dIdtP);
set(hdIdtN,'xdata',Time,'ydata',dIdtN);
set(hdIdt0V,'xdata',Time,'ydata',dIdt0V);
set(hIEnd,'xdata',Time,'ydata',IEnd);
set(hZetaP,'xdata',Time,'ydata',ZetaP);
set(hZetaN,'xdata',Time,'ydata',ZetaN);
drawnow;

x_lim = [0,tnow];
axI.XLim = x_lim;
axd1I.XLim = x_lim;
axZ.XLim = x_lim;

[chan,type,buf,elec,dateFilt,GPIB,res,pH,buf_mM,elec_mM] = ...
    convertCharsToStrings(chan,type,buf,elec,dateFilt,GPIB,res,pH,buf_mM,elec_mM);
if process == 1
    Info = table(chan,type,d,L,res,pH,sig,eps,mu,buf_mM,buf,elec_mM,elec,dateFilt,...
    GPIB,process,speed,I_range,I_autorange);
elseif process == 2 || process == 3
    Info = table(chan,tCycle,tCycle1,t_PlatHold,type,d,L,res,pH,sig,eps,mu,buf_mM,buf,elec_mM,elec,dateFilt,...
    GPIB,process,speed,I_range,I_autorange);
elseif process == 4
Info = table(chan,net_CV_cycle,zeta_EP,tCycle,tCycle1,type,d,L,res,pH,sig,eps,mu,buf_mM,buf,elec_mM,elec,dateFilt,...
    GPIB,process,speed,I_range,I_autorange);
end
Info.chan = char(Info.chan);
Info.type = char(Info.type);
Info.buf = char(Info.buf);
Info.elec = char(Info.elec);
Info.GPIB = char(Info.GPIB);
Info.dateFilt = char(Info.dateFilt);
Info.pH = char(Info.pH);
Info.buf_mM = char(Info.buf_mM);
Info.elec_mM = char(Info.elec_mM);
if process == 3
	Info.hold_V = hold_V;
    Info.hold_t = hold_t;
%     Info.mins_bias_aza_switch = t_hold;
    Info.vout_hold = vout_hold;
end

Data = table(Time,Volt,CurrP,CurrN,dIdtP,dIdtN,IBeg,IEnd,ZetaP,ZetaN,ETemp);%);
% Save / email file
c = clock;
cs = sprintf('%0.2i-%0.2i-%0.4i_%0.2i-%0.2i-%0.2i', ...
   c(2), c(3), c(1), c(4), c(5), round(c(6)));
path_fig = strcat(folder,'\',cs,'.png');
path_dat = strcat(folder,'\',cs,'.mat');
saveas(gcf, path_fig);
save(path_dat,'Data','Info');


%% Notice
% mail = 'ucsbnanolab@gmail.com'; 
% password = '100mMBorate';
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail',mail);
% setpref('Internet','SMTP_Username',mail);
% setpref('Internet','SMTP_Password',password);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% sendmail(yourmail,strcat('Your experiment "',name,'" has finished'),...
%     'Hello! The result of your experiment is attached',{path_fig,path_dat})
