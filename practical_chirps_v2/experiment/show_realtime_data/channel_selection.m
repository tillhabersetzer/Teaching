function [channel_out] = channel_selection(desired,datachannel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the desired channel group elements which are actually
% contained in datachannel
% - datachannel contains the actual channel-labels in the data
% - desired is a cell containing element-wise the desired channels
%   'E' : all EEG 
%   'M' : all MEG 
%   'M1': all Magnetometer
%   'M2': all Gradiometer with 2
%   'M3': all Gradiometer with 3
%   'MG': all Gradiometer (2 and 3)
%   'EF','ET','EP',EO','EC': EEG Frontal, Temporal, Parietal,Occipital,
%                                                              Central
%   'ELF','ERF',......,'ELO','ERO: EEG Left/Right Temporal,...Occipital
%   'MF','MT','MP',MO': MEG Frontal, Temporal, Parietal, Occipital                       
%   'MLF','MRF',......,'MLO','MRO: MEG Left/Right Temporal,...Occipital
%
% and there is a specific Option to plot only specified MEG sensortypes:
%    -M1: only Magnetometer
%    -M2: only Gradiometer with 2
%    -M3: only Gradiometer with 3
%    -MG: only Gradiometer
%    e.g. MLF-M1 : MEG Left Frontal, only Magnetometers
%
% so all in all it could be something like this: 
% desired = {'ELT','ELP','MLF-M1,....};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define channel groups
%--------------------------------------------------------------------------

EEG = {'EEG001';'EEG002';'EEG003';'EEG004';'EEG005';'EEG006';'EEG007';'EEG008';'EEG009';'EEG010';...
       'EEG011';'EEG012';'EEG013';'EEG014';'EEG015';'EEG016';'EEG017';'EEG018';'EEG019';'EEG020';...
       'EEG021';'EEG022';'EEG023';'EEG024';'EEG025';'EEG026';'EEG027';'EEG028';'EEG029';'EEG030';...
       'EEG031';'EEG032'};
   
MAGNETOMETER = ...
       {'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';... 
        'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';...
        'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';...
        'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';...
        'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';...
        'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';...
        'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';...
        'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';...
        'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';...
        'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';...
        'MEG2631';'MEG2641'};
    
GRADIOMETER2 = ...
        {'MEG0112';'MEG0122';'MEG0132';'MEG0142';'MEG0212';'MEG0222';'MEG0232';'MEG0242';'MEG0312';'MEG0322';...
         'MEG0332';'MEG0342';'MEG0412';'MEG0422';'MEG0432';'MEG0442';'MEG0512';'MEG0522';'MEG0532';'MEG0542';...
         'MEG0612';'MEG0622';'MEG0632';'MEG0642';'MEG0712';'MEG0722';'MEG0732';'MEG0742';'MEG0812';'MEG0822';...
         'MEG0912';'MEG0922';'MEG0932';'MEG0942';'MEG1012';'MEG1022';'MEG1032';'MEG1042';'MEG1112';'MEG1122';...
         'MEG1132';'MEG1142';'MEG1212';'MEG1222';'MEG1232';'MEG1242';'MEG1312';'MEG1322';'MEG1332';'MEG1342';...
         'MEG1412';'MEG1422';'MEG1432';'MEG1442';'MEG1512';'MEG1522';'MEG1532';'MEG1542';'MEG1612';'MEG1622';...
         'MEG1632';'MEG1642';'MEG1712';'MEG1722';'MEG1732';'MEG1742';'MEG1812';'MEG1822';'MEG1832';'MEG1842';...
         'MEG1912';'MEG1922';'MEG1932';'MEG1942';'MEG2012';'MEG2022';'MEG2032';'MEG2042';'MEG2112';'MEG2122';...
         'MEG2132';'MEG2142';'MEG2212';'MEG2222';'MEG2232';'MEG2242';'MEG2312';'MEG2322';'MEG2332';'MEG2342';...
         'MEG2412';'MEG2422';'MEG2432';'MEG2442';'MEG2512';'MEG2522';'MEG2532';'MEG2542';'MEG2612';'MEG2622';...
         'MEG2632';'MEG2642'};
                
GRADIOMETER3 = ...
        {'MEG0113';'MEG0123';'MEG0133';'MEG0143';'MEG0213';'MEG0223';'MEG0233';'MEG0243';'MEG0313';'MEG0323';...
         'MEG0333';'MEG0343';'MEG0413';'MEG0423';'MEG0433';'MEG0443';'MEG0513';'MEG0523';'MEG0533';'MEG0543';...
         'MEG0613';'MEG0623';'MEG0633';'MEG0643';'MEG0713';'MEG0723';'MEG0733';'MEG0743';'MEG0813';'MEG0823';...
         'MEG0913';'MEG0923';'MEG0933';'MEG0943';'MEG1013';'MEG1023';'MEG1033';'MEG1043';'MEG1113';'MEG1123';...
         'MEG1133';'MEG1143';'MEG1213';'MEG1223';'MEG1233';'MEG1243';'MEG1313';'MEG1323';'MEG1333';'MEG1343';...
         'MEG1413';'MEG1423';'MEG1433';'MEG1443';'MEG1513';'MEG1523';'MEG1533';'MEG1543';'MEG1613';'MEG1623';...
         'MEG1633';'MEG1643';'MEG1713';'MEG1723';'MEG1733';'MEG1743';'MEG1813';'MEG1823';'MEG1833';'MEG1843';...
         'MEG1913';'MEG1923';'MEG1933';'MEG1943';'MEG2013';'MEG2023';'MEG2033';'MEG2043';'MEG2113';'MEG2123';...
         'MEG2133';'MEG2143';'MEG2213';'MEG2223';'MEG2233';'MEG2243';'MEG2313';'MEG2323';'MEG2333';'MEG2343';...
         'MEG2413';'MEG2423';'MEG2433';'MEG2443';'MEG2513';'MEG2523';'MEG2533';'MEG2543';'MEG2613';'MEG2623';...
         'MEG2633';'MEG2643'};
     
GRADIOMETER = [GRADIOMETER2;GRADIOMETER3];
MEG         = [MAGNETOMETER;GRADIOMETER2;GRADIOMETER3];
 
% you can also adress Sensor-Triplets
% one row is one Sensortriplett
SENSORTRIPLETT = [MAGNETOMETER,GRADIOMETER2,GRADIOMETER3];
 
% now arrange in topograhical order
%-------------------------------------------------------------------------- 
% EEG
%----
EEG_LeftFrontal    = {'EEG001';'EEG003';'EEG004';'EEG021'};
EEG_RightFrontal   = {'EEG002';'EEG006';'EEG007';'EEG024'};
EEG_LeftTemporal   = {'EEG021';'EEG008';'EEG009';'EEG025';'EEG013'};
EEG_RightTemporal  = {'EEG024';'EEG011';'EEG012';'EEG028';'EEG017'};
EEG_LeftParietal   = {'EEG025';'EEG026';'EEG013';'EEG014';'EEG029'};
EEG_RightParietal  = {'EEG027';'EEG028';'EEG016';'EEG017';'EEG031'};
EEG_LeftOccipital  = {'EEG029';'EEG018'};
EEG_RightOccipital = {'EEG031';'EEG019'};
EEG_Frontal        = [EEG_LeftFrontal;EEG_RightFrontal;{'EEG020';'EEG005';'EEG022';'EEG023'}];
EEG_Temporal       = [EEG_LeftTemporal;EEG_RightTemporal];
EEG_Parietal       = [EEG_LeftParietal;EEG_RightParietal;{'EEG015';'EEG030'}];
EEG_Occipital      = [EEG_LeftOccipital;EEG_RightOccipital;{'EEG030';'EEG032'}];
EEG_Central        = {'EEG020';'EEG005';'EEG021';'EEG022';'EEG023';'EEG024';'EEG009';'EEG010';'EEG011';'EEG025';...
                      'EEG026';'EEG027';'EEG028';'EEG015';'EEG030';'EEG032'};
 
% MEG
%----
% define channels with channel-tripletts by using the magnetometer as a
% landmark for the triplett

% magnetometer landmark
MEG_LeftFrontal_mag    = {'MEG0121';'MEG0341';'MEG0321';'MEG0331';'MEG0641';'MEG0621';'MEG0311';'MEG0541';'MEG0611';'MEG0511';...
                          'MEG0531';'MEG0821';'MEG0521'};
MEG_RightFrontal_mag   = {'MEG1031';'MEG1241';'MEG1231';'MEG1221';'MEG1411';'MEG1011';'MEG1021';'MEG0931';'MEG1211';'MEG0941';...
                          'MEG0921';'MEG0811';'MEG0911'};
MEG_LeftTemporal_mag   = {'MEG0111';'MEG0131';'MEG0211';'MEG0221';'MEG0141';'MEG1511';'MEG0241';'MEG0231';'MEG1541';'MEG1521';...
                          'MEG1611';'MEG1621';'MEG1531'};
MEG_RightTemporal_mag  = {'MEG1421';'MEG1311';'MEG1321';'MEG1441';'MEG1431';'MEG1341';'MEG1331';'MEG2611';'MEG2621';'MEG2411';...
                          'MEG2421';'MEG2641';'MEG2631'};
MEG_LeftParietal_mag   = {'MEG0411';'MEG0421';'MEG0631';'MEG0441';'MEG0431';'MEG0711';'MEG1811';'MEG1821';'MEG0741';'MEG1631';...
                          'MEG1841';'MEG1831';'MEG2011'};
MEG_RightParietal_mag  = {'MEG1041';'MEG1111';'MEG1121';'MEG0721';'MEG1141';'MEG1131';'MEG0731';'MEG2211';'MEG2221';'MEG2241';...
                          'MEG2231';'MEG2441';'MEG2021'};
MEG_LeftOccipital_mag  = {'MEG1721';'MEG1641';'MEG1711';'MEG1731';'MEG1941';'MEG1911';'MEG1741';'MEG1931';'MEG1921';'MEG2041';...
                          'MEG2141';'MEG2111'};
MEG_RightOccipital_mag = {'MEG2121';'MEG2131';'MEG2031';'MEG2341';'MEG2331';'MEG2541';'MEG2311';'MEG2321';'MEG2511';'MEG2531';...
                          'MEG2431';'MEG2521'};

% logial index values
lf = contains(SENSORTRIPLETT(:,1),MEG_LeftFrontal_mag);
rf = contains(SENSORTRIPLETT(:,1),MEG_RightFrontal_mag);
lt = contains(SENSORTRIPLETT(:,1),MEG_LeftTemporal_mag);
rt = contains(SENSORTRIPLETT(:,1),MEG_RightTemporal_mag);
lp = contains(SENSORTRIPLETT(:,1),MEG_LeftParietal_mag);
rp = contains(SENSORTRIPLETT(:,1),MEG_RightParietal_mag);
lo = contains(SENSORTRIPLETT(:,1),MEG_LeftOccipital_mag);
ro = contains(SENSORTRIPLETT(:,1),MEG_RightOccipital_mag);

% reshape to column
MEG_LeftFrontal    = reshape(SENSORTRIPLETT(lf,:),[],1);
MEG_RightFrontal   = reshape(SENSORTRIPLETT(rf,:),[],1);
MEG_LeftTemporal   = reshape(SENSORTRIPLETT(lt,:),[],1);
MEG_RightTemporal  = reshape(SENSORTRIPLETT(rt,:),[],1);
MEG_LeftParietal   = reshape(SENSORTRIPLETT(lp,:),[],1);
MEG_RightParietal  = reshape(SENSORTRIPLETT(rp,:),[],1);
MEG_LeftOccipital  = reshape(SENSORTRIPLETT(lo,:),[],1);
MEG_RightOccipital = reshape(SENSORTRIPLETT(ro,:),[],1);
MEG_Frontal        = [MEG_LeftFrontal;MEG_RightFrontal];
MEG_Temporal       = [MEG_LeftTemporal;MEG_RightTemporal];
MEG_Parietal       = [MEG_LeftParietal;MEG_RightParietal]; 
MEG_Occipital      = [MEG_LeftOccipital;MEG_RightOccipital]; 

%% 
% desired = {'M';'T';'L'};
% sensor type: M: MEG, M1: Magnetomer, M2: Gradiometer, M3: Gradiometer, 
% E: EEG
% topography: F: Fromtal, T: Temporal, P: Parietal, O: Occipital,
% C: Central only for EEG
% side: L: left, R: right, LR: 

channel = [];
for n = 1:length(desired)
    
    element      = desired{n};
    check_sensor = 0;
    % check if a specified sensortyp is requested e.g. MLT-M1
    if contains(element,'-')
        pair = strsplit(element,'-');
        element    = pair{1};
        sensortype = pair{2}; % specified sensor type
        
        switch sensortype
            case 'M1'
                sensor = MAGNETOMETER;
            case 'M2'
                sensor = GRADIOMETER2;
            case 'M3'
                sensor = GRADIOMETER3;
             case 'MG'
                sensor = GRADIOMETER;
            otherwise 
                error(['sensortype ' sensor{1} 'is not supported']);
        end
        % respone to apply sensor filterung
        check_sensor = 1;
    end
    
    switch element
        
        case 'E'
            channel = [channel;EEG];
        case 'M' 
            channel = [channel;MEG];
        case 'M1'
            channel = [channel;MAGNETOMETER];
        case 'M2'
            channel = [channel;GRADIOMETER2];
        case 'M3'
            channel = [channel;GRADIOMETER3];
        case 'MG' 
            channel = [channel;GRADIOMETER];
        case 'ELF' 
            channel = [channel;EEG_LeftFrontal];
        case 'ERF' 
            channel = [channel;EEG_RightFrontal];
        case 'ELT' 
            channel = [channel;EEG_LeftTemporal];
        case 'ERT' 
            channel = [channel;EEG_RightTemporal];
        case 'ELP' 
            channel = [channel;EEG_LeftParietal];
        case 'ERP' 
            channel = [channel;EEG_RightParietal];
        case 'ELO' 
            channel = [channel;EEG_LeftOccipital];
        case 'ERO' 
            channel = [channel;EEG_RightOccipital];   
        case 'EF' 
            channel = [channel;EEG_Frontal];
        case 'ET' 
            channel = [channel;EEG_Temporal];    
        case 'EP'
            channel = [channel;EEG_Parietal];
        case 'EO' 
            channel = [channel;EEG_Occipital];    
        case 'EC' 
            channel = [channel;EEG_Central];
        case 'MLF' 
            channel = [channel;MEG_LeftFrontal];
        case 'MRF' 
            channel = [channel;MEG_RightFrontal];
        case 'MLT' 
            channel = [channel;MEG_LeftTemporal];
        case 'MRT' 
            channel = [channel;MEG_RightTemporal];
        case 'MLP' 
            channel = [channel;MEG_LeftParietal];
        case 'MRP' 
            channel = [channel;MEG_RightParietal];
        case 'MLO' 
            channel = [channel;MEG_LeftOccipital];
        case 'MRO' 
            channel = [channel;MEG_RightOccipital];    
        case 'MF' 
            channel = [channel;MEG_Frontal];
        case 'MT' 
            channel = [channel;MEG_Temporal];    
        case 'MP'
            channel = [channel;MEG_Parietal];
        case 'MO' 
            channel = [channel;MEG_Occipital];   
        otherwise
            error(['channel group ' element ' not found.']);
      
    end
    
     % if specific sensortyp has been requested, only these sensortype
     % should be considered
     if check_sensor
         channel = sensor(contains(sensor,channel));
     end
end

% return the elements which are embedded in datachannel
channel_out = datachannel(contains(datachannel,channel));

end
     