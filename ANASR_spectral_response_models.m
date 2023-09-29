%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% FOR ANALYSING (ALL) THE DATA FILES THAT ARE OUTPUT OF APE_Isc.m FOR A GIVEN
% MODULE TYPE SPECIFY THE INPUT DIRECTORY WHERE THE FILES ANA_dat_input.txt
% AND ANA_txt_input.txt ARE. THESE FILES CONTAIN RESPECTIVELY THE .dat AND
% .txt FILES WHICH ARE THE OUTPUT OF APESR_Isc.m FOR EACH RUN OF THE SCRIPT.
%
%                            IMPORTANT:

% In FITTINGx mode these files contain the data for which the fits will be
% perforemed.
%
% In FORECAST mode these files contain the data for which the Isc forecasts
% will be performed from fit parameters calculated from other data. The data
% files with the fit parameters are read in Isc_Forecast.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               YEAR 2006
% For aSi_Tianjin:
%INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\aSi_Tianjin\2006\';
%
% For pcSi_HYUNDAI:
%INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\pcSi_Hyundai\2006\';
%
% For CdTe_FS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               YEAR 2007
% For aSi_Tianjin:
%INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\aSi_Tianjin\2007\';
%
% For pcSi_HYUNDAI:
%INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\pcSi_Hyundai\2007\';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                              YEARS 2006 AND 2007
% For aSi_Tianjin:
INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\aSi_Tianjin\';
%
% For pcSi_HYUNDAI:
%INPUTDIR='C:\DIEGO\SR_APE_Isc_SOLIS_PVPM1000C\APE_Isc\APESR_Isc\OUT_APE_Isc_m\pcSi_Hyundai\';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   DON'T FORGET TO SPECIFY THE MODULE TYPE
%                   ALSO!!!!!!!!!!!!!!!!!!!!!
%
%
%                    PLEASE CHOOSE THE MODE OF RUN:
%
% FIT MODE READS DATA OUTPUT FROM APESR_Isc.m AND PERFORMS FITS OF Isc VS G
% AND Isc VS APE*Guf (SEVERAL CURVES, ONE FOR EACH CLEAR SKY INDEX GROUP DEFINED IN filterCSI.m)
%
% FORECAST MODE READS DATA OUTPUT FROM APESR_Isc.m AND FROM ANASR.m IN FIT
% MODE IN ORDER TO CALCULATE Isc VALUES FOR THE DATA READ ACCORDING TO THE
% FIT PARAMETERS READ:
%
                        modus='FITTINGx'
                        %modus='FORECAST'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%               FITS FOR EVALUATING THE MAGNITUDE OF THE SPECTRAL EFFECTS
%               
%
% modus_speff = 1: For fitting Isc25/G vs SE and obtaining the correcton
% parameters.
                    %modus_speff = 1
% modus_speff = 2: For using the correction parameters for fitting Isc/G vs
% G, T, APE and UD. This is for evaluate the magnitude and behaviour of the
% spectral effects. The fits of the Isc models should be included in this
% ANASR.m run.
                    %modus_speff = 2
                    
                    modus_speff = -1 % No empirical correction at all
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File with the .dat files to be processed:  
%
% For fitting in FITTING mode or
% for applying previuosly calculated fits to data in FORECAST mode. The fit
% files and paths are specified in Isc_Forecasts.m.
%
file1=[INPUTDIR,'ANASR_dat_input_DATA_ALL_aSi.txt'];                                     
[files]=textread(file1,'%s','headerlines',2);                             
%                                                                                                   
% File with the .txt files with the location name and time difference 
% (local time-UTC) for each of the .dat files in file1:
%
file2=[INPUTDIR,'ANASR_txt_input_DATA_ALL_aSi.txt'];                   
[files2]=textread(file2,'%s','headerlines',2);                                         
%                                                                         %  
%[loca TimeZone]=textread(file2,'%s %5.2f');                              %  
%location=loca{1,1};                                                      %
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%% MODULE TYPE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         PLEASE SPECIFY THE MODULE TYPE: cSi, aSi or CdTe   !!!!!!!!!!!  %
%                                                                         %
%                                                                         % 
modtype='aSi-Tianjin';                                                   %          
%modtype='CdTe-FSxxxx';                                                   %
%modtype='cSi-Hyundai';    %POLYCRYSTALLINE SILICON                        %
%                                                                         %
disp(' ')                                                                 %
disp(['The module type is ',modtype])                                     %
disp(' ')                                                                 %
%                                                                         %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SPECIFYING THE COEFFICIENTS FOR THE DIFFERENT MODULES:                  %
                                                                          %  
if modtype=='aSi-Tianjin'                                                 %
    TIscf =  +9E-4;  %+0.09% / °C                                         %      
    TUocf = -2.8E-3; %-0.28% /°C                                          %
    
    % The coefficients for the 3rd degree polynomial used for empirical
    % correction of SE effects: 
    
    a = +3.17e-9;
    b = -4.22e-7;
    c = +1.91e-5;
    d = +1.85e-3; 
    
elseif modtype=='cSi-Hyundai';                                            %
    TIscf =   +8E-4; %+0.08% / °C 
    TUocf = -3.0E-3; %-0.30% /°C
    
    % The coefficients for the 3rd degree polynomial used for empirical
    % correction of SE effects: 
    
    a = +3.17e-9;
    b = -4.22e-7;
    c = +1.91e-5;
    d = +1.85e-3; 
    
    a = +1.33e-8;
    b = -1.56e-6;
    c = +6.02e-5;
    d = +7.14e-3; 
elseif modtype=='CdTe-FSxxxx'
    TIscf =   +4E-4; %+0.04% / °C
    TUocf = -2.9E-3; %-0.29% /°C
end                                                                       %  
%                                                                         %  
% END OF PAREMETERS SPECIFICATION                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% READING OUTPUT DATA FROM APESR_Isc.m:
% Checking whether the number of .dat files equals the number of .txt
% files. They must be equal given that for each .dat file the location and
% the time zone are specified. This arose from the need to process data
% of summer time together with winter time data.
if length(files)~=length(files2)
    disp(' ')
    disp('Error: The number of .dat files (data)is different from the number of .txt files (location & time zone)')
    disp('Sequence aborted!')
    disp(' ')
    return
end

% The .txt files are read and each is stored in an individual structure.
for i=1:length(files2)
    %[loca TZ]=textread(files2{i,1},'%s %5.2f'); 
    %info(i).TimeZone=TZ;
    
    [loca]=textread(files2{i,1},'%s'); 
    info(i).location=loca;
    
    if info(i).location{1}=='MC'
        info(i).TimeZone=1
    else
        disp(' ')
        disp('      UNKNOWN LOCATION ')
        disp(' ')
        return
    end
end

% Checking that location is the same in all the .txt files.
if length(files2)>1
for i=1:length(files2)-1
    i
    if strcmp(info(i).location, info(i+1).location); 
        disp('Locations matched!')
    else
        disp(' ')
        disp('ERROR: The location of the .txt are different! ')
        disp('Sequence aborted!')
        disp(' ')
        return
    end
end
end

% Unique location and TimeZone are then defined:                                      %                       
location=info(i).location{1,1};  
TimeZone=info(i).TimeZone;
%  
% LOCATION PARAMETERS:  Latitude and longitude in degrees.                %
if strcmp(location,'MC')                                                  %   
    lat=48.36; %deg                                                       %  
    long=10.89; %deg                                                      %
    %TimeZone=1; %Local Time = UTC + TimeZone in agreement with convention used in APE_Isc.m,
    % explicitly in APE_Isc_READ_Isc.m, when the measured data is read.   %
    disp(' ')                                                             %  
    %disp(['Location name : ',location])                                  % 
    disp(location)
    disp(['Latitude  [deg]: ',num2str(lat)])                              %  
    disp(['Longitude [deg]: ',num2str(long)])                             %  
    disp(' ')                                                             %
else                                                                      %  
    disp('Location not defined! ')                                        %      
    disp('Sequence aborted')                                              %      
    return                                                                %  
end                                                                       % 

% If the location are the same then the .dat files are read and stored in
% individual structures.
for i=1:length(files2)
    
    data=load(files{i,1});
    Dat(i).year=data(:,1);
    Dat(i).doy=data(:,2);
    Dat(i).utc=data(:,3);
    Dat(i).ipt=data(:,4); % The flag witht the interpolation type.
    Dat(i).T=data(:,5);
    Dat(i).G=data(:,6);
    Dat(i).Gsolis=data(:,7);
    Dat(i).Isc=data(:,8);
    Dat(i).Uoc=data(:,9);
    Dat(i).Impp=data(:,10);
    Dat(i).Umpp=data(:,11);
    Dat(i).Pmpp=data(:,12);
    Dat(i).cix=data(:,13);
    Dat(i).ape=data(:,14);
    Dat(i).apesr=data(:,15); 
    Dat(i).otherp=data(:,16);
    Dat(i).flux=data(:,17); 
    Dat(i).fluxsr=data(:,18);
    Dat(i).uf=data(:,19); 
    Dat(i).ud=data(:,20); 
    
    disp(['Loading data file ',num2str(i)])
    disp(['Reading Data File: ',files{i,1}])
    disp([num2str(size(data)),' data points'])
    %clear data
end

% STORING DATA IN A SINGLE STRUCTURE:

k=0; % counter for writing data of al files into single structure
%lastj=0;
j=0;
for i=1:length(Dat)
    %lastj=j;
    lastk=k;
    for j=1:length(Dat(i).year)
         % WRITING TO STRUCTURE WHICH WILL CONTAIN THE DATA FOR ALL THE
         % FILES
         %if i==1
         %    k=j;
         %elseif i>1
             k=lastk+j;
         %end
         D.year(k)=Dat(i).year(j);
         D.doy(k)=Dat(i).doy(j);
         D.utc(k)=Dat(i).utc(j);
         D.ipt(k)=Dat(i).ipt(j);
         D.T(k)=Dat(i).T(j);
         D.G(k)=Dat(i).G(j);
         D.Gsolis(k)=Dat(i).Gsolis(j);
         D.Isc(k)=Dat(i).Isc(j);
         D.Uoc(k)=Dat(i).Uoc(j);
         D.Impp(k)=Dat(i).Impp(j);
         D.Umpp(k)=Dat(i).Umpp(j);
         D.Pmpp(k)=Dat(i).Pmpp(j);
         D.cix(k)=Dat(i).cix(j);
         D.ape(k)=Dat(i).ape(j);
         D.apesr(k)=Dat(i).apesr(j); 
         D.otherp(k)=Dat(i).otherp(j); 
         D.flux(k)=Dat(i).flux(j); 
         D.fluxsr(k)=Dat(i).fluxsr(j); 
         D.uf(k)=Dat(i).uf(j); 
         D.ud(k)=Dat(i).ud(j); 
         
    end
end

clear Dat % All the data is now in structure D.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ana_plots %obsolete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATION OF SUN ELEVATION (SE), AND AZIMUTH (SA); AIR MASS (AM) AND OTHERS:

% AMC_DIN5034.m is the function for the calculation of the Sun Elevation 
% and AM, the Sun Azimuth and other related values. The angle values are 
% given in degrees.

for i=1:length(D.utc)
    
    [D.TST(i), D.SE(i), D.AM(i), D.SA(i)]=AMC_DIN5034(D.year(i),D.doy(i),D.utc(i),TimeZone,lat,long);
    
end

% TST = True Solar Time. Time on site relative to fixed stars (No time zones or other alterations).
% SE = Solar Elevation (deg)
% AM = Air Mass
% SA = Solar Azimuth (deg). Measured respect to north (geographical
% celestial direction ????) clockwise.

%outfile = 'OUT_ANA_m\check_SE.dat';
%fid = fopen(outfile,'w');
%data=[D.year; D.doy; D.utc; D.TST; D.SE; D.AM; D.SA]; 
%fprintf(fid,'%d %d %6.3f %6.3f %8.4f %6.3f %8.4f\n',data); 
%fclose(fid);
%disp('WRITING TO FILE:')
%disp(outfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       FILTERS OF SUN ELEVATION AND SHADOW ZONES:

% These filters must be applied in both modes, FITTING or FORECAST.

% "minse": Minimum Solar Elevation (deg) to process a data point. 
% Under minse reflection effects in module's glass are too distorting.  The
% 5° value includes also the horizon line calculated with the program
% "Horizon". Above minse there are no objects in the light path. Exception
% are treated with the function filter_shadow.m

minse = 5; %default value
%minse = 10;

A = filter_minSE_nullnegIsc(minse, D);
clear D % The filtered data is now in the structure A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter for excluding data under certain Gmeasured value:
%Gmin = 300;
%C = filter_G_low_cutoff(Gmin, A);
%clear A
%A = C;
%clear C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% About the Shadow filter:

% On 15.05.2007 the horizon line was calculated using the program
% "Horizon". It is seen that there are obstacles producing slight shadow up 
% to 5°. In the west, centered at 90° (about 80°-100°) there is a big
% obstacle, a bulding which shadows the direct beam up to 12° in 
% April/August and May/July. Therefore filter for all the data with sun
% azimuth 80°-100°. In the reference frame of AMC_DIN5034 (North:0°, 
% [0°-360°), increasing clockwise (east)) this corresponds to azimuth 
% [260°-280°]. Data with Sun Elevation under 14° for this azimuth zone will
% be discarded for all the times of the year. This is not a problem for the
% other momths of the year because either the path or the Sun does not go
% that far in the east or the Sun Elvation is higher for this azimuth zone
% (the case of June). 

SAso = 260;  % Start of Shadow zone
SAsf = 280;  % End of Shadow zone
SEsmin = 14; % "Horizon" angle for the shadow zone

B = filter_shadow(SAso, SAsf, SEsmin, A);

clear A
A = B;
clear B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               FILTER EXTREME TEMPERATURES
%
% Because there were errors in the measurement of T maybe due to a short
% circuit, these values must be dropped.

lowcut  = -20; %°C Data below is discarded.
highcut = +70; %°C Data above is discarded.

% These are the same values used in filterT.m called in FITTING mode. This
% filter (filter_extremeT) must be applied in both modes, FITTING and
% FORECAST.

B = filter_extremeT(lowcut, highcut, A);

clear A
A = B;
clear B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       OTHER IMPORTANT FIELDS:
%
% INCLUSION OF FIELDS "location", "latitude", "longitude", "modtype", "Guf" and "Gud" IN THE
% STRUCTURE A WITH ALL THE DATA POINTS.
% It is necessary to calculate Umpp because the PVPM software does not
% do it for all the points (no rule detectable at first view).

            A.location=location;
            A.latitude=lat;
            A.longitude=long;
            A.modtype=modtype;
            A.Guf=A.G.*A.uf;
            A.Gufsolis=A.Gsolis.*A.uf;
            A.Gud=A.G.*A.ud;
            A.Gudsolis=A.Gsolis.*A.ud;
            A.Umppca=A.Pmpp./A.Impp; 
            A.Gdiff=A.Gsolis-A.G;

% INCLUSION OF FIELD "filters" WITH A RECORD OF THE FILTERS APPLIED TO THE
% DATA.
%A.filters.SE_lt=minse; %lt=less than
%A.filters.Isc_lte=0;   %lte=less than or equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TEMPERATURE CORRECTION: Isc -> Isc25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   CALCULATION OF THE Isc AT STC TEMPERATURE OF 25°C

%           Isc(T2) = Isc(T1)*(1+Tcoeff*(T2-T1))

for i=1:length(A.Isc)
    
    A.Isc25(i) = A.Isc(i)*(1+TIscf*(25-A.T(i)));
    A.Uoc25(i) = A.Uoc(i)*(1+TUocf*(25-A.T(i)));
    %A.G(i)     = A.G(i)*(1+TIscf*(25-A.T(i))); % ??????????
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       APPLICATION OF EMPIRICAL CORRECTION FOR SUN ELEVATION EFFECTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULATION OF THE Isc AT STC AM OF 1.5 =  SE = 41.8°

%       [Isc/G](SE2) = [Isc/G](SE1)+[f(SE2)-f(SE1)]    
%           Isc(SE2) = Isc(SE1)+G[f(SE2)-f(SE1)]
%
%  f(SE) is the fit function Isc/G(SE), which tends to be linear. A 3rd
%  degree polynomial is used as fitting function for not ignoring the fast
%  decay at low SE values.

if modus_speff == 2
    disp('***************************************************************** ')
    disp('                  APPLYING EMPIRICAL CORRECTION ')
    disp('***************************************************************** ')
    SEref = 41.8; % deg
    for i=1:length(A.Isc)
        % The 3rd degree polynomila model:
        f1 = a*A.SE(i)^3 + b*A.SE(i)^2 + c*A.SE(i)+d; 
        f2 = a*SEref^3 + b*SEref^2 + c*SEref+d;
        A.Isc25(i) = A.Isc25(i)+A.G(i)*(f2-f1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the TITLE "tit" which will be used for the plots of the FIT and FORECAST sections.

YEARS(1).year=A.year(1);
YEARS(1).doy(1)=A.doy(1);
x=1;
for i=1:length(A.year)-1
    if A.year(i)~=A.year(i+1)
       YEARS(x).doy(x+1)=A.doy(i);
       YEARS(x+1).year=A.year(i+1);
       YEARS(x+1).doy(x)=A.doy(i+1);
       x=x+1;
       %A.doy(i)
       %YEARS(1)
       %YEARS(2)
    elseif i==length(A.year)-1 && x==1
        YEARS(x).doy(x+1)=A.doy(i+1);
        %A.doy(i)
        %A.doy(i+1)
        %YEARS(1)
    elseif i==length(A.year)-1 && x>1
        YEARS(x).doy(x)=A.doy(i+1);
        %A.doy(i)
        %A.doy(i+1)
        %YEARS(1)
        %YEARS(2)
    end
end

tit=['Location: ',A.location,'. Module: ',A.modtype,'. DATES: '];
for i=1:length(YEARS)
    jours='DOYs: ';
    for j=1:length(YEARS(i).doy)
        if j==length(YEARS(i).doy)
            jours=[jours,num2str(YEARS(i).doy(j))];
        else
            jours=[jours,num2str(YEARS(i).doy(j)),'-'];
        end
    end
    tit=[tit, '  ',num2str(YEARS(i).year),' ',jours];
end

%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALL TO SCRIPT WITH SEVERAL PLOTS OF ALL THE DATA POINTS

plots_before_after_filters(A, [pwd,'\OUT_ANA_m\All_Data\'], tit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if modus=='FITTINGx'
    Isc_Fits
    anasr_files_write
elseif modus=='FORECAST'
    Isc_Forecast
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15.06.2007
%
%               FITS FOR EVALUATING THE MAGNITUDE OF THE SPECTRAL EFFECTS
%
% modus_speff = 1: For fitting Isc25/G vs SE and obtaining the correcton
% parameters.
                    %modus_speff = 1;
% modus_speff = 2: For using the correction parameters for fitting Isc/G vs
% G, T, APE and UD. This is for evaluate the magnitude and behaviour of the
% spectral effects. The fits of the Isc models should be included in this
% ANASR.m run.
                    %modus_speff = 2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if modus_speff == 1
    
    % For Isc25/G:

    % Fit of Isc25/G vs SE is used for correcting the Isc25 values for SE. 
    % Forecasts turned off. Fits for filtering and plotting. In this way
    % plots of the data corrected only for temperature are corrected and
    % can be compared with the ones after applying the SE correction.
    % The fit output parameters are a, b, c, d (different for each module 
    % type) and are specified at the beginning of this script together with
    % the temperature coefficients of each module:
    DOfit_IscdG_var_POLY3(A, 'Isc25', '[A]', 'SunElevation',  'SE',  '[deg]',tit, [pwd,'\OUT_ANA_m\SE_correction\Fit_Isc25dG_SE'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif modus_speff == 2
     
    % For Isc25/G:
        
    %  The results of the fit Isc25/G vs SE are used in a second run of ANASR.m
    %  Then we have the Isc25 values corrected by temperature and by the
    %  empirical SE  correction. Then we perform the fits Isc25/G vs G, T, 
    %  APE and UD. These fits should help to recognize the magnitude of the
    %  spectral effects dur to clouds.

    DOfit_IscdG_var_LINEAL(A, 'Isc25', '[A]', 'G',  'G',  '[W/m^2]',tit, [pwd,'\OUT_ANA_m\SpEff\G'])
    DOfit_IscdG_var_LINEAL(A, 'Isc25', '[A]', 'T',  'T',  '[°C]',tit, [pwd,'\OUT_ANA_m\SpEff\G'])
    DOfit_IscdG_var_LINEAL(A, 'Isc25', '[A]', 'APE','ape','[eV]',  tit, [pwd, '\OUT_ANA_m\SpEff\G'])
    DOfit_IscdG_var_LINEAL(A, 'Isc25', '[A]', 'UD','ud',  '[1]',  tit,  [pwd, '\OUT_ANA_m\SpEff\G'])

% For Isc25/Gsolis:
    DOfit_IscdGsolis_var_LINEAL(A, 'Isc25', '[A]', 'Gsolis',  'Gsolis',  '[W/m^2]',tit, [pwd,'\OUT_ANA_m\SpEff\Gsolis'])
    DOfit_IscdGsolis_var_LINEAL(A, 'Isc25', '[A]', 'T',  'T',  '[°C]',tit, [pwd,'\OUT_ANA_m\SpEff\Gsolis'])
    DOfit_IscdGsolis_var_LINEAL(A, 'Isc25', '[A]', 'APE','ape','[eV]',  tit, [pwd, '\OUT_ANA_m\SpEff\Gsolis'])
    DOfit_IscdGsolis_var_LINEAL(A, 'Isc25', '[A]', 'UD','ud',  '[1]',  tit,  [pwd, '\OUT_ANA_m\SpEff\Gsolis'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

