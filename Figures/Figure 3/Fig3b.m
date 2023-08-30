% Plot time series of global vertically integrated vapor d18O in iCESM to help interpret
% synthesis paper results
% October 2021
% Sam Stevenson

runnames={'b.ie12.B1850CN.f19_g16.001','b.ie12.B1850C5CN.f19_g16.LME.002','b.ie12.B1850C5CN.f19_g16.LME.003'};
% runnames={'b.ie12.B1850CN.f19_g16.001'};
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];
startyr=850;
regbox=[-90,90,0,360];
refper=[1961,1990];
windlen=30;

% Parameters for computing d18O
Rstd=2005.2e-6; % O18/O16 ratio in standard (here, VSMOW)
Mh2o=18.015;    % mean molecular mass of water
Mh218o=20.013;  % molecular mass of H218Oh
Mh216o=18.009;  % molecular mass of H216O

% Coordinates
unc=netcdf(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/LHFLX/',runnames{1},'.cam.h0.LHFLX.085001-184912.nc'));
lat=unc{'lat'}(:);
lon=unc{'lon'}(:);
mylat=find(lat >= regbox(1) & lat <= regbox(2));
mylon=find(lon >= regbox(3) & lon <= regbox(4));
myw=find(lon < 0);
mye=find(lon >= 0);
gw=unc{'gw'}(:);
latwgt=gw./mean(gw);

% Locations of Iso2k proxies
proxylats=[-84,76.617,-77.78,-75.1,-75,-3.08,-3.08,-3.08,-3.08,60.59,65.18,72.3,75.1,-77.32,80.7,80.7,80.7,71.27,80.7,75.1,80.7,71.27, ...
    71.12,72.58,78.86,78.4,-72.82,-66.77,78.0006,78.8333,75.33,72.6,67.25,60.4,68.3,-10.7,39.8,48.9,47.1,-13.9333,44.8333,47.09,25.28, ...
    32.083,30.45,-5.94,33.3,41.416,20.75,46.95,47.09,16.2086,17.4,-12.37,32.1,-11.27,17.4,36.57,-19.4,58.1969];
proxylons=[43,-36.4033,158.72,123.39,0.01,37.35,37.35,37.35,37.35,-140.5,-43.83,-37.4,-42.32,39.7,-73.1,-73.1,-73.1,-26.73,-73.1,-42.32, ...
    -73.1,-26.73,-37.32,-37.64,17.42,-80.4,159.18,112.81,-36.3978,-36.5,-82.5,-38.5,-66.75,-134.8,18.7,-76.06,-107.3,-117.3,11.0166, ...
    -70.8333,-92.25,11.67,108.08,-105.1667,110.416,-77.31,105,31.934,-89.47,10.55,11.67,-89.0735,-99.2,-41.57,104.26,-75.79,-99.2, ...
    -118.78,17.883,25.6275];

figure(1)
clf
for rr=1:length(runnames)
    runname=runnames{rr}
    
    % Read in vapor H216O, H218O
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/H216OV/',runname,'.cam.h0.H216OV.085001-184912.nc');
    ts=nc_varget(atmfile,'H216OV',[0 end-1 min(mylat)-1 min(mylon)-1],[-1 1 length(mylat) length(mylon)]);
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/H216OV/',runname,'.cam.h0.H216OV.185001-200512.nc');
    ts=cat(1,ts,nc_varget(atmfile,'H216OV',[0 end-1 min(mylat)-1 min(mylon)-1],[-1 1 length(mylat) length(mylon)]));
    
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/H218OV/',runname,'.cam.h0.H218OV.085001-184912.nc');
    ts18=nc_varget(atmfile,'H218OV',[0 0 min(mylat)-1 min(mylon)-1],[-1 -1 length(mylat) length(mylon)]);
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/H218OV/',runname,'.cam.h0.H218OV.185001-200512.nc');
    ts18=cat(1,ts18,nc_varget(atmfile,'H218OV',[0 end-1 min(mylat)-1 min(mylon)-1],[-1 1 length(mylat) length(mylon)]));

    % Vertically "integrate" H216OV, H218OV
    %h216ov_int=squeeze(nansum(ts,2));
    %h218ov_int=squeeze(nansum(ts18,2));    
    % Lowest model level
    h216ov_int=ts;
    h218ov_int=ts18;
    
    clear ts
    clear ts18
    
    % Compute d18O vapor
    h218ov_int=h218ov_int.*Rstd.*(Mh218o./Mh216o);
    vt=h218ov_int+h216ov_int;
    o16frac=h216ov_int./vt;
    o18frac=h218ov_int./vt;
    d18ov=1000*((Mh216o./Mh218o)*(o18frac./o16frac)-Rstd)./Rstd; 

    
    % Assign output variables   
    if rr == 1
        time1=nc_varget(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/LHFLX/',runname,'.cam.h0.LHFLX.085001-184912.nc'),'time');
        time2=nc_varget(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/LHFLX/',runname,'.cam.h0.LHFLX.185001-200512.nc'),'time');
        atmtime=cat(1,time1,time2+365000.);
        [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[850 1 1]);
        uyr=unique(pyr);
    
        d18ovarr=zeros(length(runnames),length(uyr));
    end
    
    % Annual mean
    d18ovannmn=zeros(length(uyr),length(mylat),length(mylon));
    for yy=1:length(uyr)
        myy=find(pyr == uyr(yy));
        
        d18ovannmn(yy,:,:)=squeeze(nanmean(d18ov(myy,:,:),1));
    end
    
%     % Subset over proxy locations
%     for pp=1:length(proxylats)
%        prlat=find(lat >= proxylats(pp)-1 & lat <= proxylats(pp)+1);
%        prlon=find(lon >= proxylons(pp)-1 & lon <= proxylons(pp)+1);
%        proxyts(rr,pp,:)=squeeze(nanmean(nanmean(d18ovannmn(:,prlat,prlon),2),3));
%        proxyd18op(rr,pp,:)=squeeze(nanmean(nanmean(d18oannmn(:,prlat,prlon),2),3));
%     end
        
    % Calculate matrix of weights
    wgtfac=repmat(latwgt',[size(d18ovannmn,1) 1 size(d18ovannmn,3)]);

    % Area-weighted global average
    d18ovarr(rr,:)=squeeze(nanmean(nanmean(d18ovannmn.*wgtfac,2),3));
    
    % Save time series for later
    savets=cat(2,uyr',squeeze(d18ovarr(rr,:))');
    save(strcat('/glade/scratch/samantha/d18ov_lowestlevel_timeseries_areawgt_iLMEmem',num2str(rr),'.txt'),'savets','-ASCII')
end

% Subtract ensemble mean over reference period
myref=find(uyr >= refper(1) & uyr <= refper(2));
d18ovannmn=d18ovannmn-nanmean(d18ovannmn(myref,:,:),1);
d18ovarr=d18ovarr-nanmean(d18ovarr(:,myref),2);
globlhmn=squeeze(nanmean(nanmean(d18ovarr,2),3));

d18ovrunmn=zeros(size(d18ovarr));
for rr=1:length(runnames)
    d18ovrunmn(rr,:)=runmean(d18ovarr(rr,:),floor(windlen/2));
end


% Plot results
good=find(uyr > 860 & uyr <= 1995);
h=[];
figure(1)
clf
%yyaxis left
%h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxyd18ovrunmn(:,:,good),1),2)),'LineWidth',2,'Color','k')];
h=[h,plot(uyr(good),d18ovrunmn(:,good),'--','LineWidth',2,'Color','k')];
%hold all
% yyaxis right
% h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxyd18oprunmn(:,:,good),1),2)),'LineWidth',2,'Color','r')];
% h=[h,plot(uyr(good),nanmean(d18oprunmn(:,good),1),'-.','LineWidth',2,'Color','r')];
set(gca,'FontSize',24)
title('d18O vapor')
saveas(gcf,strcat('/glade/scratch/samantha/iCESM/plots/d18ovapor_',num2str(windlen),'yrrunmntimeseries.fig'),'fig')
