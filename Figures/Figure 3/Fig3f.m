% Plot time series of global precipitation in iCESM to help interpret
% synthesis paper results
% April 2022
% Sam Stevenson

runnames={'b.ie12.B1850CN.f19_g16.001','b.ie12.B1850C5CN.f19_g16.LME.002','b.ie12.B1850C5CN.f19_g16.LME.003'};
% runnames={'b.ie12.B1850CN.f19_g16.001'};
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];
startyr=850;
regbox=[-90,90,0,360];
refper=[1961,1990];
windlen=30;

% Coordinates
unc=netcdf(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECC/',runnames{1},'.cam.h0.PRECC.085001-184912.nc'));
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
    
    % Read in convective, large-scale precip
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECC/',runname,'.cam.h0.PRECC.085001-184912.nc');
    precc=nc_varget(atmfile,'PRECC',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECC/',runname,'.cam.h0.PRECC.185001-200512.nc');
    precc=cat(1,precc,nc_varget(atmfile,'PRECC',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]));

    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECL/',runname,'.cam.h0.PRECL.085001-184912.nc');
    precl=nc_varget(atmfile,'PRECL',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    atmfile=strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECL/',runname,'.cam.h0.PRECL.185001-200512.nc');
    precl=cat(1,precl,nc_varget(atmfile,'PRECL',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]));
    
    prect=(precc+precl)*1000*86400;
    
    % Assign output variables   
    if rr == 1
        time1=nc_varget(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECC/',runname,'.cam.h0.PRECC.085001-184912.nc'),'time');
        time2=nc_varget(strcat('/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/atm/proc/tseries/monthly/PRECC/',runname,'.cam.h0.PRECC.185001-200512.nc'),'time');
        atmtime=cat(1,time1,time2+365000.);
        [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[850 1 1]);
        uyr=unique(pyr);
    
        lharr=zeros(length(runnames),length(uyr));
    end
    
    % Annual mean
    lhannmn=zeros(length(uyr),length(mylat),length(mylon));
    for yy=1:length(uyr)
        myy=find(pyr == uyr(yy));
        
        lhannmn(yy,:,:)=squeeze(nanmean(prect(myy,:,:),1));
    end
        
    % Calculate matrix of weights
    wgtfac=repmat(latwgt',[size(lhannmn,1) 1 size(lhannmn,3)]);

    % Area-weighted global average
    lharr(rr,:)=squeeze(mean(mean(lhannmn.*wgtfac,2,'omitnan'),3,'omitnan'));
    
    % Save time series for later
    savets=cat(2,uyr',squeeze(lharr(rr,:))');
    save(strcat('/glade/scratch/samantha/precttimeseries_areawgt_iLMEmem',num2str(rr),'.txt'),'savets','-ASCII')
end

% Subtract ensemble mean over reference period
myref=find(uyr >= refper(1) & uyr <= refper(2));
lhannmn=lhannmn-nanmean(lhannmn(myref,:,:),1);
lharr=lharr-nanmean(lharr(:,myref),2);
globlhmn=squeeze(nanmean(nanmean(lharr,2),3));

lhrunmn=zeros(size(lharr));
for rr=1:length(runnames)
    lhrunmn(rr,:)=runmean(lharr(rr,:),floor(windlen/2));
end


% Plot results
good=find(uyr > 860 & uyr <= 1995);
h=[];
figure(1)
clf
%yyaxis left
%h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxylhrunmn(:,:,good),1),2)),'LineWidth',2,'Color','k')];
h=[h,plot(uyr(good),lhrunmn(:,good),'--','LineWidth',2,'Color','k')];
%hold all
% yyaxis right
% h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxyd18oprunmn(:,:,good),1),2)),'LineWidth',2,'Color','r')];
% h=[h,plot(uyr(good),nanmean(d18oprunmn(:,good),1),'-.','LineWidth',2,'Color','r')];
set(gca,'FontSize',24)
title('Precipitation (mm/day)')
saveas(gcf,strcat('/glade/scratch/samantha/iCESM/plots/global_prect_areawgt_',num2str(windlen),'yrrunmntimeseries.fig'),'fig')
