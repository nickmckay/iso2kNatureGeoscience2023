% Relate global temperature to precipitation d18O in iCESM, see if it makes
% any sense in comparison with Iso2k
% September 2019
% Sam Stevenson
% Update January 2022: implement area-weighted average for global values

runnames={'b.ie12.B1850CN.f19_g16.001','b.ie12.B1850C5CN.f19_g16.LME.002','b.ie12.B1850C5CN.f19_g16.LME.003'};
%runnames={'b.ie12.B1850C5CN.f19_g16.LME.GHG.001','b.ie12.B1850C5CN.f19_g16.LME.Orbital.001','b.ie12.B1850C5CN.f19_g16.LME.SSI_VSK_L.001', ...
%    'b.ie12.B1850C5CN.f19_g16.LME.VOLC_GRA.002'};
%runnames={'b.ie12.B1850C5CN.f19_g16.LME.VOLC_GRA.001'};
%runnames={'b.ie12.B1850C5CN.f19_g16.LME.GHG-20th.003','b.ie12.B1850C5CN.f19_g16.LME.GHG-20th.004'};
%runnames={'b.ie12.B1850C5CN.f19_g16.LME.GHG.004'};
%runnames={'b.ie12.B1850C5CN.f19_g16.LME.GHG.003', ...
%    'b.ie12.B1850C5CN.f19_g16.LME.O3AER.001','b.ie12.B1850C5CN.f19_g16.LME.O3AER.003', ...
%    'b.ie12.B1850C5CN.f19_g16.LME.O3AER.004','b.ie12.B1850C5CN.f19_g16.LME.O3AER.005'};
precvars={'PRECRC_H216Or','PRECRL_H216OR','PRECSC_H216Os','PRECSL_H216OS','PRECRC_H218Or','PRECRL_H218OR','PRECSC_H218Os','PRECSL_H218OS'};
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];
regbox=[-90,90,0,360];
refper=[1961,1990];
windlen=30;

rootdir='/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/';
%rootdir='/glade/campaign/univ/ucsb0006/CESM-CAM5-LME/';
strtyr=850;

% Coordinates
unc=netcdf(strcat(rootdir,'atm/proc/tseries/monthly/TS/',runnames{1},'.cam.h0.TS.185001-200512.nc'));
lat=unc{'lat'}(:);
lon=unc{'lon'}(:);
mylat=find(lat >= regbox(1) & lat <= regbox(2));
mylon=find(lon >= regbox(3) & lon <= regbox(4));
myw=find(lon < 0);
mye=find(lon >= 0);
gw=unc{'gw'}(:);
latwgt=gw./mean(gw);


% Locations of Iso2k proxies
% "original" list from 2022
%proxylats=[-84,76.617,-77.78,-75.1,-75,-3.08,-3.08,-3.08,-3.08,60.59,65.18,72.3,75.1,-77.32,80.7,80.7,80.7,71.27,80.7,75.1,80.7,71.27, ...
%    71.12,72.58,78.86,78.4,-72.82,-66.77,78.0006,78.8333,75.33,72.6,67.25,60.4,68.3,-10.7,39.8,48.9,47.1,-13.9333,44.8333,47.09,25.28, ...
%    32.083,30.45,-5.94,33.3,41.416,20.75,46.95,47.09,16.2086,17.4,-12.37,32.1,-11.27,17.4,36.57,-19.4,58.1969];
%proxylons=[43,-36.4033,158.72,123.39,0.01,37.35,37.35,37.35,37.35,-140.5,-43.83,-37.4,-42.32,39.7,-73.1,-73.1,-73.1,-26.73,-73.1,-42.32, ...
%    -73.1,-26.73,-37.32,-37.64,17.42,-80.4,159.18,112.81,-36.3978,-36.5,-82.5,-38.5,-66.75,-134.8,18.7,-76.06,-107.3,-117.3,11.0166, ...
%    -70.8333,-92.25,11.67,108.08,-105.1667,110.416,-77.31,105,31.934,-89.47,10.55,11.67,-89.0735,-99.2,-41.57,104.26,-75.79,-99.2, ...
%    -118.78,17.883,25.6275];

% updated list from Georgy 4/2023: P_isotope records with end date >= 850
proxylats=[65.21,-84,45.9244,-74.97,-74.5,-75,-74.85,-75,-75,-74.96,-75.08,-75.17, ...
    -75.25,-75.17,-74.75,-74.67,-75.08,-75.25,-74.21,73.9402,76.617,80,77.2533, ...
    76.0039,-73.1,-77.78,28.38,28.38,-75.1,-75,-75.58,-75.17,81,-79.57,-70.62, ...
    -72.8,-3.08,-16.62,27.98,60.59,79.83,-70.24,38.1,65.18,72.3,75.1,35.81,-77.32, ...
    28.03,78.86,-77.32,-73.59,80.7,71.27,-46.58,49.8067,80.7,75.1,80.52,80.7,77.17, ...
    71.27,-71.34,71.12,65.18,72.6,72.58,70.3,71.27,70.63,70.64,71.76,71.15,-77.33, ...
    78.86,78.4,-70.86,43.05,-72.82,-64.2,-77.52,-72.9,45.8411,45.9297,80.52,-66.77, ...
    -79.38,-78.43,-78.08,-77.68,-78.33,-82,-78.12,-77.06,-82,-83.5,-86.5,-80.62, ...
    -73.72,-74.88,-89.93,-74.57,-13.9333,-78.47,46.55,-67.43,29.04,77.45,-66.04, ...
    75.2504,78.0006,78.8333,79.3414,78,76.6594,76.6594,75.0016,76,-79.36,-69.95, ...
    -77.11,-76.41,-72.81,43.3483,75.33,75.33,46.8842,66.38,60.58,72.59,-70.68,80.7, ...
    -9,-9,-73.6,-75.92,72.6,35.2833,67.25,-18.1,-74.65,-70.25,-74.86,-74.97,-74.49, ...
    -74.41,-75,-75,-74.59,-75.75,-75.93,-75.22,46.5481,72.4,72.7,69.23,67.95,60.4, ...
    8.85,68.3,60.95,54.35,63.12,-10.7,-10.0183,-3.3167,39.8,48.9,70.459,65.6107, ...
    -8.0089,34.7778,0.0833,38.86,38.4397,6.5,37.9603,68.6374,37,45.9166,79.7733, ...
    58.1969,69.0763,47.1,33.05,-5.5883,35.5,35.5,-3.8833,14.2704,13.3322,32.58, ...
    44.8333,-27.22,47.09,25.28,32.083,38.9,17.17,-27.22,30.45,37.976,-5.94,33.3, ...
    41.416,19,20.75,46.95,30.92,42.098,47.09,16.2086,17.4,-12.37,32.1,-8.533,-5.93, ...
    -11.27,17.4,36.57,-15.5,-19.4,-5.92,-5.92,43.47,43.47,10,36.62,46.35,46.5,68.1, ...
    52.5,53.2833,10.2,-12.6,-22,30.3083,19.08,19.08,10.43,50.23,38.8,-11.4,16.65, ...
    40.875,41.7883,41.4567,68.4,27.983,49,-10.0833,54.56,52.2217,51.8399,29.45, ...
    57.815,-50.5167];

proxylons=[-138.32,43,7.8675,3.92,1.96,0.04,-8.5,-6.5,-4.51,-1.5,2.5,5,-6,-1, ...
    1,4,6.5,6.5,-9.75,-37.6299,-36.4033,-41.1374,-49.2167,-43.492,165.4,158.72,85.72, ...
    85.72,123.39,0.01,-3.43,6.5,64,-45.72,-8.37,159.06,37.35,-67.77,86.92,-140.5, ...
    24,4.8,96.4,-43.83,-37.4,-42.32,90.76,39.7,86.96,17.42,39.7,-70.36,-73.1,-26.73, ...
    -73.32,86.56,-73.1,-42.32,94.82,-73.1,-61.13,-26.73,11.59,-37.32,-43.83,-38.5, ...
    -37.64,-44.58,-26.73,-35.82,-39.62,-35.85,-35.84,162.53,17.42,-80.4,11.54,94.32, ...
    159.18,-57.69,167.68,169.08,6.8478,7.877,94.82,112.81,-111.24,-111.92,-120.08,-124, ...
    -124.48,-110.01,-95.65,-89.14,-110.01,-104.99,-107.99,-122.63,7.94,1.6,144.39,-86.9, ...
    -70.8333,106.83,8.07,93.38,90.2,-51.06,-64.08,-37.6248,-36.3978,-36.5,-45.9116,-44, ...
    -46.4837,-46.4837,-42.0004,-46,-161.7,95.62,95.07,102.17,79.93,42.4267,-82.5,-82.5, ...
    10.8256,-46.18,-140.58,-38.46,-64.87,-73.1,-77.5,-77.5,-12.43,-84.25,-38.5,81.4833, ...
    -66.75,-68.9,12.8,4.82,-2.55,3.92,1.97,7.22,0.01,8.01,-3.44,3.29,7.22,11.35,8.0461, ...
    126,143.5,86.57,32.48,-134.8,-70.87,18.7,-148.15,-6.69,12.32,-76.06,34.1878,37.7, ...
    -107.3,-117.3,-70.086,-37.6935,113.3128,-120.0392,37.5333,73.26,75.0572,-1.4167,48.5553, ...
    -50.98,100,8.8225,10.7378,25.6275,-146.9297,11.0166,5,11.2217,9.5,9.5,119.45,39.4478, ...
    39.3646,35.03,-92.25,-49.16,11.67,108.08,-105.1667,-92.3,54.3,-49.16,110.416,-80.4,-77.31, ...
    105,31.934,82,-89.47,10.55,90.07,-123.4072,11.67,-89.0735,-99.2,-41.57,104.26,120.433, ...
    -77.3,-75.79,-99.2,-118.78,167,17.883,-77.35,-77.35,-91.97,-91.97,-85,74.98,8.6,8.77,60, ...
    -118,107.6333,-85.35,-69.2,-66,91.5167,82.33,82.33,76.93,89.04,-105,-68.716,-95,-124.0683, ...
    -124.0767,-124.0467,-133.8,90,86,-66.3,-71.2,-4.228,-4.1515,96.43,15.26,-70.1167];

figure(1)
clf
for rr=1:length(runnames)
    runname=runnames{rr}
    
    % Read in temperature
    atmfile=strcat(rootdir,'atm/proc/tseries/monthly/TS/',runname,'.cam.h0.TS.085001-184912.nc');
    ts=nc_varget(atmfile,'TS',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)])-273.15;
    %ts=[];
    atmfile=strcat(rootdir,'atm/proc/tseries/monthly/TS/',runname,'.cam.h0.TS.185001-200512.nc');
    ts=cat(1,ts,nc_varget(atmfile,'TS',[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)])-273.15);
        
%     % Plot map of mean temperature to make sure iCESM isn't crazy
%     tsmn=squeeze(nanmean(ts,1));
%     figure(1)
%     clf
%     cla
%     m_proj('mercator','lon',[regbox(3) regbox(4)],'lat',[regbox(1) regbox(2)])
%     ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
%     m_pcolor(lon,lat(mylat),tsmn);
%     shading flat
%     hold all
%     m_coast('color',[0 0 0]);
%     m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
%     set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
%     title('Mean temperature: 850-2005','FontSize',20)
%     colorbar
%     strcat('/glade/scratch/samantha/iCESM/plots/meanTS_850-1850_',runname,'.fig')
%     saveas(gcf,strcat('/glade/scratch/samantha/iCESM/plots/meanTS_850-1850_',runname,'.fig'))

    % Read in precip variables, calculate precip d18O    
    if rr == 1
        time1=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H216Or/',runname,'.cam.h0.PRECRC_H216Or.085001-184912.nc'),'time');
        %time1=[];
        time2=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H216Or/',runname,'.cam.h0.PRECRC_H216Or.185001-200512.nc'),'time');
        %atmtime=cat(1,time1,time2);
        atmtime=cat(1,time1,time2+365000.);
        %atmtime=time2;
        [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[strtyr 1 1]);
        uyr=unique(pyr);
        min(uyr)
        max(uyr)
    
        tsarr=zeros(length(runnames),length(uyr));
        d18oparr=zeros(length(runnames),length(uyr));
        proxyts=zeros(length(runnames),length(proxylats),length(uyr));
        proxyd18op=zeros(length(runnames),length(proxylats),length(uyr));
    end
    
    % Calculate total precipitation
    vtmp1=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H216Or/',runname,'.cam.h0.PRECRC_H216Or.085001-184912.nc'),precvars{1},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    %vtmp1=[];
    vtmp2=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H216Or/',runname,'.cam.h0.PRECRC_H216Or.185001-200512.nc'),precvars{1},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    PREC_H216O=cat(1,vtmp1,vtmp2);
    vtmp1=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H218Or/',runname,'.cam.h0.PRECRC_H218Or.085001-184912.nc'),precvars{5},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    %vtmp1=[];
    vtmp2=nc_varget(strcat(rootdir,'atm/proc/tseries/monthly/PRECRC_H218Or/',runname,'.cam.h0.PRECRC_H218Or.185001-200512.nc'),precvars{5},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
    PREC_H218O=cat(1,vtmp1,vtmp2);

    for ii=2:4
        atmfile=strcat(rootdir,'atm/proc/tseries/monthly/',precvars{ii},'/',runname,'.cam.h0.',precvars{ii},'.085001-184912.nc');
        vtmp1=nc_varget(atmfile,precvars{ii},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
        %vtmp1=[];
        atmfile=strcat(rootdir,'atm/proc/tseries/monthly/',precvars{ii},'/',runname,'.cam.h0.',precvars{ii},'.185001-200512.nc');
        vtmp2=nc_varget(atmfile,precvars{ii},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
        PREC_H216O=PREC_H216O+cat(1,vtmp1,vtmp2);
    end
    for ii=6:8
        atmfile=strcat(rootdir,'atm/proc/tseries/monthly/',precvars{ii},'/',runname,'.cam.h0.',precvars{ii},'.085001-184912.nc');
        vtmp1=nc_varget(atmfile,precvars{ii},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
        %vtmp1=[];
        atmfile=strcat(rootdir,'atm/proc/tseries/monthly/',precvars{ii},'/',runname,'.cam.h0.',precvars{ii},'.185001-200512.nc');
        vtmp2=nc_varget(atmfile,precvars{ii},[0 min(mylat)-1 min(mylon)-1],[-1 length(mylat) length(mylon)]);
        PREC_H218O=PREC_H218O+cat(1,vtmp1,vtmp2);
    end
    
    % Calculate precip d18O
    d18op=1000*(PREC_H218O./PREC_H216O - 1);
    
    % Annual mean
    pannmn=zeros(length(uyr),length(mylat),length(mylon));
    d18oannmn=zeros(length(uyr),length(mylat),length(mylon));
    tsannmn=zeros(length(uyr),length(mylat),length(mylon));
    for yy=1:length(uyr)
        myy=find(pyr == uyr(yy));
        ptot=PREC_H218O(myy,:,:)+PREC_H216O(myy,:,:);
        
        tsannmn(yy,:,:)=squeeze(nanmean(ts(myy,:,:),1));
        pannmn(yy,:,:)=squeeze(nanmean(ptot,1));
        d18oannmn(yy,:,:)=nanmean(squeeze(d18op(myy,:,:)).*ptot./repmat(pannmn(yy,:,:),[12 1 1]),1);    % precip-weighting following Vachon et al. (2007)
    end
    
    % Subset over proxy locations
    for pp=1:length(proxylats)
       prlat=find(lat >= proxylats(pp)-1 & lat <= proxylats(pp)+1);
       prlon=find(lon >= proxylons(pp)-1 & lon <= proxylons(pp)+1);
       proxyts(rr,pp,:)=squeeze(nanmean(nanmean(tsannmn(:,prlat,prlon),2),3));
       proxyd18op(rr,pp,:)=squeeze(nanmean(nanmean(d18oannmn(:,prlat,prlon),2),3));
    end
        
    % Take global average
%     tsarr(rr,:)=squeeze(nanmean(nanmean(tsannmn,2),3));
%     d18oparr(rr,:)=squeeze(nanmean(nanmean(d18oannmn,2),3));

    % Calculate matrix of weights
    wgtfac=repmat(latwgt(mylat)',[size(tsannmn,1) 1 size(tsannmn,3)]);

    % Area-weighted global average
    tsarr(rr,:)=squeeze(nanmean(nanmean(tsannmn.*wgtfac,2),3));
    d18oparr(rr,:)=squeeze(nanmean(nanmean(d18oannmn.*wgtfac,2),3));    

    % Save data in text files
    % Global
    savets=cat(2,uyr',squeeze(tsarr(rr,:))');
    save(strcat('/glade/scratch/samantha/globalTtimeseries_areawgt_',runnames{rr},'.txt'),'savets','-ASCII')
    savets=cat(2,uyr',squeeze(d18oparr(rr,:))');
    save(strcat('/glade/scratch/samantha/globald18Optimeseries_areawgt_',runnames{rr},'.txt'),'savets','-ASCII')
    
    % Iso2k proxy site average
    ts_sitetmp=squeeze(nanmean(proxyts(rr,:,:),2));
    d18op_sitetmp=squeeze(nanmean(proxyd18op(rr,:,:),2));
    savets=cat(2,uyr',ts_sitetmp);
    save(strcat('/glade/scratch/samantha/Iso2ksiteTtimeseries_',runnames{rr},'.txt'),'savets','-ASCII')
    savets=cat(2,uyr',d18op_sitetmp);
    save(strcat('/glade/scratch/samantha/Iso2ksited18Optimeseries_',runnames{rr},'.txt'),'savets','-ASCII')

end

% Subtract ensemble mean over reference period
% myref=find(uyr >= refper(1) & uyr <= refper(2));
% proxyts=proxyts-nanmean(nanmean(proxyts(:,:,myref),1),3);
% proxyd18op=proxyd18op-nanmean(nanmean(proxyd18op(:,:,myref,:),1),3);
% tsannmn=tsannmn-nanmean(tsannmn(myref,:,:),1);
% d18oannmn=d18oannmn-nanmean(d18oannmn(myref,:,:),1);
% tsarr=tsarr-nanmean(tsarr(:,myref),2);
% d18oparr=d18oparr-nanmean(d18oparr(:,myref),2);
% globtsmn=squeeze(nanmean(nanmean(tsarr,2),3));
% globd18opmn=squeeze(nanmean(nanmean(d18oparr,2),3));

tsrunmn=zeros(size(tsarr));
d18oprunmn=zeros(size(d18oparr));
proxytsrunmn=zeros(size(proxyts));
proxyd18oprunmn=zeros(size(proxyd18op));
for rr=1:length(runnames)
    tsrunmn(rr,:)=runmean(tsarr(rr,:),floor(windlen/2));
    d18oprunmn(rr,:)=runmean(d18oparr(rr,:),floor(windlen/2));
    for pp=1:length(proxylats)
        proxytsrunmn(rr,pp,:)=runmean(proxyts(rr,pp,:),floor(windlen/2));
        proxyd18oprunmn(rr,pp,:)=runmean(proxyd18op(rr,pp,:),floor(windlen/2));
    end
end

% Determine relationship between T, d18Op pre- and post-industrial
PIyr=find(uyr > 860 & uyr < 1850);
poyr=find(uyr >= 1850 & uyr <= 1995);
% PIcoef=zeros(2,length(runnames),length(proxylats));
% pocoef=zeros(2,length(runnames),length(proxylats));
% for rr=1:length(runnames)
%     for pp=1:length(proxylats)
%     PIcoef(:,rr,pp)=regress(squeeze(proxyd18oprunmn(rr,pp,PIyr)),cat(2,ones(length(PIyr),1),squeeze(proxytsrunmn(rr,pp,PIyr))));
%     pocoef(:,rr,pp)=regress(squeeze(proxyd18oprunmn(rr,pp,poyr)),cat(2,ones(length(poyr),1),squeeze(proxytsrunmn(rr,pp,poyr))));
%     end
% end
% squeeze(nanmean(nanmean(PIcoef(2,:,:),2),3))
% squeeze(nanmean(nanmean(pocoef(2,:,:),2),3))

proxytsensmn=squeeze(nanmean(nanmean(proxytsrunmn,1),2));
globtsensmn=nanmean(tsrunmn,1);
proxyd18oensmn=squeeze(nanmean(nanmean(proxyd18oprunmn,1),2));
globd18oensmn=nanmean(d18oprunmn,1);

PIcoef_proxy=regress(proxyd18oensmn(PIyr),cat(2,ones(length(PIyr),1),proxytsensmn(PIyr)));
PIcoef_glob=regress(globd18oensmn(PIyr)',cat(2,ones(length(PIyr),1),globtsensmn(PIyr)'));
pocoef_proxy=regress(proxyd18oensmn(poyr),cat(2,ones(length(poyr),1),proxytsensmn(poyr)));
pocoef_glob=regress(globd18oensmn(poyr)',cat(2,ones(length(poyr),1),globtsensmn(poyr)'));


% Plot results
good=find(uyr > 860 & uyr <= 1995);
h=[];
figure(1)
clf
yyaxis left
h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxytsrunmn(:,:,good),1),2)),'LineWidth',2,'Color','k')];
hold all
h=[h,plot(uyr(good),nanmean(tsrunmn(:,good),1),'--','LineWidth',2,'Color','k')];
yyaxis right
h=[h,plot(uyr(good),squeeze(nanmean(nanmean(proxyd18oprunmn(:,:,good),1),2)),'LineWidth',2,'Color','r')];
h=[h,plot(uyr(good),nanmean(d18oprunmn(:,good),1),'-.','LineWidth',2,'Color','r')];
set(gca,'FontSize',24)
title('Temperature, precip \delta^{18}O')
saveas(gcf,strcat('/glade/scratch/samantha/iCESM/plots/proxyloc_plusglobal_td18op_',num2str(windlen),'yrrunmntimeseries.fig'),'fig')
