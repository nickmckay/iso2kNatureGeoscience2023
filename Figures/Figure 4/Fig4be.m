% Code to make a PCA figure following the Iso2k conventions, look at
% associated SLP signatures: for 20th century ONLY
% July 2020
% Sam Stevenson

runnames={'b.ie12.B1850CN.f19_g16.001','b.ie12.B1850C5CN.f19_g16.LME.002','b.ie12.B1850C5CN.f19_g16.LME.003'};
% runnames={'b.ie12.B1850CN.f19_g16.001'};
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];
startyr=850;
refper=[1961,1990];
levthr=0.1;
yrint=[1850,2005];  % years to consider
avgrege=[-5,5,210,270];	% Averaging regions, following Tokinaga et al. (2012)
avgregw=[-5,5,90,150];

% Parameters for computing soil water d18O
Rstd=2005.2e-6; % O18/O16 ratio in standard (here, VSMOW)
Mh2o=18.015;    % mean molecular mass of water
Mh218o=20.013;  % molecular mass of H218O
Mh216o=18.009;  % molecular mass of H216O

% Coordinates
unc=netcdf(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runnames{1},'.clm2.h0.H2OSOI.185001-200512.nc'));
lat=unc{'lat'}(:);
lon=unc{'lon'}(:);
lev=unc{'levgrnd'}(:);
mylev=find(lev <= levthr);     % 0-1cm depth

binlen=3;
dodtr=0;        % 1 = detrend all data; 0 = detrend no data
pcnum=1;        % Which mode to use for regressions?
blklen=10;      % length of blocks for block bootstrap
nboot=100;      % number of bootstrap samples
ns=floor(blklen/binlen);
    
Fs=0.0278;              % sampling frequency, 3 years in units of 1/month
wlen=3;                 % window length, samples 
noverlap=1;             % over which to overlap windows


% Locations of Iso2k proxies
% From Matt Fischer: "Iso2k_20C_Sites"
% SITES FOR 1850-2004
% proxylats=[-3.2556,27.85,0.933,-8.2573,-1.5,-16.8167,-21.2378,-21.0333,-23.15,13.598,-16.8167,32.467,-16.8167,-4.1916,-0.13,30.6486,17.933, ...
%     27.1059,-15.94,19.287,7.2859,-28.4589,-12.087,-0.4084,29.4333,7.98,-15,25.3903,24.9167,-4.6062,-17.5,-22.48,-22.1,-28.4617,38.34,-10.7, ...
%     50.83,37.87,10.75,30.45,33.3,41.416,20.75,16.2086,17.4,17.4,36.57,-15.5,36.62,46.35,46.5,53.2833,-12.6,-22,30.3083,19.08,10.43,29.85, ...
%     19.9,38.8,21.67,27.983,54.56,52.2217,51.8399,45.7333,48.3833,29.45,30.42,31.15,32.2167,-50.5167,28.183,29.633];
% proxylons=[40.1433,34.32,173,115.5757,124.833,179.2333,200.1722,55.25,43.58,144.836,179.2333,295.3,179.2333,151.9772,98.52,295.0112,292.999, ...
%     142.1941,166.04,110.656,134.2503,113.749,96.875,268.766,34.9667,277.95,167,279.8285,279.25,55.4244,210.1667,166.45,153,113.7683,34.46, ...
%     283.94,243.61,240.84,295.3,110.416,105,31.934,270.53,270.9265,260.8,260.8,241.22,167,74.98,8.6,8.77,107.6333,290.8,294,91.5167,82.33, ...
%     76.93,81.93,101.2,255,104.1,90,288.8,-4.228,-4.1515,0.3,2.6667,96.43,95.07,97.033,77.2167,289.8833,85.183,79.85];
% proxyint=["Temperature","Temperature","Temperature","EffectiveMoisture","EffectiveMoisture","Temperature","Temperature","Temperature", ...
%     "Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature", ...
%     "EffectiveMoisture","EffectiveMoisture","EffectiveMoisture","Temperature","Temperature","Temperature","Temperature","EffectiveMoisture", ...
%     "Temperature","EffectiveMoisture","EffectiveMoisture","Temperature","Temperature","Temperature","Temperature","Temperature","EffectiveMoisture", ...
%     "P_isotope","EffectiveMoisture","EffectiveMoisture","Temperature","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope", ...
%     "P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope", ...
%     "EffectiveMoisture","EffectiveMoisture","P_isotope","EffectiveMoisture","P_isotope","P_isotope","P_isotope","P_isotope","EffectiveMoisture", ...
%     "EffectiveMoisture","P_isotope","EffectiveMoisture","EffectiveMoisture","EffectiveMoisture","P_isotope","EffectiveMoisture","EffectiveMoisture"];


% % SITES FOR 1900-2004
% proxylats=[-3.2556,27.85,-21.905,0.933,-5.217,-8.2573,-1.5,5.87,-16.8167,-21.2378,-3.6206,-21.0333,-23.15,13.598,-16.8167,-8.0167,-16.8167,-4.1916, ...
%     -0.13,30.6486,11.77,17.933,27.1059,-3.2,16.2,5.87,5.87,-15.94,19.287,-10.2,-10.2,7.2859,-28.4589,-23.3572,-20.2792,-12.087,19.7,1.4167,29.4333, ...
%     7.98,25.3903,24.9167,-4.6062,-17.5,-22.48,-0.533,-28.4617,10.2773,24.9167,38.34,-10.7,50.83,10.75,30.45,33.3,41.416,20.75,16.2086,17.4,17.4,36.57, ...
%     -15.5,23.67,23.92,36.62,46.35,46.5,53.2833,10.2,-12.6,-22,30.3083,19.08,10.43,29.85,50.23,19.9,38.8,-11.4,21.67,27.983,49,42.6411,-10.0833,54.56, ...
%     52.2217,51.8399,45.7333,48.3833,29.45,30.42,31.15,32.2167,-50.5167,28.183,29.633];
% proxylons=[40.1433,34.32,113.965,173,145.817,115.5757,124.833,197.87,179.2333,200.1722,143.6819,55.25,43.58,144.836,179.2333,39.5,179.2333,151.9772, ...
%     98.52,295.0112,293.25,292.999,142.1941,40.1,298.51,197.87,197.87,166.04,110.656,123.5167,123.5167,134.2503,113.749,43.6195,185.3577,96.875,279.94, ...
%     173.0333,34.9667,277.95,279.8285,279.25,55.4244,210.1667,166.45,166.9283,113.7683,250.7869,279.25,34.46,283.94,243.61,295.3,110.416,105,31.934, ...
%     270.53,270.9265,260.8,260.8,241.22,167,284.25,283.17,74.98,8.6,8.77,107.6333,274.65,290.8,294,91.5167,82.33,76.93,81.93,89.04,101.2,255,291.284, ...
%     104.1,90,86,1.0025,293.7,288.8,-4.228,-4.1515,0.3,2.6667,96.43,95.07,97.033,77.2167,289.8833,85.183,79.85];
% proxyint=["Temperature","Temperature","Temperature","Temperature","Temperature","EffectiveMoisture","EffectiveMoisture","Temperature","Temperature", ...
%     "Temperature","EffectiveMoisture","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","Temperature", ...
%     "Temperature","Temperature","Temperature","Temperature","Temperature","Temperature","EffectiveMoisture","Temperature","EffectiveMoisture", ...
%     "EffectiveMoisture","EffectiveMoisture","Temperature","EffectiveMoisture","Temperature","Temperature","Temperature","Temperature","Temperature", ...
%     "EffectiveMoisture","Temperature","EffectiveMoisture","EffectiveMoisture","EffectiveMoisture","Temperature","Temperature","Temperature","Temperature", ...
%     "Temperature","Temperature","EffectiveMoisture","EffectiveMoisture","P_isotope","EffectiveMoisture","Temperature","P_isotope","P_isotope","P_isotope", ...
%     "P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","Temperature","Temperature","P_isotope","P_isotope","P_isotope","P_isotope", ...
%     "P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","P_isotope","EffectiveMoisture","P_isotope","EffectiveMoisture","P_isotope","P_isotope", ...
%     "EffectiveMoisture","P_isotope","P_isotope","EffectiveMoisture","P_isotope","P_isotope","P_isotope","P_isotope","EffectiveMoisture","EffectiveMoisture", ...
%     "P_isotope","EffectiveMoisture","EffectiveMoisture","EffectiveMoisture","P_isotope","EffectiveMoisture","EffectiveMoisture"];

% useint = "EffectiveMoisture";
% inds=find(proxyint == useint);
% inds=inds([1,2,4,6:end]);       % exclude points which aren't land in iCLM
% proxylats=proxylats(inds);
% proxylons=proxylons(inds);
% proxyint=proxyint(inds);


% % From Georgy: 1850-2004, no-carbonate
% proxylats=[19.9,21.67,45.7333,48.3833,29.633,31.15,32.2167,29.85,28.183,30.42];
% proxylons=[101.2,104.1,0.3,2.6667,79.85,97.033,77.2167,81.93,85.183,95.07];
% inds=1:length(proxylats);
% useint="NoCarbonate";

% % From Georgy: 1850-2004, carbonate-only
% proxylats=[25.3903,24.9167,50.83,37.87,38.34,46,37.87,-8.2573,7.98,-15.94,-1.5,7.2859,19.287];
% proxylons=[-80.1715,-80.75,-116.39,-119.16,34.46,-94.7,-119.16,115.5757,-82.05,166.04,124.833,134.2503,110.656];
% inds=1:length(proxylats);
% useint="CarbonateOnly";

% % From Georgy: 1850-2004, Asian monsoon only
% proxylats=[19.9,21.67,29.633,31.15,32.2167,29.85,19.287,28.183,30.42];
% proxylons=[101.2,104.1,79.85,97.033,77.2167,81.93,110.656,85.183,95.07];
% useint="AsianMonsoonOnly";

% From Georgy: 1850-2004, no Asian monsoon
proxylats=[25.3903,24.9167,50.83,37.87,38.34,46,37.87,-8.2573,7.98,-15.94,-1.5,7.2859,45.7333,48.3833];
proxylons=[-80.1715,-80.75,-116.39,-119.16,34.46,-94.7,-119.16,115.5757,-82.05,166.04,124.833,134.2503,0.3,2.6667];
useint="NoAsianMonsoon";

proxylons(proxylons < 0)=proxylons(proxylons < 0)+360;

% From Matt F's "Methods" document, on steps for block bootstrap:
% Steps:
% i) Split the time series into regular 10-year blocks (see below)
% ii) Randomly choose the first block
% iii) Find blocks which are similar (in their mean) to the current block, and  randomly sample from those blocks
% iv)  The second block is the successor block to the randomly chosen block in iii).    
% v) Continue as in iii)

% In step i), 10-years was chosen because this preserves interannual variance within blocks, and decadal persistence between blocks.
% In step iii), the similar blocks are chosen by nearest neighbor using ranked means. For blocks which have no successor (e.g. if the time series has a gap), the successor block is randomly selected from all other blocks.



figure(1)
clf
for rr=1:length(runnames)
    runname=runnames{rr}
    
    % Read in SLP
    atmfile=strcat('/glade/scratch/samantha/iCESM/atm/PSL/',runname,'.cam.h0.PSL.185001-200512.nc');

    slp=nc_varget(atmfile,'PSL')/100.;
    if rr == 1

        atmtime=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.185001-200512.nc'),'time');
        [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[1850 1 1]);
        myt=find(pyr >= yrint(1) & pyr <= yrint(2));
        pyr=pyr(myt);
        
        uyr=unique(pyr);
        binyr=min(uyr):binlen:(max(uyr)-binlen+1);
        pbin=zeros(length(binyr),length(proxylats));
        d18obin=zeros(length(binyr),length(proxylats));
        slpbin=zeros(length(binyr),size(slp,2),size(slp,3));
        regslope=zeros(length(runnames),size(slp,2),size(slp,3));
        dslpregslope=zeros(length(runnames),size(slp,2),size(slp,3));
        Vs=zeros(length(runnames),min(15,length(proxylats)));
        ECts=zeros(length(runnames),length(binyr));
        EOFs=zeros(length(runnames),length(proxylats));
        ECtsboot=zeros(nboot,length(runnames),length(binyr));
        EOFsboot=zeros(nboot,length(runnames),length(proxylats));
        Vboot=zeros(nboot,length(runnames),min(15,length(proxylats)));
        
        dslpbin=zeros(length(binyr),1);
    end
    slp=slp(myt,:,:);
    
    if dodtr == 1
        slpdtr=zeros(size(slp));
        for la=1:size(slp,2)
            for lo=1:size(slp,3)
                slpdtr(:,la,lo)=detrend(slp(:,la,lo));
            end
        end
        slp=slpdtr;
    end
    
    % Compute Walker strength time series using the Tokinaga et al. metric:
    mylate=find(lat >= avgrege(1) & lat <= avgrege(2));
    mylone=find(lon >= avgrege(3) & lon <= avgrege(4));
    mylatw=find(lat >= avgregw(1) & lat <= avgregw(2));
    mylonw=find(lon >= avgregw(3) & lon <= avgregw(4));
    slpdiff=squeeze(nanmean(nanmean(slp(:,mylate,mylone),2),3))-squeeze(nanmean(nanmean(slp(:,mylatw,mylonw),2),3));
    
    % Annual mean
    dslpannmn=zeros(length(uyr),1);
    for yy=1:length(uyr)
        myy=find(pyr == uyr(yy));
        dslpannmn(yy)=squeeze(nanmean(slpdiff(myy)));
    end

    % Do binning
    for yy=1:length(binyr)
        myy=find(uyr >= binyr(yy) & uyr <= binyr(yy)+binlen);
        dslpbin(yy)=squeeze(nanmean(dslpannmn(myy)));
    end
    
    for pp=1:length(proxylats)
        mylat=find(lat >= proxylats(pp)-2 & lat <= proxylats(pp)+2);
        mylon=find(lon >= proxylons(pp)-2 & lon <= proxylons(pp)+2);
%         if pp == 2 
%             mylat=find(lat >= proxylats(pp)-3 & lat <= proxylats(pp)+3);
%             mylon=find(lon >= proxylons(pp)-3 & lon <= proxylons(pp)+3);
%         end

        % Read in variables, calculate soil d18O    
        h2osoi=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.185001-200512.nc'),'H2OSOI',[min(myt)-1 min(mylev)-1 min(mylat)-1 min(mylon)-1],[length(myt) length(mylev) length(mylat) length(mylon)]);
        h2osoi_h218o=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI_H218O/',runname,'.clm2.h0.H2OSOI_H218O.185001-200512.nc'),'H2OSOI_H218O',[min(myt)-1 min(mylev)-1 min(mylat)-1 min(mylon)-1],[length(myt) length(mylev) length(mylat) length(mylon)]);

        % Calculate soil d18O, make regional average
        d18osoil=squeeze(nanmean(nanmean(nanmean(1000*(h2osoi_h218o./h2osoi - 1),2),3),4));

        % Annual mean
        pannmn=zeros(length(uyr),1);
        d18oannmn=zeros(length(uyr),1);
        slpannmn=zeros(length(uyr),size(slp,2),size(slp,3));
        for yy=1:length(uyr)
            myy=find(pyr == uyr(yy));
            ptot=squeeze(nanmean(nanmean(h2osoi_h218o(myy,:,:)+h2osoi(myy,:,:),2),3));

            slpannmn(yy,:,:)=squeeze(nanmean(slp(myy,:,:),1));
            pannmn(yy,:)=squeeze(nanmean(ptot,1));
            d18oannmn(yy,:)=nanmean(squeeze(d18osoil(myy,:)).*ptot./repmat(pannmn(yy,:),[12 1 1]),1);    % precip-weighting following Vachon et al. (2007)
        end

        % Do binning
        for yy=1:length(binyr)
            myy=find(uyr >= binyr(yy) & uyr <= binyr(yy)+binlen);
            pbin(yy,pp)=squeeze(nanmean(pannmn(myy,:),1));
            slpbin(yy,:,:)=squeeze(nanmean(slpannmn(myy,:,:),1));
            d18obin(yy,pp)=squeeze(nanmean(d18oannmn(myy,:),1));
        end
        
        if dodtr == 1
            % Detrend d18O to correspond with what Matt F did
            d18obin(:,pp)=detrend(d18obin(:,pp));
        end
    end
    
    % SVD analysis
    d18obin(abs(d18obin) > 1e10)=0;
    d18obin(isnan(d18obin))=0;
    [V,EOFtmp,EC,error]=EOF(d18obin,15);
    V/sum(V)        % percent variance in each mode

    
    % Make sure signs match based on EOF loadings
    ltmp=find(abs(EOFtmp(:,pcnum)) == max(abs(EOFtmp(:,pcnum))));
    
    if EOFtmp(ltmp(1),pcnum) < 0
        EOFtmp(:,pcnum)=-1*EOFtmp(:,pcnum);
        EC(:,pcnum)=-1*EC(:,pcnum);
    end
    
    ECts(rr,:)=EC(:,pcnum);
    EOFs(rr,:)=EOFtmp(:,pcnum);
    Vs(rr,:)=V/sum(V);
    
    % Save PC, EOF data for later
%     savets=cat(2,binyr',EC(:,1:10));
    savets=cat(2,binyr',EC);
    save(strcat('/glade/scratch/samantha/PCtimeseries_iLMEmem',num2str(rr),'-',useint,'.txt'),'savets','-ASCII')

    savets=cat(2,binyr',dslpbin);
    save(strcat('/glade/scratch/samantha/dSLPtimeseries_iLMEmem',num2str(rr),'-',useint,'.txt'),'savets','-ASCII')
    
    EOFmat=EOFtmp;
    PCts=EC;
    save(strcat('/glade/scratch/samantha/EOF_PC_iLMEmem',num2str(rr),'-',useint,'.mat'),'PCts','V','EOFmat')

    % Create bootstrap time series
    bootts=zeros([nboot,size(d18obin)]);
    nblk=floor(size(d18obin,1)/ns);    % number of possible blocks
    
    % Get set of mean values for sorting on
    blkmnts=zeros(nblk,length(proxylats));
    for bb=1:nblk
        blkmnts(bb,:)=nanmean(d18obin((bb-1)*ns+1:bb*ns,:),1);
    end
    
    for bb=1:nboot
        % Pick random starting block
        bind=randi(nblk,[1 length(proxylats)]);
        
        % Assign that random starting block as the first set of 
        % values in the bootstrap time series 
        for pp=1:length(proxylats)
            bootts(bb,1:ns,pp)=d18obin((bind(pp)-1)*ns+1:bind(pp)*ns,pp);
        
            % Loop through to generate the rest of the time series
            for ss=2:nblk
                bootmn=squeeze(nanmean(bootts(bb,(ss-2)*ns+1:(ss-1)*ns,:),2));

                % Sort blocks by similarity of mean states
                [tmp,ind]=sort(abs(blkmnts(:,pp)-bootmn(pp)));
                
                % choose a random decently close block
%                 btmp=randi([2,10],1);
                btmp=randi([2,5],1);
                bootts(bb,(ss-1)*ns+1:ss*ns,pp)=d18obin((ind(btmp)-1)*ns+1:ind(btmp)*ns,pp); 
            end
        end
        
        % Do SVD analysis on the bootstrap-resampled data
        [V,EOFtmp,ECboot,error]=EOF(squeeze(bootts(bb,:,:)),15);
%         V/sum(V)        % percent variance in each mode
        ECtsboot(bb,rr,:)=ECboot(:,pcnum);
        EOFsboot(bb,rr,:)=EOFtmp(:,pcnum);     
        Vboot(bb,rr,:)=V/sum(V);
    end
    
%     % Do gridpoint regression
%     for la=1:size(slpbin,2)
%         for lo=1:size(slpbin,3)
%             [tmp,tmpint]=regress(slpbin(:,la,lo),cat(2,ones(length(binyr),1),EC(:,pcnum)));
%             regslope(rr,la,lo)=tmp(2);
%         end
%     end
    % Do gridpoint correlation
    for la=1:size(slpbin,2)
        for lo=1:size(slpbin,3)
            regslope(rr,la,lo)=corr(slpbin(:,la,lo),EC(:,pcnum));
            dslpregslope(rr,la,lo)=corr(slpbin(:,la,lo),dslpbin);
        end
    end
end

% Plot ensemble member results
for rr=1:length(runnames)
    % SLP regression slopes
    figure(1)
    clf
    m_proj('robinson','lon',[0 360],'lat',[-90 90])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lon,lat,squeeze(regslope(rr,:,:)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    [C,hc]=m_contour(lon,lat,squeeze(dslpregslope(rr,:,:)),[0 0.2 0.4 0.6 0.8 1],'Color','k','LineWidth',2.1);
    clabel(C,hc,'Color','k','FontSize',16)
    [C,hc]=m_contour(lon,lat,squeeze(dslpregslope(rr,:,:)),[-1 -0.8 -0.6 -0.4 -0.2],'--','Color','k','LineWidth',2.1);
    clabel(C,hc,'Color','k','FontSize',16)

    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title(strcat('SLP correlated with precip \delta^{18}O PC',num2str(pcnum)),'FontSize',20)
    colorbar  
    if dodtr == 1
        strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig'),'fig')
    else
        strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig'),'fig')
    end
    
    % PC time series
    figure(1)
    clf
    plot(binyr,ECts(rr,:),'Color','r','LineWidth',2)
    set(gca,'FontSize',24)
    hold all
    boottmp=prctile(ECtsboot,90,1);
    plot(binyr,squeeze(boottmp(:,rr,:)),'--','Color','k','LineWidth',1)
    title(strcat('\delta^{18}O PC',num2str(pcnum),' time series'))
    if dodtr == 1
        strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig'),'fig')
    else
        strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig'),'fig')
    end

    % PC spectra
    figure(1)
    clf
    nfft=length(ECts(rr,:));
    h = spectrum.welch('Hann',wlen,100*noverlap/wlen);
    hpsd = psd(h,ECts(rr,:),'NFFT',nfft,'Fs',Fs,'ConfLevel',0.9);
    Pw = hpsd.Data; 
    Fw = hpsd.Frequencies;
    CI = hpsd.ConfInterval;

    h1=loglog(Fw,Pw,'Color','k','LineWidth',2);
    hold all
    loglog(Fw,CI,'--','Color','k','LineWidth',2);
    set(gca,'FontSize',24,'box','on','LineWidth',2)
    xlabel('Freq. (cycles/month)')
    ylabel('Power')
    title(strcat('\delta^{18}O PC',num2str(pcnum),' spectrum'))
    if dodtr == 1
        strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'spec_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'spec_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr.fig'),'fig')
    else
        strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'spec_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'spec_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'.fig'),'fig')
    end    
    
    % Loadings
    figure(1)
    clf
    m_proj('robinson','lon',[0 360],'lat',[-90 90])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    shading flat
    hold all
    [x,y] = m_ll2xy(proxylons,proxylats);
    for pp=1:length(proxylats)
        scatter(x(pp),y(pp),320,EOFs(rr,pp),'o','filled','MarkerEdgeColor','k','LineWidth',3)
    end
    m_coast('color',[0 0 0]);
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title(strcat('Soil \delta^{18}O EOF',num2str(pcnum),' loadings'),'FontSize',20)
    colorbar  
    if dodtr == 1
        strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_dtr_codewdslp.fig'),'fig')
    else
        strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig')
        saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_mem',num2str(rr),'_codewdslp.fig'),'fig')
    end

end

% Plot ensemble mean results
% Regression plus loadings
figure(1)
clf
m_proj('robinson','lon',[0 360],'lat',[-90 90])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
m_pcolor(lon,lat,squeeze(nanmean(regslope,1)));
shading flat
hold all
[x,y] = m_ll2xy(proxylons,proxylats);
for pp=1:length(proxylats)
    scatter(x(pp),y(pp),320,squeeze(nanmean(EOFs(:,pp),1)),'o','filled','MarkerEdgeColor','k','LineWidth',3)
end
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
[C,hc]=m_contour(lon,lat,squeeze(nanmean(dslpregslope,1)),[0 0.2 0.4 0.6 0.8 1],'Color','k','LineWidth',2.1);
%clabel(C,hc,'Color','k','FontSize',16)
[C,hc]=m_contour(lon,lat,squeeze(nanmean(dslpregslope,1)),[-1 -0.8 -0.6 -0.4 -0.2],'--','Color','k','LineWidth',2.1);
%clabel(C,hc,'Color','k','FontSize',16)
title(strcat('SLP correlated with soil \delta^{18}O PC',num2str(pcnum)),'FontSize',20)
colorbar  
if dodtr == 1
    strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_',useint,'sites_withEOFloadings_wdslpreg_dtr.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_',useint,'sites_withEOFloadings_wdslpreg_dtr.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_',useint,'sites_withEOFloadings_wdslpreg.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_',useint,'sites_withEOFloadings_wdslpreg.fig'),'fig')
end

% PC time series
figure(1)
clf
plot(binyr,ECts)
set(gca,'FontSize',24)
title(strcat('\delta^{18}O PC',num2str(pcnum),' time series'))
legend({'1','2','3'})
if dodtr == 1
    strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme',useint,'sites_dtr.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme',useint,'sites_dtr.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme',useint,'sites.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme',useint,'sites.fig'),'fig')
end

% Loadings
figure(1)
clf
m_proj('robinson','lon',[0 360],'lat',[-90 90])
ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
shading flat
hold all
[x,y] = m_ll2xy(proxylons,proxylats);
for pp=1:length(proxylats)
    scatter(x(pp),y(pp),320,squeeze(nanmean(EOFs(:,pp),1)),'o','filled','MarkerEdgeColor','k','LineWidth',3)
end
m_coast('color',[0 0 0]);
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
title(strcat('Soil \delta^{18}O EOF',num2str(pcnum),' loadings'),'FontSize',20)
colorbar  
if dodtr == 1
    strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_dtr.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_dtr.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme.fig'),'fig')
end

% Significance of modes
sigl=prctile(squeeze(Vboot(1,:,:)),90,1);
figure(1)
clf
plot(1:length(Vs),Vs);
hold all
plot(1:length(Vs),sigl,'--','Color','k')
if dodtr == 1
    strcat('/glade/scratch/samantha/plots/d18OsoilEOFvar_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_dtr.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOFvar_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme_dtr.fig'),'fig')
else
    strcat('/glade/scratch/samantha/plots/d18OsoilEOFvar_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOFvar_0-',num2str(levthr),'cm_',useint,'sites_',num2str(yrint(1)),'-',num2str(yrint(2)),'_ilme.fig'),'fig')
end