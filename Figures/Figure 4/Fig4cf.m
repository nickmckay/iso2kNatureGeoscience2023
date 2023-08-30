% Code to make a PCA figure following the Iso2k conventions, look at
% associated SLP signatures
% April 2020
% Sam Stevenson

runnames={'b.ie12.B1850CN.f19_g16.001','b.ie12.B1850C5CN.f19_g16.LME.002','b.ie12.B1850C5CN.f19_g16.LME.003'};
% runnames={'b.ie12.B1850CN.f19_g16.001'};
cm=1-[[1*(1:-.1:0)' 1*(1:-.1:0)' 0*(1:-.1:0)'];[0*(0:.1:1)' 1*(0:.1:1)' 1*(0:.1:1)']];
startyr=850;
refper=[1961,1990];
windlen=30;
use20th=0;
levthr=0.1;
avgrege=[-5,5,210,270];	% Averaging regions, following Tokinaga et al. (2012)
avgregw=[-5,5,90,150];
pcnum=1;

% Parameters for computing soil water d18O
Rstd=2005.2e-6; % O18/O16 ratio in standard (here, VSMOW)
Mh2o=18.015;    % mean molecular mass of water
Mh218o=20.013;  % molecular mass of H218O
Mh216o=18.009;  % molecular mass of H216O

% Coordinates
unc=netcdf(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runnames{1},'.clm2.h0.H2OSOI.085001-184912.nc'));
lat=unc{'lat'}(:);
lon=unc{'lon'}(:);
lev=unc{'levgrnd'}(:);
mylev=find(lev <= levthr);     % 0-1cm depth

if use20th == 0
    binlen=30;  % bin every xx years
    
    % Locations of Iso2k proxies
    % "effectiveMoisture_PCA_850-1840_30yrbins"
    % NOTE: USE 30 YEAR BINS WHEN COMPARING WITH THIS DATA
    proxylats=[54.685,50.83,19.8333,40.1,16.98,37.87,20.61,38.34,38.34,48.2,46,60.1,11.35,45,47.0457, ...
        37,48.5,39.7,48.8,11.9,38.8667,19.9763,41.0658,37.87,20.63,-3.53,55.5,1.4,-5.2,40.5, ...
        61.4817,31.15,13.9];
    proxylons=[-122.617,-116.39,-88.75,-119.6,-89.666,-119.16,-89.715,34.46,34.46,-114.4,-94.7, ...
        -133.8,39.72,-110.5999,-113.1426,100,-119.6,-107.4,-118.2,-85.9167,93.95,76.5077, ...
        20.6728,-119.16,-87.61,119.2,-13.9,119.08,117.48,4.03,-19.536,97.033,-16.6];
    proxytype=["LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment", ...
        "LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment", ...
        "LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment", ...
        "LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment","LakeSediment", ...
        "LakeSediment","MarineSediment","MarineSediment","MarineSediment","MarineSediment", ...
        "MarineSediment","MarineSediment","Wood","MolluskShells"];
else
    binlen=10;
    
    % Locations of Iso2k proxies
    % NOTE: USE 10 YEAR BINS WHEN COMPARING WITH THIS DATA
    proxylats=[];
    proxylons=[];
    proxytype=[];
end
proxylons(proxylons < 0)=proxylons(proxylons < 0)+360;
syms=proxytype;
syms(proxytype == "GlacierIce")='d';
syms(proxytype == "LakeSediment")='s';
syms(proxytype == "Speleothem")='d';
syms(proxytype == "Wood")='^';

figure(1)
clf
for rr=1:length(runnames)
    runname=runnames{rr}
    
    % Read in SLP
    if use20th == 0
        atmfile=strcat('/glade/scratch/samantha/iCESM/atm/PSL/',runname,'.cam.h0.PSL.085001-184912.nc');
    else
        atmfile=strcat('/glade/scratch/samantha/iCESM/atm/PSL/',runname,'.cam.h0.PSL.185001-200512.nc');
    end
    slp=nc_varget(atmfile,'PSL')/100.;
    if rr == 1
        if use20th == 0
            atmtime=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.085001-184912.nc'),'time');
            [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[850 1 1]);
        else
            atmtime=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.185001-200512.nc'),'time');
            [pyr,pmon,pdy]=datenumnoleap(atmtime-29,[1850 1 1]);
        end
        uyr=unique(pyr);
        binyr=min(uyr):binlen:(max(uyr)-binlen+1);
        pbin=zeros(length(binyr),length(proxylats));
        d18obin=zeros(length(binyr),length(proxylats));
        slpbin=zeros(length(binyr),size(slp,2),size(slp,3));
        regslope=zeros(length(runnames),size(slp,2),size(slp,3));
        ECts=zeros(length(runnames),length(binyr));
        EOFs=zeros(length(runnames),length(proxylats));
        Vs=zeros(length(runnames),25);
        dslpregslope=zeros(length(runnames),size(slp,2),size(slp,3));
        dslpbin=zeros(length(binyr),1);
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
        if pp == 29 
            mylat=find(lat >= proxylats(pp)-3 & lat <= proxylats(pp)+3);
            mylon=find(lon >= proxylons(pp)-3 & lon <= proxylons(pp)+3);
        end

        % Read in variables, calculate soil d18O    
        if use20th == 0
            h2osoi=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.085001-184912.nc'),'H2OSOI',[0 min(mylev)-1 min(mylat)-1 min(mylon)-1],[-1 length(mylev) length(mylat) length(mylon)]);
            h2osoi_h218o=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI_H218O/',runname,'.clm2.h0.H2OSOI_H218O.085001-184912.nc'),'H2OSOI_H218O',[0 min(mylev)-1 min(mylat)-1 min(mylon)-1],[-1 length(mylev) length(mylat) length(mylon)]);
        else
            h2osoi=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI/',runname,'.clm2.h0.H2OSOI.185001-200512.nc'),'H2OSOI',[0 min(mylev)-1 min(mylat)-1 min(mylon)-1],[-1 length(mylev) length(mylat) length(mylon)]);
            h2osoi_h218o=nc_varget(strcat('/glade/scratch/samantha/iCESM/lnd/H2OSOI_H218O/',runname,'.clm2.h0.H2OSOI_H218O.185001-200512.nc'),'H2OSOI_H218O',[0 min(mylev)-1 min(mylat)-1 min(mylon)-1],[-1 length(mylev) length(mylat) length(mylon)]);
        end

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
    end
    
    % SVD analysis
    d18obin(abs(d18obin) > 1e10)=0;
    d18obin(isnan(d18obin))=0;
    [V,EOFtmp,EC,error]=EOF(d18obin,25);
%     [U,S,V]=svd(d18O);
    V/sum(V)        % percent variance in each mode
    ECts(rr,:)=EC(:,pcnum);
    EOFs(rr,:)=EOFtmp(:,pcnum);
    Vs(rr,:)=V/sum(V);
    

    
%     % Do gridpoint regression
%     for la=1:size(slpbin,2)
%         for lo=1:size(slpbin,3)
%             [tmp,tmpint]=regress(slpbin(:,la,lo),cat(2,ones(length(binyr),1),EC(:,pcnum)));
%             regslope(rr,la,lo)=tmp(2);
%         end
%     end

    % Make sure signs match based on EOF loadings
    %if EOFtmp(32,pcnum) > 0    %PC1
    if EOFtmp(28,pcnum) > 0
        EOFtmp(:,pcnum)=-1*EOFtmp(:,pcnum);
        EC(:,pcnum)=-1*EC(:,pcnum);
    end

    
    savets=cat(2,binyr',EC);
    save(strcat('/glade/scratch/samantha/PCtimeseries_iLMEmem',num2str(rr),'_850-1850.txt'),'savets','-ASCII')

    savets=cat(2,binyr',dslpbin);
    save(strcat('/glade/scratch/samantha/dSLPtimeseries_iLMEmem',num2str(rr),'_850-1850.txt'),'savets','-ASCII')
    
    EOFmat=EOFtmp;
    PCts=EC;
    save(strcat('/glade/scratch/samantha/EOF_PC_iLMEmem',num2str(rr),'_850-1850.mat'),'PCts','V','EOFmat')
    

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
    % SLp regression slopes
    figure(1)
    clf
    m_proj('robinson','lon',[0 360],'lat',[-90 90])
    ylabel('Lat (^{\circ}N)','FontSize',18,'FontWeight','bold')
    m_pcolor(lon,lat,squeeze(regslope(rr,:,:)));
    shading flat
    hold all
    m_coast('color',[0 0 0]);
    [C,hc]=m_contour(lon,lat,squeeze(dslpregslope(rr,:,:)),[0 0.2 0.4 0.6 0.8 1],'Color','k','LineWidth',2.1);
    clabel(C,hc,'Color','k','FontSize',16)
    [C,hc]=m_contour(lon,lat,squeeze(dslpregslope(rr,:,:)),[-1 -0.8 -0.6 -0.4 -0.2],'--','Color','k','LineWidth',2.1);
%     clabel(C,hc,'Color','k','FontSize',16)
    m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
    set(gca,'DataAspectRatio',[1 0.75 1],'FontSize',24)
    title(strcat('SLP correlated with precip \delta^{18}O PC',num2str(pcnum)),'FontSize',20)
    colorbar  
    strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_850-1850_ilme_mem',num2str(rr),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_850-1850_ilme_mem',num2str(rr),'.fig'),'fig')
    
    % PC time series
    figure(1)
    clf
    plot(binyr,ECts(rr,:))
    set(gca,'FontSize',24)
    title('\delta^{18}O PC1 time series')
    strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_850-1850_ilme_mem',num2str(rr),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_850-1850_ilme_mem',num2str(rr),'.fig'),'fig')

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
    strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_effectivemoisturesites_850-1850_ilme_mem',num2str(rr),'.fig')
    saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_effectivemoisturesites_850-1850_ilme_mem',num2str(rr),'.fig'),'fig')


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
[C,hc]=m_contour(lon,lat,squeeze(nanmean(dslpregslope,1)),[0 0.2 0.4 0.6 0.8 1],'Color','k','LineWidth',2.1);
%clabel(C,hc,'Color','k','FontSize',16)
[C,hc]=m_contour(lon,lat,squeeze(nanmean(dslpregslope,1)),[-1 -0.8 -0.6 -0.4 -0.2],'--','Color','k','LineWidth',2.1);
% clabel(C,hc,'Color','k','FontSize',16)
m_grid('xaxis','bottom','tickdir','out','linewidth',3,'FontSize',24);
title(strcat('SLP correlated with soil \delta^{18}O PC',num2str(pcnum)),'FontSize',20)
colorbar  
strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_850-1850_ilme_withEOFloadings_wdslpcorr_ensmn_signfix.fig')
saveas(gcf,strcat('/glade/scratch/samantha/plots/slpcorr_d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_850-1850_ilme_withEOFloadings_wdslpcorr_ensmn_signfix.fig'),'fig')


% PC time series
figure(1)
clf
plot(binyr,ECts)
set(gca,'FontSize',24)
title('\delta^{18}O PC1 time series')
legend({'1','2','3'})
strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_850-1850_ilme.fig')
saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilPC',num2str(pcnum),'_0-',num2str(levthr),'cm_850-1850_ilme.fig'),'fig')

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
strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_effectivemoisturesites_850-1850_ilme_ensmn_signfix.fig')
saveas(gcf,strcat('/glade/scratch/samantha/plots/d18OsoilEOF',num2str(pcnum),'_0-',num2str(levthr),'cm_effectivemoisturesites_850-1850_ilme_ensmn_signfix.fig'),'fig')

