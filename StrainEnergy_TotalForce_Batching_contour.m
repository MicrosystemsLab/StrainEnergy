%Aleksandra Denisin (adenisin@stanford.edu), Pruitt Laboratory
%3/11/2017
%Code written for analysis of PIV-FTTC Traction Force
%We have PIV info and TFM info from Tseng script

function [OutputMatrix]=StrainEnergy_TotalForce_Batching_contour
close all
clear all

conversion =0.10317; %um/pixel

%read traction force output file from ImageJ plugins from Tseng, ask user to enter the filename...
searchitemPIV = input('Type the base of your PIV files (i.e., *30kPa_PIV*): ');
PIVfiles = dir(searchitemPIV);
searchitemTFM = input('Type the base of your TFMfiles (i.e., *30kPa_TFM*): ');
TFMfiles = dir(searchitemTFM);

filenameROI = input('Enter the ROI file name:');

ROI=ReadImageJROI(filenameROI);
%%
%need to import in ROI of cell
%ROI importer: http://dylan-muir.com/articles/read_ij_roi/

%initiate variables w zeros
enddim = zeros(size(PIVfiles,1),1);
cellarea = zeros(size(PIVfiles,1),1); %m^2 cell area
totaltractionstress = zeros(size(PIVfiles,1),1);
averagetractionstress = zeros(size(PIVfiles,1),1);
stdevtractionstress =zeros(size(PIVfiles,1),1);
totalcellforce =zeros(size(PIVfiles,1),1); %N
totalstrainenergy = zeros(size(PIVfiles,1),1); %J
OutputMatrix = [];
for z = 1:size(PIVfiles,1)
    %% Load PIV data into matrices
    DataName = PIVfiles(z).name;
    DataNameShort = DataName(1:end-4);
    PIVdata = load(PIVfiles(z).name);
    x = PIVdata(:,1);
    y = PIVdata(:,2);
    Px = PIVdata(:,3);
    Py = PIVdata(:,4);
    Pmag = PIVdata(:,5); %norm of the displacement
    %magcheck = sqrt( F(:,1).^2 + F(:,2).^2 );
    
    enddim(z) = PIVdata(end, 1)+40; %end dimensions of our figure, should be around 816
    mask = poly2mask(ROI{1,z}.mnCoordinates(:,1),ROI{1,z}.mnCoordinates(:,2), enddim(z), enddim(z));
    
    %For u and v, need to make a matrix
    spacingP = (x(2,1)-x(1,1));
    udisp = zeros((max(x)-x(1,1))/spacingP+1,(max(y)-y(1,1))/spacingP+1);
    vdisp = zeros((max(x)-x(1,1))/spacingP+1,(max(y)-y(1,1))/spacingP+1);
    
    k = 1;
    %this is to take our spacing datapoints from columns into spatial distributions
    for j = 1:size(udisp,2)
        for i = 1:size(udisp,1)
            udisp(i,j) = x(k);
            vdisp(i,j) = y(k);
            k=k+1;
        end
    end
    Pxx = zeros((max(x)-x(1,1))/spacingP+1,(max(y)-y(1,1))/spacingP+1);
    Pyy = zeros((max(x)-x(1,1))/spacingP+1,(max(y)-y(1,1))/spacingP+1);
    
    m = 1;
    %this is to take our TFM datapoints from columns into spatial distributions
    for j = 1:size(Pxx,2)
        for i = 1:size(Pxx,1)
            Pxx(i,j) = Px(m);
            Pyy(i,j) = Py(m);
            m=m+1;
        end
    end
    DispMagMatrix= sqrt( Pxx.^2 + Pyy.^2 );
    DispMagMatrixTrans = DispMagMatrix';
    %need to convert to meters (conversion is um/pixel)
    DispMagMatrixTransinM = DispMagMatrixTrans.*conversion*1E-6;
    %imagesc(DispMagMatrixTrans);
    %%
    %%Load TFM data
    TFMdata = load(TFMfiles(z).name);
    X = TFMdata(:,1);
    Y = TFMdata(:,2);
    Sx = TFMdata(:,3);
    Sy = TFMdata(:,4);
    mag = TFMdata(:,5); %norm of the tractions
    %magcheck = sqrt( F(:,1).^2 + F(:,2).^2 );
    
    %For u and v, need to make a matrix
    spacingS = (X(2,1)-X(1,1));
    u = zeros((max(X)-X(1,1))/spacingS+1,(max(Y)-Y(1,1))/spacingS+1);
    v = zeros((max(X)-X(1,1))/spacingS+1,(max(Y)-Y(1,1))/spacingS+1);
    
    m = 1;
    %this is to take our spacing datapoints from columns into spatial distributions
    for j = 1:size(u,2)
        for i = 1:size(u,1)
            u(i,j) = X(m);
            v(i,j) = Y(m);
            m=m+1;
        end
    end
    
    Sxx = zeros((max(X)-X(1,1))/spacingS+1,(max(Y)-Y(1,1))/spacingS+1);
    Syy = zeros((max(X)-X(1,1))/spacingS+1,(max(Y)-Y(1,1))/spacingS+1);
    
    n = 1;
    %this is to take our TFM datapoints from columns into spatial distributions
    for j = 1:size(Sxx,2)
        for i = 1:size(Sxx,1)
            Sxx(i,j) = Sx(n);
            Syy(i,j) = Sy(n);
            n=n+1;
        end
    end
    Smagmatrix= sqrt( Sxx.^2 + Syy.^2 );
    %need to transform matrix to match how ROIs are taken from ImageJ
    SmagmatrixTrans = Smagmatrix';
    
%     figure(1)
%     figure1 = imagesc(SmagmatrixTrans);
%    colorbar
%    figureHandle = gcf;    %# make all text in the figure to size 14
%     set(findall(figureHandle,'type','text'),'fontSize',14)
%     
    %%
    %calculate the strain energy density (units = J/m^2)
    %w(r) = 1/2*tractionstress(r)*displacement(r)
    %tractionstress[Pa]*displacement[m] = [N/m^2]*m = N/m = J/m^2
    % because J = N*m
    
    w = 0.5*SmagmatrixTrans.*DispMagMatrixTransinM;
    
    %plot map of strain energy density (J/m^2)
    %figure(2)
    %figure2 = imagesc(w);
    %colorbar
    %%
    %Interpolate both traction stress and strain energy
    %our PIV analysis had a window size of 16
    %so we multiply our matrix results to expand all pixels in the 'original'
    %image which is generally size of the cropped image (800 pix) plus a border
    %of 16 (816)
    
    %Interpolate traction data for every pixel dimensions
    [MeshX,MeshY] = meshgrid(min(x):spacingS:max(max(x)));
    [XI,YI] = meshgrid(1:1:enddim(z));
    StressInterpolated= interp2(MeshX,MeshY,SmagmatrixTrans,XI,YI);
    
    %Interpolate strain energy calculation for every pixel dimensions
    StrainEnergyInterpolated= interp2(MeshX,MeshY,w,XI,YI);
    
    %we want to apply the ROI
    StressInterpolatedROI = StressInterpolated.*double(mask);
    StressInterpolatedROI(StressInterpolatedROI==0)=NaN;
    StrainEnergyInterpolatedROI = StrainEnergyInterpolated.*double(mask);
    StrainEnergyInterpolatedROI(StrainEnergyInterpolatedROI==0)=NaN;
    
    figure(3)
    figure3 = imagesc(StressInterpolatedROI);
    title('Stress in Pa')
    colormap(jet)
    colorbar
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14) % make all text in the figure to size 14
    print(strcat('Stress_',DataNameShort),'-dpng')
    
    figure(4)
    figure4 = imagesc(StrainEnergyInterpolatedROI);
    title('Strain Energy Density in J/m^2')
    colormap(jet)
    colorbar
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14) % make all text in the figure to size 14
    print(strcat('StrainE_',DataNameShort),'-dpng')
    
      %Try doing 'linescans' through different areas of the cell
    SE = strel('square',50); %structuring element for eroding
    maskerode = imerode(mask', SE);
    
    figure(5)
    figure5 = imagesc(StrainEnergyInterpolatedROI);
    title('Strain Energy Density in J/m^2')
    colormap(jet)
    colorbar
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14) % make all text in the figure to size 14
    %print(strcat('StrainE_',DataNameShort),'-dpng')
    hold on
    B = bwboundaries(maskerode);
     plot(B{1,1}(:,1),B{1,1}(:,2),'r','LineWidth',5)
    A = improfile(StrainEnergyInterpolatedROI,B{1,1}(:,1),B{1,1}(:,2));
   
    figure(6)
    figure6 = plot(A);
    title('Strain Energy Profile in J/m^2')
    
    SE2 = strel('square',75); %structuring element for eroding
    maskerode2 = imerode(mask', SE2);
    
    figure(7)
    figure7 = imagesc(StrainEnergyInterpolatedROI);
    title('Strain Energy Density in J/m^2')
    colormap(jet)
    colorbar
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14) % make all text in the figure to size 14
    %print(strcat('StrainE_',DataNameShort),'-dpng')
    hold on
    B2 = bwboundaries(maskerode2);
    plot(B2{1,1}(:,1),B2{1,1}(:,2),'r','LineWidth',2)
     hold on
     plot(B2{1,1}(1:10,1),B2{1,1}(1:10,2),'w', 'LineWidth',2)%,'LineWidth',5)
    A2 = improfile(StrainEnergyInterpolatedROI,B2{1,1}(:,1),B2{1,1}(:,2));
   
    figure(8)
    figure8 = plot(A2);
    title('Strain Energy Profile in J/m^2')
    
    
    %export measurements
    %Total Cell-Generated Force = sum(traction stress in the cell area) * cell area
    %units: Pa*(m^2) = N/m^2*m^2 =N
    cellarea(z) = nansum(nansum(mask))*conversion*conversion*1E-6*1E-6; %m^2 cell area
    totaltractionstress(z) = nansum(StressInterpolatedROI(:)); %N/m^2
    averagetractionstress(z) = nanmean(StressInterpolatedROI(:)); %N/m^2
    stdevtractionstress(z)=nanstd(StressInterpolatedROI(:));%N/m^2
    totalcellforce (z)=totaltractionstress(z)*cellarea(z); %N 
    totalstrainenergy(z) = nansum(StrainEnergyInterpolatedROI(:))*cellarea(z);  %J
end
OutputMatrix = [];
OutputMatrix = horzcat(cellarea,averagetractionstress,stdevtractionstress,totaltractionstress, totalcellforce,totalstrainenergy);

disp('finished');
end


