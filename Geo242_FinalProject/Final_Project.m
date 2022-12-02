clear all; close all; clc

load('Auto.mat'); %load('Depth_Picks1.mat');
%Auto.mat contains the depth values of the Lithosphere-Asthenosphere
%Boundary from an autopicking function

load('CCP_Sp_20.6251_California_Sp_2_100s_Crust_1.0_SL14AK135_SL14AK135_23.09.22_300km_deg_1.mat')
%the CCP file is needed to mask the parts of the region that do not
%have sufficient data coverage at 70km depth and to obtain
%latitude/longitude information 

States = shaperead('usastatehi', 'UseGeoCoords', true); 
%information about state boundaries so that they can be plotted on the map 

lon = model_lons; lat = model_lats'; %latitude/longitudes from CCPs
lontmp = tmplon(:,1); lattmp = tmplat(1,:)'; 
%latitude/longitude vectors corresponding to original matrix dimensions
lon_fine = [lontmp(1):0.05:lontmp(end)]; %new longitude vector for interpolation
lat_fine = [lattmp(1):0.05:lattmp(end)]'; %new latitude vector for interpolation 

events_70km=weighted_n_events(:,141); %pulls weighted events at depth 70km, comes from CCP file 
events_fine= (griddata(lon,lat,events_70km,tmplon(:,1),tmplat(1,:)))'; 
%takes the vector with information about the number of events at each point
%and turns it into a grid that is the same size as the matrix with the
%autopicked boundary depths 

%use the newly defined grid of weighted events to remove values when
%log10(# of events) <= 0.25
for m = 1:length(tmplat(1,:))
    for n = 1:length(tmplon(:,1))
        if log10(events_fine(n,m))<= 0.25;  LABdp(n,m) = NaN; end
        %values removed by replacing the actual value with a NaN
    end
end
labdepth1 = double(LABdp); %makes sure the matrix of boundary depths is a double, not a single, 
%as is required by the interpolation methods
myColorMap = fliplr(jet); %set the color map option, flip it so that red is shallow values, blue is deep


%% Three Dimensional Interpolations
%Just to see what the boundary picks look like as a 3D surface

%Attempt to make it a surface - interpn & surf functions
S = interpn(lontmp,lattmp,labdepth1,lon_fine,lat_fine);
%double the size of the grid, interpolate using linear method (default) for
%gridded data in ngrid format
figure
%plot the data as a surface
surf(lon_fine,lat_fine,-S.','EdgeColor','interp'); colormap jet; colorbar; hold on
title('Lithosphere-Asthenosphere Boundary Depth')
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end %plot the states
ylim([30 44]); xlim([-125 -111]); 


% Gridded Interpolation & surf functions 
GI = griddedInterpolant({lontmp,lattmp},labdepth1); %different interpolation function with linear method as default
grinterp = GI({lon_fine,lat_fine});
figure
surf(lon_fine,lat_fine,-grinterp.','EdgeColor','interp'); colormap jet; colorbar; hold on
title('Lithosphere-Asthenosphere Boundary Depth')
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-125 -111]); 

%% interp2 Interpolation
%First use the interp2 function to do interpolation, use three different
%interpolation methods: linear, nearest neighbor, cubic 

% Interp2 - linear
I2lin = interp2(lattmp,lontmp,labdepth1,lat_fine,lon_fine);
figure; %subplot(2,2,1)
hlin = pcolor(lon_fine,lat_fine,I2lin.'); colormap(myColorMap); colorbar; hold on 
set(hlin,'EdgeColor','none')
title('Linear Interpolation')
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]); 

% Interp2 - Nearest Neighbor
I2nn = interp2(lattmp,lontmp,labdepth1,lat_fine,lon_fine,'nearest');
figure %subplot(2,2,2)
hnn = pcolor(lon_fine,lat_fine,I2nn.'); colormap(myColorMap); colorbar; hold on 
set(hnn,'EdgeColor','none');title('Nearest Neighbor Interpolation')
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]);

% Interp2 - Cubic
I2c = interp2(lattmp,lontmp,labdepth1,lat_fine,lon_fine,'cubic');
figure %subplot(2,2,3)
hc = pcolor(lon_fine,lat_fine,I2c.'); colormap(myColorMap); colorbar; hold on 
set(hc,'EdgeColor','none'); title('Cubic Interpolation')
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]);

%% imresize Interpolation
%Now use the imresize function to do interpolation, use three different
%interpolation methods: nearest neighbor, bilinear, bicubic 
nnresize = imresize(labdepth1,[287 279],'nearest');
figure
hnnr = pcolor(lon_fine,lat_fine,nnresize.');hold on
set(hnnr,'EdgeColor','none'); colormap(myColorMap); colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]); title('Nearest Neighbor Interpolation');

blresize = imresize(labdepth1,[287 279],'bilinear');
figure
hbl = pcolor(lon_fine,lat_fine,blresize.');hold on
set(hbl,'EdgeColor','none'); colormap(myColorMap); colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]); title('Bilinear Interpolation');

bcresize = imresize(labdepth1,[287 279],'bicubic');
figure
hbc = pcolor(lon_fine,lat_fine,bcresize.');hold on
set(hbc,'EdgeColor','none'); colormap(myColorMap); colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 43.9]); xlim([-126 -111.2]); title('Bicubic Interpolation');

%% Look at differences in 2D interpolations

%linear
I2lin_avg = nanmean(nanmean(I2lin));
%I2 nearest neighbor
I2nn_avg = nanmean(nanmean(I2nn));
%cubic
I2c_avg = nanmean(nanmean(I2c));
%imresize nearest neighbor
nnresize_avg = nanmean(nanmean(nnresize));
%bilinear
blresize_avg = nanmean(nanmean(blresize));
%bicubic
bcresize_avg = nanmean(nanmean(bcresize));

%comparing differences in models - start w/ linear
for j=1:length(lon_fine)
    for k=1:length(lat_fine)
        linear.I2nn(j,k) = I2lin(j,k) - I2nn(j,k);
        linear.I2c(j,k) = I2lin(j,k) - I2c(j,k);
        linear.nnresize(j,k) = I2lin(j,k) - nnresize(j,k);
        linear.blresize(j,k) = I2lin(j,k) - blresize(j,k);
        linear.bcresize(j,k) = I2lin(j,k) - bcresize(j,k);
    end 
end 
%Calculation is taking simple difference of each of the models from the
%linear interpolation. If a value is positive, the linear interpolation had
%a deeper depth value for the bounday; if a value is negative, the other
%interpolation method has a deeper value for the boundary depth 

figure
subplot(2,3,1)
l1=pcolor(lon_fine,lat_fine,linear.I2nn.'); hold on
set(l1,'EdgeColor','none');  colorbar;  
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Interp2 Linear - Interp2 NN');
subplot(2,3,2)
l2=pcolor(lon_fine,lat_fine,linear.I2c.'); hold on 
set(l2,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Interp2 Linear - Interp2 Cubic');
subplot(2,3,3)
l3=pcolor(lon_fine,lat_fine,linear.nnresize.'); hold on 
set(l3,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Interp2 Linear - Imresize NN');
subplot(2,3,4)
l4=pcolor(lon_fine,lat_fine,linear.blresize.'); hold on 
set(l4,'EdgeColor','none');  colorbar;
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Interp2 Linear - Imresize Bilinear');
subplot(2,3,6)
l5=pcolor(lon_fine,lat_fine,linear.bcresize.'); hold on 
set(l5,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Interp2 Linear - Imresize Bicubic');

%comparison - Imresize nearest neighbor 
for j=1:length(lon_fine)
    for k=1:length(lat_fine)
        nnr.I2nn(j,k) = nnresize(j,k) - I2nn(j,k);
        nnr.I2c(j,k) = nnresize(j,k) - I2c(j,k);
        nnr.I2lin(j,k) = nnresize(j,k) - I2lin(j,k);
        nnr.blresize(j,k) = nnresize(j,k) - blresize(j,k);
        nnr.bcresize(j,k) = nnresize(j,k) - bcresize(j,k);
    end 
end 
%Calculation is taking simple difference of each of the models from the
%imsize nearest neighbor interpolation. If a value is positive, the imresize 
%nearest neighbor interpolation had a deeper depth value for the bounday; 
%if a value is negative, the other interpolation method has a deeper value 
%for the boundary depth 

figure
subplot(2,3,1)
l1=pcolor(lon_fine,lat_fine,nnr.I2nn.'); hold on
set(l1,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Imresize NN - Interp2 NN');
subplot(2,3,2)
l2=pcolor(lon_fine,lat_fine,nnr.I2c.'); hold on 
set(l2,'EdgeColor','none');  colorbar;  
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Imresize NN - Interp2 Cubic');
subplot(2,3,3)
l3=pcolor(lon_fine,lat_fine,nnr.I2lin.'); hold on 
set(l3,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Imresize NN - Interp2 Linear');
subplot(2,3,4)
l4=pcolor(lon_fine,lat_fine,nnr.blresize.'); hold on 
set(l4,'EdgeColor','none');  colorbar; 
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]); xlim([-126 -111.2]); title('Imresize NN - Imresize Bilinear');
subplot(2,3,6)
l5=pcolor(lon_fine,lat_fine,nnr.bcresize.'); hold on 
set(l5,'EdgeColor','none');  colorbar;  
for i = 1:length(States); plot(States(i).Lon, States(i).Lat, 'k-'); end
ylim([30 44]);xlim([-126 -111.2]); title('Imresize NN - Imresize Bicubic');