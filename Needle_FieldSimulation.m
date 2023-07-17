%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code to simulate field distortions produced by a needle
%%% in 3D space and fields produced by a simple coil
%%% configuration by Biot Savarts to compensate that field in the same space

%%% Parameter Settings
Center_Freq = 127;% spectrometer frequency in MHz
% Center_Freq = 63.5;% spectrometer frequency in MHz

chi_Water = -9.05*10^-6; % (Vol Susc : Schenk)
FOV = 40; % FOV (in mm)
FOV_wraparoundpad = FOV+FOV/2;

Res = 0.1;  %Resolution (im mm/voxel)
%  Res = 1;  %Resolution (im mm/voxel)

Res_CurElem = 0.4;  %Resolution (im mm)
gamma = 4*pi*1e-4;   % Induction constant gamma = mu0 = 4*pi*1e-7  Tm/A or  4*pi*1e-4 Tmm/A 
Current = 1;  % Amperes

%%%% Needle/Stylet Design: Volume Susceptiblities = Mass Susceptibly * Density (Kg/m3)

% Nitinol = 245*10^-6;  (Vol Susc : Schenk)
% Stainless_steel = 1500*10^-6; Koch 2013 
% Stainless_steel_316 = 1500*10^-6; Koch 2013  
% Titanium = 182*10^-6; (Vol Susc : Schenk)
% Cobalt_Chromium = 900*10^-6; (Vol Susc : Koch)
% Brass = -16*10^-6;
% chi_Air = 0.37*10^-6;  %(Vol Susc : Schenk)
% Bismuth = -165*10^-6;
% Pyrolytic_carbon = -204*10^-6;

%OD-ID
%10G: 3.404- 2.692
%12G: 2.769- 2.159
%14G: 2.108 - 1.6

  chi_Material  = 182*10^-6;% 
%   chi_Material  = 1500*10^-6;% 
%   chi_Material  = -204*10^-6;% 

%%%%%% 8G Needle
%    needle_od = 3.7;  % in mm
%    needle_id = 3;  % % in mm

%%%%%% 10G Needle
%     needle_od = 4.2;  % in mm
%     needle_id = 3.4;  % % in mm

% needle_od = 6.35;  % in mm
% needle_id = 4.57;  % % in mm

%%%%% 3 mm/2 mm Needle
% needle_od = 3;  % in mm
% needle_id = 2;  % % in mm

%%%%% 4 mm/3 mm Needle
needle_od = 4;  % in mm
needle_id = 3;  % % in mm
wire_dia = 0.4;
tip_clr = 1.3;
% Gauge_str = '18_gauge.stl';
res_pts = 14;


needle_wallthk = (needle_od-needle_id)/2;  % % in mm
needle_length = 100; % mm;
% needle_crop_length = 40 ; % mm 
needle_crop_length = FOV ; % mm 
needle_centroid = [0;0;0];  
bevel_angle = 30;

%%% For a closed Needle, the wall thickness at the tip is the same as
%%% needle_wallthk.  Therefore, wall thickness  in the Needle_X direction
%%% is needle_wallthk/sind(bevel_angle). For example, if wall thickness is
%%% 0.35 mm and bevel angle is 30, Needle_X wall thickness is  0.35/sind(30) = 0.7 mm
%%% With a tip clearance of 3.4/2/tand(30) = 2.95 mm, total clearance for
%%% coil just inside the needle tip is 2.95+0.7+wire_dia/2 = 3.85 mm.

% Needle Setup

clear Needle_Grid Needle_Grid_cropped Needle_Grid_cropped_padded  Needle_Grid_cropped_padded_rotated
clear Needle_Coord  Needle_Coord_MagRot
clear Mask Mask_Magnet_rotated


ND = round(FOV/Res); % Number of points

%%% 10 Gauge 3.4/2.7 mm
%  [Needle_Grid] = VOXELISE(needle_length/Res,needle_od/Res,needle_od/Res,'STL_Files/Needle.stl','xyz');
%   [Needle_Grid] = VOXELISE(round(needle_length/Res),linspace(0,needle_od,36),linspace(0,needle_od,36),'STL_Files/Needle_10G_30DegBevel.stl','xyz');
% Needle_Grid = Needle_Grid(:,2:35,2:35);

%%%% 3mm/2mm Needle
%%% [Needle_Grid] = VOXELISE(needle_length/Res,needle_od/Res,needle_od/Res,'STL_Files/Needle_3mm2mm_30degBevel.stl','xyz');
% [Needle_Grid] = VOXELISE(round(needle_length/Res),linspace(0,needle_od,32),linspace(0,needle_od,32),'STL_Files/Needle_3mm2mm_30degBevel.stl','xyz');
% Needle_Grid = Needle_Grid(:,2:end-1,2:end-1);

%%%% 4mm/3mm Needle
%%[Needle_Grid] = VOXELISE(needle_length/Res,needle_od/Res,needle_od/Res,'STL_Files/Needle_3mm2mm_30degBevel.stl','xyz');
% [Needle_Grid] = VOXELISE(round(needle_length/Res),linspace(0,needle_od,42),linspace(0,needle_od,42),'STL_Files/Needle_4mm3mm_30degBevel.stl','xyz');
% Needle_Grid = Needle_Grid(:,2:end-1,2:end-1);

%%% Different Gauges : Use the nearest Even number of voxels
%%%   [Needle_Grid] = VOXELISE(needle_od/Res,needle_od/Res,needle_length/Res,'STL_Files/12_gauge.stl','xyz');
%  [Needle_Grid] = VOXELISE(linspace(0,needle_od,res_pts),linspace(0,needle_od,res_pts),round(needle_length/Res),strcat('STL_Files/',Gauge_str),'xyz');
%  Needle_Grid = permute(Needle_Grid,[3,1,2]); Needle_Grid = flip(Needle_Grid,1);
%  Needle_Grid = Needle_Grid(:,2:end-1,2:end-1);


%%% Implant
% [Needle_Grid] = VOXELISE(400,400,400,'STL_Files/knee_metal_decimated.stl','xyz');





[XN,YN,ZN] = size(Needle_Grid);
Needle_Grid_cropped = Needle_Grid(XN-(needle_crop_length/Res-1):XN,:,:);
[XN_Cropped,YN_Cropped,ZN_Cropped] = size(Needle_Grid_cropped);


%%%%% Temoporary change for flat needle #$#######################
% Needle_Grid_cropped(end-80:end,:,:) = Needle_Grid_cropped(1:81,:,:);
%%%%% Temoporary change ***************########################


%%%% Pad Y and Z Dimensions to ND
Needle_Grid_cropped_padded = padarray(Needle_Grid_cropped,[0 ,(ND-YN_Cropped)/2 ,(ND-ZN_Cropped)/2 ]);



%%% Pad to a larger grid to avoid wraparound; In the needle x direction, pad only at the bottom. Store ND for restoring later.
needle_pad_bottom = FOV_wraparoundpad-needle_crop_length; % mm
y_pad = (FOV_wraparoundpad-FOV)/2;
z_pad = (FOV_wraparoundpad-FOV)/2;
Needle_Grid_cropped_padded = padarray(Needle_Grid_cropped_padded,[needle_pad_bottom/Res,0,0],'post');
Needle_Grid_cropped_padded = padarray(Needle_Grid_cropped_padded,[0,y_pad/Res,z_pad/Res]);

ND_Orig = ND;
ND = round(FOV_wraparoundpad/Res);


%%%% Create Mask
Mask = zeros(size(Needle_Grid_cropped_padded));
Mask_Tip = zeros(size(Needle_Grid_cropped_padded));


for n = 1: size(Needle_Grid_cropped_padded,2)
[row,col] = find(squeeze(Needle_Grid_cropped_padded(:,n,:)));
min_row = min(row);
max_row = max(row);
min_col = min(col);
max_col = max(col);

Mask(min_row:max_row,n,min_col:max_col) = 1;
% Mask_Tip(max_row-round(needle_wallthk/Res):max_row,n,min_col:max_col) = 1;
% if bevel_angle ~= 
Mask_Tip(max_row-round(needle_wallthk/Res/sind(bevel_angle)):max_row,n,min_col:max_col) = 1;
end


%%%%% Temporary change for flat needle closed on both ends ########################
%   Needle_Grid_cropped_padded(1:200,:,:) = Needle_Grid_cropped_padded(401:end,:,:);
%   Mask(1:200,:,:) = Mask(401:end,:,:);
% 
% Needle_Grid_cropped_padded(395:400,:,:) = Mask(395:400,:,:);
%  Needle_Grid_cropped_padded(200:205,:,:) = Mask(200:205,:,:);
% %%%%% Temoporary change ***************########################


%%%%% Temoporary change for Closed Bevelled needle ########################
% Needle_Grid_cropped_padded = Needle_Grid_cropped_padded | Mask_Tip;
%%%%% Temoporary change ***************########################



%%%%%% Needle Grid is defined in the UNROATED Scanner Geometry
Needle_Grid_cropped_padded_rotated = zeros(size(Needle_Grid_cropped_padded));
[r,c,d] = ind2sub(size(Needle_Grid_cropped_padded),find(Needle_Grid_cropped_padded));
Needle_Coord = [r,c,d];
Needle_Coord = Needle_Coord - (ND+1)/2;
clear r c d

Mask_Magnet_rotated = zeros(size(Mask));
[r,c,d]= ind2sub(size(Mask),find(Mask));
Mask_Coord =[r,c,d];
Mask_Coord = Mask_Coord - (ND+1)/2;
clear r c d


%%% ---------------------------IMPORTANT -------------------------------%%%
%%% --------------------------------------------------------------------%%%
%%%  We Need three angles to define the orientation of the needle in the 
%%%  Magnet coordinate frame . The reference (Default) Needle orientation
%%%  is along the X (Up-Down) axis with the tip pointed down. If there is
%%%  a bevel, the shorter end is along the + Y axis.

%%% Angles are defined in the LEFT-HANDED system, i.e, Clockwise rotations
%%% are + and counterclockwise rotations are -.

%%% For eg, for a needle along the Z (FH) axis, with tip along +Z (Foot)
%%% Ang_X = 0;
%%% Ang_Y = 90;
%%% Ang_Z = 0;

%%% For a needle along the Z (FH) axis, with tip along -Z (Head)
%%% Ang_X = 0;
%%% Ang_Y = -90;
%%% Ang_Z = 0;

%%% For a needle along the Y (LR) axis, with tip along +Y (Patient Left)
%%% Ang_X = 0;
%%% Ang_Y = 0;
%%% Ang_Z = -90;

%%% For a needle along the Y (LR) axis, with tip along -Y (Patient Right)
%%% Ang_X = 0;
%%% Ang_Y = 0;
%%% Ang_Z = 90;

%%% For a needle along the X (AP) axis, with tip along -X (Patient Posterior)
%%% Ang_X = 0;
%%% Ang_Y = 0;
%%% Ang_Z = 0;

Ang_X = 0;
Ang_Y = -90;
Ang_Z = 0;


% RotX = [ 1, 0, 0; 0 cosd(Ang_X)  -sind(Ang_X); 0 sind(Ang_X) cosd(Ang_X)];
% RotY = [cosd(Ang_Y), 0 , sind(Ang_Y); 0,1 ,0 ; -sind(Ang_Y) , 0 , cosd(Ang_Y)];
% RotZ = [cosd(Ang_Z), -sind(Ang_Z), 0; sind(Ang_Z), cosd(Ang_Z) , 0; 0, 0,1];
% Mag_MagRot_Transform = RotZ * RotY * RotX (Premultiplication Factor)

%%%% Left Handed Rotation
% Mag_MagRot_Transform = [cosd(Ang_Y)*cosd(Ang_Z), cosd(Ang_X)*sind(Ang_Z)+ sind(Ang_X)*sind(Ang_Y)*cosd(Ang_Z), sind(Ang_X)*sind(Ang_Z)- cosd(Ang_X)*sind(Ang_Y)*cosd(Ang_Z) ;...
%           -cosd(Ang_Y)*sind(Ang_Z), cosd(Ang_X)*cosd(Ang_Z)- sind(Ang_X)*sind(Ang_Y)*sind(Ang_Z), sind(Ang_X)*cosd(Ang_Z)+ cosd(Ang_X)*sind(Ang_Y)*sind(Ang_Z);...    
%                 sind(Ang_Y)  ,                    -sind(Ang_X)*cosd(Ang_Y)  ,                                 cosd(Ang_X)*cosd(Ang_Y)];

%%%% Right Handed Rotation
Mag_MagRot_Transform = [cosd(Ang_Y)*cosd(Ang_Z), -cosd(Ang_X)*sind(Ang_Z)+ sind(Ang_X)*sind(Ang_Y)*cosd(Ang_Z), sind(Ang_X)*sind(Ang_Z)+cosd(Ang_X)*sind(Ang_Y)*cosd(Ang_Z) ;...
           cosd(Ang_Y)*sind(Ang_Z), cosd(Ang_X)*cosd(Ang_Z)+sind(Ang_X)*sind(Ang_Y)*sind(Ang_Z), -sind(Ang_X)*cosd(Ang_Z)+ cosd(Ang_X)*sind(Ang_Y)*sind(Ang_Z);...    
                 -sind(Ang_Y)  ,                    sind(Ang_X)*cosd(Ang_Y)  ,                                 cosd(Ang_X)*cosd(Ang_Y)];
           
%%%% First Transform to UnRotated Magnet Axes
Needle_Mag_Transform = [ -1 0 0; 0 1 0; 0 0 1];  
Needle_Coord_Mag = Needle_Coord * Needle_Mag_Transform;
Mask_Coord_Mag = Mask_Coord * Needle_Mag_Transform;

%%%% Then perform Rotations to Orient Needle along any direction (NOTE : We
%%%% are Premultiplying, there transpose is required)

Needle_Coord_MagRot =  Mag_MagRot_Transform * Needle_Coord_Mag' ;
Mask_Coord_MagRot =  Mag_MagRot_Transform * Mask_Coord_Mag'  ;

Needle_Coord_MagRot = Needle_Coord_MagRot';
Mask_Coord_MagRot = Mask_Coord_MagRot';

% Needle_Coord_MagRot = Needle_Coord_Mag * Mag_MagRot_Transform ;
% Mask_Coord_MagRot = Mask_Coord_Mag * Mag_MagRot_Transform ;


for n = 1 : length(Needle_Coord_MagRot)
    
    source_pt = round([ Needle_Coord(n,1),Needle_Coord(n,2),Needle_Coord(n,3)]+(ND+1)/2);
    target_pt = round([ Needle_Coord_MagRot(n,1),Needle_Coord_MagRot(n,2),Needle_Coord_MagRot(n,3)]+(ND+1)/2);
        
    Needle_Grid_cropped_padded_rotated(target_pt(1),target_pt(2),target_pt(3)) = Needle_Grid_cropped_padded(source_pt(1),source_pt(2),source_pt(3));    
end  
Needle_Grid_cropped_padded_rotated = single(Needle_Grid_cropped_padded_rotated);


for n = 1 : length(Mask_Coord_MagRot)
    
    source_pt = round([ Mask_Coord(n,1),Mask_Coord(n,2),Mask_Coord(n,3)]+(ND+1)/2);
    target_pt = round([ Mask_Coord_MagRot(n,1),Mask_Coord_MagRot(n,2),Mask_Coord_MagRot(n,3)]+(ND+1)/2);
        
    Mask_Magnet_rotated(target_pt(1),target_pt(2),target_pt(3)) = Mask(source_pt(1),source_pt(2),source_pt(3));    
end   

%%% Overwrite Mask and ensure ODD size
Mask = Mask_Magnet_rotated(2:ND,2:ND,2:ND);
% Mask = Mask_Magnet_rotated;


% Plotting
ds = 8;
figure; subplot(1,3,1) ; scatter3(Needle_Coord(1:ds:end,3),Needle_Coord(1:ds:end,2),Needle_Coord(1:ds:end,1),'.');
xlim([-200,200]);ylim([-200,200]);zlim([-200,200]);
xlabel('Needle_Z');ylabel('Needle_Y');zlabel('Needle_X'); 
title(' Needle in the Needle Coordinate System'); 

subplot(1,3,2) ; scatter3(Needle_Coord_Mag(1:ds:end,3),Needle_Coord_Mag(1:ds:end,2),Needle_Coord_Mag(1:ds:end,1),'.');
xlim([-200,200]);ylim([-200,200]);zlim([-200,200]);
xlabel('Magnet_Z');ylabel('Magnet_Y');zlabel('Magnet_X'); 
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
title(' Needle in the UnRotated Magnet Coordinate System'); 


subplot(1,3,3) ; scatter3(Needle_Coord_MagRot(1:ds:end,3),Needle_Coord_MagRot(1:ds:end,2),Needle_Coord_MagRot(1:ds:end,1),'.');
xlim([-200,200]);ylim([-200,200]);zlim([-200,200]);
xlabel('Magnet_Z');ylabel('Magnet_Y');zlabel('Magnet_X'); 
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
title(' Needle in the Rotated Magnet Coordinate System'); 


clear Mask_Coord Mask_Coord_Mag Mask_Coord_MagRot
clear Needle_Coord Needle_Coord_Mag 



%%  Simulate Needle B0  %%%%%%%%%%%%%
%%%%% NOTE : Dimension 1 is X : Up-Down in the Scanner,with Up +
%%%%% NOTE : Dimension 2 is Y : Left-Right in the Scanner
%%%%% NOTE : Dimension 3 is Z : Foot-Head in the Scanner
%%%%% NOTE : X-Y Plane is Axial : eg (Chi(:,:,200,4))
%%%%% NOTE : X-Z Plane is Sagittal : eg (Chi(:,200,:,4))
%%%%% NOTE : Y-Z Plane is Coronal : eg (Chi(200,:,:,4))

x = linspace(-FOV_wraparoundpad/2,FOV_wraparoundpad/2, ND);
y = linspace(-FOV_wraparoundpad/2, FOV_wraparoundpad/2, ND);
z = linspace(-FOV_wraparoundpad/2, FOV_wraparoundpad/2, ND);
[Y, X, Z] = meshgrid(y,-x,-z);


%%%%%%%%%%%%%% Define Chi Distriution  %%%%%%%%%%%

id = find(Needle_Grid_cropped_padded_rotated == 1);
S_Object = ones(size(X))* chi_Water; %%%% Assume Water as surrounding Medium
S_Object(id) = chi_Material;
% mask =  Needle_Grid_cropped_padded; %%% Mask 

%%%% If Channel inside is AIR
% Air = Mask_Magnet_rotated-Needle_Grid_cropped_padded_rotated;
% id_Air = find(Air == 1);
% S_Object(id_Air) = chi_Air;


E=single(S_Object); 
clear R_object S_object;
clear Needle_Grid Needle_Grid_cropped Needle_Grid_cropped_padded;


%%%%%%%%%%%%%% Performing B0 Estimation  %%%%%%%%%%%

% ensure odd matrix size
if mod(ND,2)~=1 
    E=E(2:ND,2:ND,2:ND); 
    X=X(2:ND,2:ND,2:ND); 
    Y=Y(2:ND,2:ND,2:ND); 
    Z=Z(2:ND,2:ND,2:ND); 
end
ND_new = size(E,1);


%%%%%% Define 3D Spatial Domain Chi Matrix %%%%%%%%%%%%%%%%
Chi(:,:,:,4) = E;
Chi(:,:,:,1:3) = cat(4,X,Y,Z);
Chi=single(Chi);


%%%%%% Define 3D Chi_K space %%%%%%%%%%%%%%%%
Chi_K = zeros(ND_new,ND_new,ND_new,4);
Chi_Kx = linspace(-(ND_new-1)/2,(ND_new-1)/2,ND_new);

%  [ Chi_K_Y, Chi_K_X, Chi_K_Z] = meshgrid(Chi_Kx,Chi_Kx,-Chi_Kx);
 [ Chi_K_Y, Chi_K_X, Chi_K_Z] = meshgrid(Chi_Kx,-Chi_Kx,-Chi_Kx);

Chi_K(:,:,:,1:3) = cat(4,Chi_K_X, Chi_K_Y, Chi_K_Z);
Chi_K=single(Chi_K);



%%%%% calculate B-field %%%%
pause(0.001);
Chi_K(:,:,:,4)=fftshift(fftn(Chi(:,:,:,4)));
B(:,:,:,1:3)= Chi(:,:,:,1:3);
% clear Chi;


% B(k) = Chi(k)*(1/3 - kz^2/(kz^2+kx^2+ky^2)) 
B_K(:,:,:,1:3)= Chi_K(:,:,:,1:3); % now physics!
% B_K(:,:,:,4)= Chi_K(:,:,:,4).*(1/3-Chi_K(:,:,:,1).^2./(Chi_K(:,:,:,1).^2+Chi_K(:,:,:,2).^2+Chi_K(:,:,:,3).^2)); 
  B_K(:,:,:,4)= Chi_K(:,:,:,4).*(1/3-Chi_K(:,:,:,3).^2./(Chi_K(:,:,:,1).^2+Chi_K(:,:,:,2).^2+Chi_K(:,:,:,3).^2)); 
% singularity --> substitude with dc-offset
B_K((size(B_K,1)+1)/2,(size(B_K,2)+1)/2,(size(B_K,3)+1)/2,4) = 0;
clear Chi_K;


pause(0.001);
B(:,:,:,4)=(ifftn(ifftshift(B_K(:,:,:,4))));
Deltaf = Center_Freq*10^6*real(B(:,:,:,4)); % frequency shift distribution in Hz


%%%%% Plotting %%%%%%%%%%%%%
min_Df = min(Deltaf(:));
max_Df = max(Deltaf(:));


%%%%% Restore Original Grid Size of 40 mm
Final_crop_range = ND/2-ND_Orig/2:ND/2+ND_Orig/2-1;
Deltaf_Final = Deltaf(Final_crop_range,Final_crop_range,Final_crop_range);
Mask_Final = Mask(Final_crop_range,Final_crop_range,Final_crop_range);
Needle_Grid_Final = Needle_Grid_cropped_padded_rotated(Final_crop_range,Final_crop_range,Final_crop_range);
ND = ND_Orig;


% vuThreePaneViewer(Deltaf);
clim = [-1000,1000];

figure; colormap jet; 
subplot(2,2,1); imagesc(Deltaf_Final(:,:,ND/2), clim);
set(gca, 'YDir','normal');title('AXIAL MIDSLICE, LOOKING FROM PATIENT FEET');
axis image; xlabel('Magnet Y');ylabel('Magnet X');

subplot(2,2,2); imagesc(squeeze(Deltaf_Final(:,ND/2,:)),clim);
set(gca, 'YDir','normal'); set(gca, 'XDir','reverse');title('SAGITTAL MIDSLICE, LOOKING FROM PATIENT LEFT');axis image
xlabel('Magnet Z');ylabel('Magnet X');hold on;
quiver(100,350,-50,0,'Color',[0,0,0],'Linewidth',3,'MaxHeadSize',0.8)

subplot(2,2,3); imagesc(squeeze(Deltaf_Final(ND/2,:,:)), clim);
set(gca, 'XDir','reverse'); title('CORONAL MIDSLICE, LOOKING FROM TOP'); 
axis image;xlabel('Magnet Z');ylabel('Magnet Y'); hold on;
quiver(100,50,-50,0,'Color',[0,0,0],'Linewidth',3,'MaxHeadSize',0.8)
colorbar


clear Chi_K_X Chi_K_Y Chi_K_Z B_K
clear B Chi
clear X_Object Y_Object Z_Object E S_Object
clear Mask Mask_Magnet_rotated mask 
clear X Y Z x y z
clear ND_new
% clear Needle_Grid_cropped_padded_rotated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Induction curve (Make sure this is in the magnet coordinate system)

if ND ~= ND_Orig
    ND = ND_Orig;
end

%%%% Define Loops in Needle Coordinate system First
%%%% :: Needle_X is along Needle length with Bevelled Tip Positive 
%%%  :: Theta_N is the radial angle in the Needle coordinate system.
%%%  :: Theta_N = 0 when along Needle_Y
%%%  :: Theta_N = 90 when along Needle_Z

%%% NOTE : We use only a single angle instead of three in the Needle 
%%% coordinate system as we are placing conductors mostly parallel to the 
%%% needle axis.

% needle_od = 3.4;  % in mm
% needle_id = 2.7;  % in mm
% needle_centroid = [0;0;0]; 


%  wire_dia = 0.4; % 26G wire im mm .
% wire_dia = 0.25; % 30G wire im mm (Used for 3mm/2mm Simulation)
% insert_od = needle_id - 0.5;  %%%% Clearance
insert_od = needle_id ;  %%%% Clearance

needle_top_edge = needle_centroid(1)-FOV/2;
needle_bottom_edge = needle_crop_length - FOV_wraparoundpad/2 ;

spearloop_flag = 0;
Ang_Loop_Flag = 0;
y_cutoff_flag = 0;
speardoublewireloop_flag = 0;

% prompt = ' Please enter a bevel angle (deg, 90 for Flat Tip) :  ';
% bevel_ang = input(prompt); % Bevel Angle of coil
bevel_ang = 40;

rad_full = (insert_od/2-wire_dia/2);
% prompt = ' Please enter a radius fraction (1 for 100%) :  ';
% rad_fraction = input(prompt);

rad_fraction = 1;
rad = rad_full*rad_fraction;

% prompt = ' Please enter a tip clearance (mm) :  ';
% tip_clr = input(prompt); % Clearance of coil from needle tip in mm

prompt = ' Please enter a CN angle (degrees) :  ';
theta_N = input(prompt);    

prompt = ' Please enter Number of Turns :  ';
nturns = input(prompt);  

prompt = ' Please choose With (1) or Without (0) Leads :  ';
leads = input(prompt);  

       
switch menu('Choose a coil geometry: Angle is defined with respect to +Y axis in Needle Frame',...
        'XN_NTurnLoop',...
        'XN_NTurnAngledLoop',...
        'XN_HelicalLoop',...        
        'XN_SquareLoop',...  
        'XN_CurvedLoop',...
        'XN_SplitSpearLoop Leg 1',...
        'XN_SplitSpearLoop Leg 2',...
        'XN_Half loop',...
        'XN_Half Rounded loop',...
        'XN_Angled Flatloop',...
        'XN_Antiparallel Helmholtz loop',...
        'Q Coil',...
        'F Coil',...
        'Double Split Loop Leg 1',...
        'Double Split Loop Leg 2')
    
    case 1  % N Turn Loop' ------------------------------------------//
                          
        %%%% Arc or Full Circle
        
        %%% For Closed Needles
%         total_wall_thk_y = (needle_wallthk/cosd(bevel_ang)) + needle_wallthk;  %%%%( Wall thickness in Y direction at Tip clearance)
          
         %%% For Open Needles
         total_wall_thk_y = needle_wallthk;  %%%% ( Wall thickness in Y direction at Tip clearance)
         
         
         air_space = tand(bevel_ang)*tip_clr - total_wall_thk_y;                        
         y_cutoff =  air_space-rad;
         y_cutoff_flag = 1;
         
         npts = (abs(nturns)*20)+1;
                           
          if nturns >= 1 
             phi_w = linspace(0, nturns*2*pi,npts);
              pitch = wire_dia; %%% Pitch of the turns in mm ( distance between turns)
%              pitch = 0; %%% Pitch of the turns in mm ( distance between turns)
             tip_clr = ones(size(phi_w,2),1)*tip_clr + linspace(0,pitch*nturns,size(phi_w,2))';
          else
              if air_space < 2*rad  %%%% ie if Coil is in the bevel region, the chop off coil
                  phi_cutoff = asin(y_cutoff/rad);
                  phi_w = linspace(0-phi_cutoff, pi+phi_cutoff,npts);
              else
                  phi_w = linspace(0, nturns*2*pi,npts);
              end                 
          end
           
          
         L =[  needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr   rad.*cos(phi_w')   rad.*sin(phi_w') ];
                      
          % Extra Diagonal Connection at the beginning of the loop
%           L = [L(1,1),L(1,2),L(1,3)-rad*2;L];
             
         
    case 2  %  Angled Loop'  ------------------------------------------//
        
        %   --> X Axis         Y Axis 
        %--------------------\   ^
        %     Needle          \  |     
        %----------------------\
        
        
        a = rad/sind(bevel_ang);  %%% Major Axis along Needle Y (Across bevel slant)
        b = rad;  %%% Minor Axis along Needle Z ( Across Symmetric needle)
        npts = (abs(nturns)*20)+1;
        Ang_Loop_Flag = 1;
         
        %%%%% Add pi/2 angle here instead of Theta_N tranform at the end.
         phi_w = linspace(0, nturns*2*pi,npts)+(theta_N*pi/180); 
               
        if nturns >= 1
         pitch = wire_dia; %%% Pitch of the turns in mm ( distance between turns)
         tip_clr = ones(size(phi_w,2),1)*tip_clr + linspace(0,pitch*abs(nturns),size(phi_w,2))';
        end

        %%% Define Ellipse
        r = (ones(1,npts)*a*b)./((b^2*cos(phi_w).^2) + (a^2*sin(phi_w).^2)).^0.5;           
        L =[ zeros(size(phi_w,2),1)  r'.*cos(phi_w')  r'.*sin(phi_w')];
        
        BAng_X = 0;
        BAng_Y = 0;
        BAng_Z = 90-bevel_ang;


        L_Bevel_Transform = [cosd(BAng_Y)*cosd(BAng_Z), cosd(BAng_X)*sind(BAng_Z)+ sind(BAng_X)*sind(BAng_Y)*cosd(BAng_Z), sind(BAng_X)*sind(BAng_Z)- cosd(BAng_X)*sind(BAng_Y)*cosd(BAng_Z) ;...
          -cosd(BAng_Y)*sind(BAng_Z), cosd(BAng_X)*cosd(BAng_Z)- sind(BAng_X)*sind(BAng_Y)*sind(BAng_Z), sind(BAng_X)*cosd(BAng_Z)+ cosd(BAng_X)*sind(BAng_Y)*sind(BAng_Z);...    
                sind(BAng_Y)  ,                    -sind(BAng_X)*cosd(BAng_Y)  ,                                 cosd(BAng_X)*cosd(BAng_Y)];

        L = L * L_Bevel_Transform;   
        L(:,1) = L(:,1) + (needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr);
        
                      
        % Extra Diagonal Connection       
        %          L(end+1,:) = [L(end,1)-2, L(end,2),L(end,3)];
        %          L(end+1,:) = [L(end,1), L(end,2),L(end,3)-rad*2];

        % Extra Diagonal Connection at the beginning of the loop
        if mod(nturns,1) == 0 
            if theta_N == 0
                    L = [max(L(:,1))-pitch,L(1,2)-rad*2,L(1,3);L];
            elseif theta_N ==90
                 L = [L(1,1),L(1,2),L(1,3)-rad*2;L];
            end
        end

   
                
    case 3  % XN_HelicalLoop'  ------------------------------------------//
     
         wire_dia = 0.4; % 32G wire im mm    
         %%%% N Full Turn s 
         NTurns = 1.5;        
         pitch = 1; %%% Pitch of the turns in mm ( distance between turns)
         x_step = linspace((insert_od-wire_dia)/tand(bevel_ang)+wire_dia+pitch,0,20);              
         phi_w = linspace(0, 2*pi*NTurns,20);
         
         L =[  needle_top_edge     rad*cos(phi_w(1))  rad*sin(phi_w(1)) ;
              (needle_bottom_edge-tip_clr)*ones(size(phi_w,2),1)-x_step'     rad*cos(phi_w')  rad*sin(phi_w')  ;  
               needle_top_edge    rad*cos(phi_w(end))   rad*sin(phi_w(end))  ];
           
                                                  
    case 4   %  XN_SquareLoop  ------------------------------------------//         
                
        L = [  needle_bottom_edge-tip_clr   rad    0;   % x y z;
               needle_bottom_edge-tip_clr  -rad   0;  ];    % x y z; 
           
               
    case 5   %  XN_Curved loop  ------------------------------------------//
        
         rho_w = linspace(0, pi,10);
         rad = insert_od/2-wire_dia/2 ; %%% Radius of wire bend in mm      
         L =[  needle_top_edge  rad 0;
               needle_bottom_edge-tip_clr + rad*sin(rho_w')  rad*cos(rho_w')  zeros(10,1)  ;  
               needle_top_edge  -rad  0  ];
        
        
    case 6   % Split Spear Loop : Leg 1  --------------------------------//
              

        nturns = 0.5;       
        a = rad/sind(bevel_ang);  %%% Major Axis along Needle Y (Across bevel slant)
        b = rad;  %%% Minor Axis along Needle Z ( Across Symmetric needle)
        npts = (abs(nturns)*20)+1;
        Ang_Loop_Flag = 1;
        spearloop_flag = 1;
         
        %%%%% Add pi/2 angle here instead of Theta_N tranform at the end.
         phi_w = linspace(0, nturns*2*pi,npts)+(theta_N*pi/180); 

        %%% Define Ellipse
        r = (ones(1,npts)*a*b)./((b^2*cos(phi_w).^2) + (a^2*sin(phi_w).^2)).^0.5;           
        L =[ zeros(size(phi_w,2),1) r'.*cos(phi_w')   r'.*sin(phi_w')];
        
        BAng_X = 0;
        BAng_Y = 0;
        BAng_Z = 90-bevel_ang;

        L_Bevel_Transform = [cosd(BAng_Y)*cosd(BAng_Z), cosd(BAng_X)*sind(BAng_Z)+ sind(BAng_X)*sind(BAng_Y)*cosd(BAng_Z), sind(BAng_X)*sind(BAng_Z)- cosd(BAng_X)*sind(BAng_Y)*cosd(BAng_Z) ;...
          -cosd(BAng_Y)*sind(BAng_Z), cosd(BAng_X)*cosd(BAng_Z)- sind(BAng_X)*sind(BAng_Y)*sind(BAng_Z), sind(BAng_X)*cosd(BAng_Z)+ cosd(BAng_X)*sind(BAng_Y)*sind(BAng_Z);...    
                sind(BAng_Y)  ,                    -sind(BAng_X)*cosd(BAng_Y)  ,                                 cosd(BAng_X)*cosd(BAng_Y)];

        L = L * L_Bevel_Transform;   
        L(:,1) = L(:,1) + (needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr);
               
        
             
  case 7   % Split Spear Loop : Leg 2  ---------------------------------//
              
        nturns = -0.5;       
        a = rad/sind(bevel_ang);  %%% Major Axis along Needle Y (Across bevel slant)
        b = rad;  %%% Minor Axis along Needle Z ( Across Symmetric needle)
        npts = (abs(nturns)*20)+1;
        Ang_Loop_Flag = 1;
        spearloop_flag = 2;
         
        %%%%% Add pi/2 angle here instead of Theta_N tranform at the end.
         phi_w = linspace(0, nturns*2*pi,npts)+(theta_N*pi/180); 
       

        %%% Define Ellipse
        r = (ones(1,npts)*a*b)./((b^2*cos(phi_w).^2) + (a^2*sin(phi_w).^2)).^0.5;           
        L =[ zeros(size(phi_w,2),1) r'.*cos(phi_w')   r'.*sin(phi_w')];
        
        BAng_X = 0;
        BAng_Y = 0;
        BAng_Z = 90-bevel_ang;

        L_Bevel_Transform = [cosd(BAng_Y)*cosd(BAng_Z), cosd(BAng_X)*sind(BAng_Z)+ sind(BAng_X)*sind(BAng_Y)*cosd(BAng_Z), sind(BAng_X)*sind(BAng_Z)- cosd(BAng_X)*sind(BAng_Y)*cosd(BAng_Z) ;...
          -cosd(BAng_Y)*sind(BAng_Z), cosd(BAng_X)*cosd(BAng_Z)- sind(BAng_X)*sind(BAng_Y)*sind(BAng_Z), sind(BAng_X)*cosd(BAng_Z)+ cosd(BAng_X)*sind(BAng_Y)*sind(BAng_Z);...    
                sind(BAng_Y)  ,                    -sind(BAng_X)*cosd(BAng_Y)  ,                                 cosd(BAng_X)*cosd(BAng_Y)];

        L = L * L_Bevel_Transform;   
        L(:,1) = L(:,1) + (needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr);
        
              
        
    case 8  % XN_Half Pointed tip loop'  -------------------------------//
     
        %%%% N Full Turn s 
         NTurns = 0.5;
         x_step1 = linspace((insert_od-wire_dia)/tand(bevel_ang),0,5);
         x_step = [x_step1 x_step1(end-1:-1:1)];
              
         phi_w = linspace(0, 2*pi*NTurns,length(x_step));
         L =[  needle_top_edge             rad*cos(phi_w(1))  rad*sin(phi_w(1)) ;
              (needle_bottom_edge-tip_clr)*ones(size(phi_w,2),1)-x_step'         rad*cos(phi_w')  rad*sin(phi_w')  ;  
               needle_top_edge              rad*cos(phi_w(end))   rad*sin(phi_w(end))  ];
           
    case 9  % XN_Half Rounded Tip loop'  ------------------------------------//
     
         %%%% N Full Turns 
         NTurns = 0.5;
         rho_w = linspace(0, pi,10);
         curve_radius = insert_od/2-wire_dia/2 ; %%% Radius of wire bend in mm
                 
         phi_w = linspace(0, 2*pi*NTurns,length(rho_w));
         L =[  needle_top_edge                                        rad*cos(phi_w(1))   rad*sin(phi_w(1)) ;
              (needle_bottom_edge-tip_clr)+ curve_radius*sin(rho_w')  rad*cos(phi_w')     rad*sin(phi_w')  ;  
               needle_top_edge                                        rad*cos(phi_w(end)) rad*sin(phi_w(end))  ];
           
           
    case 10   % Triangular Angled Flat Loop  ----------------------------------------//
              
         %%% Tip loop will follow bevel angle at a clearance distance
%          x_step = linspace((insert_od-wire_dia)/tand(bevel_ang),0,10);
%          y_step = linspace(insert_od/2-wire_dia/2,-(insert_od/2-wire_dia/2),10);
%          
%          L =[ (needle_bottom_edge-tip_clr)* ones(10,1)-x_step'  y_step'  zeros(10,1) ;   ];
         
         
         x_step = linspace((2*(rad-wire_dia/2))/tand(bevel_ang),0,10);
         y_step = linspace((rad-wire_dia/2),-(rad-wire_dia/2),10);
         
         L =[ (needle_bottom_edge-tip_clr)* ones(10,1)-x_step'  y_step'  zeros(10,1) ;   ];
         
         
    
    case 11  %XN_Antiparallel Helmholtz loop  -----------------------------//
         
        
         dist = 1 ; % Distance between Helmholtz loops  
         nturns_p = 1.5;
         nturns_ap = 2.5;
         phi_w = linspace(0,  nturns_p*2*pi,10);
         phi_w_ap = linspace(phi_w(end),phi_w(end)-nturns_ap*2*pi,10);
         
          L =[ needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr   rad.*cos(phi_w')   rad.*sin(phi_w');
              needle_bottom_edge*ones(size(phi_w_ap,2),1)-(tip_clr+dist)  rad*cos(phi_w_ap')  rad*sin(phi_w_ap')];
                  
%          L =[  needle_top_edge    rad*cos(phi_w(1))  rad*sin(phi_w(1)) ;
%               (needle_bottom_edge-tip_clr)*ones(size(phi_w,2),1)    0.5*rad*cos(phi_w')  0.5*rad*sin(phi_w') ;
%               (needle_bottom_edge-tip_clr-dist)*ones(size(phi_w_ap,2),1)  rad*cos(phi_w_ap')  rad*sin(phi_w_ap');
%               needle_top_edge     rad*cos(phi_w_ap(end))   rad*sin(phi_w_ap(end))  ];
          
     case 12  % Q coil  ------------------------------------------//
             
%           L =[ needle_bottom_edge-6-tip_clr  0   -rad;            
%               needle_bottom_edge-6-tip_clr  0   -rad+1;              
%              needle_bottom_edge-5-tip_clr  0   -rad+1;             
%              needle_bottom_edge-5-tip_clr  0    -rad; 
%              needle_bottom_edge-tip_clr 	0   -rad; 
%              needle_bottom_edge-tip_clr    0    rad ;
%               needle_bottom_edge-5-tip_clr  0   rad;
%               needle_bottom_edge-5-tip_clr  0   rad-1;
%               needle_bottom_edge-6-tip_clr    0   rad-1
%               needle_bottom_edge-6-tip_clr    0   rad ];
          
           L =[  needle_bottom_edge-5-tip_clr  -rad   -rad+1;             
             needle_bottom_edge-5-tip_clr  -rad    -rad; 
             needle_bottom_edge-tip_clr 	-rad   -rad; 
             needle_bottom_edge-tip_clr    -rad    rad ;
              needle_bottom_edge-5-tip_clr  -rad   rad;
              needle_bottom_edge-5-tip_clr  -rad   rad-1; ];
          
     case 13  %F coil  ------------------------------------------//
             
         L =[ needle_bottom_edge-10   -rad   rad ;
              needle_bottom_edge-10  -rad   -rad;
              needle_bottom_edge-15  -rad   -rad;
              needle_bottom_edge-15   -rad   rad
              needle_bottom_edge-10   -rad   rad];
          
          
    case 14   % Double Split Loop : Leg 1  --------------------------------//
              
        nturns = 0.5;       
        a = rad/sind(bevel_ang);  %%% Major Axis along Needle Y (Across bevel slant)
        b = rad;  %%% Minor Axis along Needle Z ( Across Symmetric needle)
        npts = (abs(nturns)*20)+1;
        Ang_Loop_Flag = 1;
        speardoublewireloop_flag = 1;
         
        %%%%% Add pi/2 angle here instead of Theta_N tranform at the end.
         phi_w = linspace(0, nturns*2*pi,npts)+(theta_N*pi/180); 

        %%% Define Ellipse
        r = (ones(1,npts)*a*b)./((b^2*cos(phi_w).^2) + (a^2*sin(phi_w).^2)).^0.5;           
        L =[ zeros(size(phi_w,2),1) r'.*cos(phi_w')   r'.*sin(phi_w')];
        
        BAng_X = 0;
        BAng_Y = 0;
        BAng_Z = 90-bevel_ang;

        L_Bevel_Transform = [cosd(BAng_Y)*cosd(BAng_Z), cosd(BAng_X)*sind(BAng_Z)+ sind(BAng_X)*sind(BAng_Y)*cosd(BAng_Z), sind(BAng_X)*sind(BAng_Z)- cosd(BAng_X)*sind(BAng_Y)*cosd(BAng_Z) ;...
          -cosd(BAng_Y)*sind(BAng_Z), cosd(BAng_X)*cosd(BAng_Z)- sind(BAng_X)*sind(BAng_Y)*sind(BAng_Z), sind(BAng_X)*cosd(BAng_Z)+ cosd(BAng_X)*sind(BAng_Y)*sind(BAng_Z);...    
                sind(BAng_Y)  ,                    -sind(BAng_X)*cosd(BAng_Y)  ,                                 cosd(BAng_X)*cosd(BAng_Y)];

        L = L * L_Bevel_Transform;   
        L(:,1) = L(:,1) + (needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr);
        
        L(1,3) =  L(1,3)+ wire_dia/4;
        L(end,3) =  L(end,3)+ wire_dia/4;
        
   case 15   % Double Split Loop : Leg 2  --------------------------------//
              
        nturns = 0.5;       
        a = rad/sind(bevel_ang);  %%% Major Axis along Needle Y (Across bevel slant)
        b = rad;  %%% Minor Axis along Needle Z ( Across Symmetric needle)
        npts = (abs(nturns)*20)+1;
        Ang_Loop_Flag = 1;
        speardoublewireloop_flag = 1;
         
        %%%%% Add pi/2 angle here instead of Theta_N tranform at the end.
         phi_w = linspace(0, nturns*2*pi,npts)+(theta_N*pi/180); 

        %%% Define Ellipse
        r = (ones(1,npts)*a*b)./((b^2*cos(phi_w).^2) + (a^2*sin(phi_w).^2)).^0.5;           
        L =[ zeros(size(phi_w,2),1) r'.*cos(phi_w')   r'.*sin(phi_w')];
        
        BAng_X = 0;
        BAng_Y = 0;
        BAng_Z = 90-bevel_ang;

        L_Bevel_Transform = [cosd(BAng_Y)*cosd(BAng_Z), cosd(BAng_X)*sind(BAng_Z)+ sind(BAng_X)*sind(BAng_Y)*cosd(BAng_Z), sind(BAng_X)*sind(BAng_Z)- cosd(BAng_X)*sind(BAng_Y)*cosd(BAng_Z) ;...
          -cosd(BAng_Y)*sind(BAng_Z), cosd(BAng_X)*cosd(BAng_Z)- sind(BAng_X)*sind(BAng_Y)*sind(BAng_Z), sind(BAng_X)*cosd(BAng_Z)+ cosd(BAng_X)*sind(BAng_Y)*sind(BAng_Z);...    
                sind(BAng_Y)  ,                    -sind(BAng_X)*cosd(BAng_Y)  ,                                 cosd(BAng_X)*cosd(BAng_Y)];

        L = L * L_Bevel_Transform;   
        L(:,1) = L(:,1) + (needle_bottom_edge*ones(size(phi_w,2),1)-tip_clr);
        
        L(1,3) =  L(1,3)+ wire_dia/4;
        L(end,3) =  L(end,3)+ wire_dia/4;
        
        L(:,3) = -L(:,3);
                          
end

figure; subplot(1,3,1); plot3(L(:,3),L(:,2),L(:,1),'r','linewidth',3); grid on;
xlim([-FOV/2-10 FOV/2+10]);ylim([-FOV/2 FOV/2]);zlim([-FOV/2 FOV/2]);
xlabel('Needle_Z');ylabel('Needle_Y');zlabel('Needle_X'); 
title(' Coil in the Initial Needle Coordinate System');


%%% Rotate to required radial angle. (Except for when the loop is angled at
%%% the tip, Rotating such a loop does not work very well).
if  Ang_Loop_Flag ~= 1    
RotN = [ 1 0 0; 0 cosd(theta_N)  sind(theta_N) ;  0 -sind(theta_N) cosd(theta_N) ];                       
L = L * RotN;
end


%%%%% Chop off the coil if tip clearance is less than bevel zone
if tip_clr(1) < 6 && y_cutoff_flag == 1
    
   id = find(L(:,2) < y_cutoff);
   L = L(id,:);
   
%   %%% Bring Leads to the center
%      if L(1,2) ~= 0          
%          %%% For Closed Needle
%     %    Lead_X = min(needle_bottom_edge-3.85,needle_bottom_edge-tip_clr(1));
% 
%          %%% For Open Needle
%          clr  = needle_od/2/tand(bevel_ang);
%          Lead_X = min(needle_bottom_edge-clr,needle_bottom_edge-tip_clr(1));
% 
%          L = [ [Lead_X,0,rad] ;L; [Lead_X,0,-rad]]; 
%      end             
end
   

%%%% Add Leads if needed
if leads == 1        
      %%% Full Loop with Leads          
      L = [ [needle_top_edge,L(1,2),L(1,3)] ;L; [needle_top_edge,L(end,2),L(end,3)]];                     
end


%  L =[  needle_top_edge  rad 0;
%               (needle_bottom_edge-tip_clr)* ones(10,1)-x_step'  y_step'  zeros(10,1) ;  
%                needle_top_edge  -rad  0  ];
           
           
          

subplot(1,3,2); plot3(L(:,3),L(:,2),L(:,1),'r','linewidth',3); grid on;
xlim([-FOV/2-10 FOV/2+10]);ylim([-FOV/2 FOV/2]);zlim([-FOV/2 FOV/2]);
xlabel('Needle_Z');ylabel('Needle_Y');zlabel('Needle_X'); 
title(' Coil in the Final Needle Coordinate System');



%%%%%%%%%% TRANSFORM COIL GEOMETRY TO SCANNER GEOMETRY  %%%%%%%%%%%%%%%
%%%% Transform Coil Geometry to Scanner Geometry so that Coil Field in 
%%%% the 'Z' direction (W) may be selected directly after BiotSavarts.

%%% Transform to Unrotated Magnet coordinate system
Needle_Mag_Transform = [ -1 0 0; 0 1 0; 0 0 1];
L_Mag = L * Needle_Mag_Transform;

 
%%% Transform to Rotated Magnet coordinate system
% Mag_MagRot_Transform = [cosd(Ang_Y)*cosd(Ang_Z), cosd(Ang_X)*sind(Ang_Z)+ sind(Ang_X)*sind(Ang_Y)*cosd(Ang_Z), sind(Ang_X)*sind(Ang_Z)- cosd(Ang_X)*sind(Ang_Y)*cosd(Ang_Z) ;...
%           -cosd(Ang_Y)*sind(Ang_Z), cosd(Ang_X)*cosd(Ang_Z)- sind(Ang_X)*sind(Ang_Y)*sind(Ang_Z), sind(Ang_X)*cosd(Ang_Z)+ cosd(Ang_X)*sind(Ang_Y)*sind(Ang_Z);...    
%                 sind(Ang_Y)  ,                    -sind(Ang_X)*cosd(Ang_Y)  ,                                 cosd(Ang_X)*cosd(Ang_Y)];


% L_Final = L_Mag * Mag_MagRot_Transform;


Mag_MagRot_Transform = [cosd(Ang_Y)*cosd(Ang_Z), -cosd(Ang_X)*sind(Ang_Z)+ sind(Ang_X)*sind(Ang_Y)*cosd(Ang_Z), sind(Ang_X)*sind(Ang_Z)+cosd(Ang_X)*sind(Ang_Y)*cosd(Ang_Z) ;...
           cosd(Ang_Y)*sind(Ang_Z), cosd(Ang_X)*cosd(Ang_Z)+sind(Ang_X)*sind(Ang_Y)*sind(Ang_Z), -sind(Ang_X)*cosd(Ang_Z)+ cosd(Ang_X)*sind(Ang_Y)*sind(Ang_Z);...    
                 -sind(Ang_Y)  ,                    sind(Ang_X)*cosd(Ang_Y)  ,                                 cosd(Ang_X)*cosd(Ang_Y)];

L_Final =  Mag_MagRot_Transform * L_Mag';
L_Final = L_Final';


subplot(1,3,3); plot3(L_Final(:,3),L_Final(:,2),L_Final(:,1),'r','linewidth',4); grid on;
xlim([-FOV/2-10 FOV/2+10]);ylim([-FOV/2 FOV/2]);zlim([-FOV/2 FOV/2]);
xlabel('Magnet_Z');ylabel('Magnet_Y');zlabel('Magnet_X'); 
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
title(' Coil in the Magnet Coordinate System'); 


figure;
ds = 32;
scatter3(Needle_Coord_MagRot(1:ds:end,3)/10,Needle_Coord_MagRot(1:ds:end,2)/10,Needle_Coord_MagRot(1:ds:end,1)/10,'.','MarkerEdgeColor',[.9 .9 .9],...
    'MarkerFaceColor',[.9 .9 .9]);
xlabel('Magnet_Z');ylabel('Magnet_Y');zlabel('Magnet_X'); 
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
title(' Needle and Coil in the Magnet Coordinate System'); hold on
plot3(L_Final(:,3),L_Final(:,2),L_Final(:,1),'r','linewidth',10); grid on;
xlim([-FOV/2 FOV/2]);ylim([-FOV/2 FOV/2]);zlim([-FOV/2 FOV/2]);


%% Numerical integration of Biot-Savart law
clear W dW
clear x_range y_range z_range Ind_X Ind_Y Ind_Z
tic 

x = linspace(-FOV/2, FOV/2, ND);
y = linspace(-FOV/2, FOV/2, ND);
z = linspace(-FOV/2, FOV/2, ND); 
[Y, X, Z] = meshgrid(y,-x,-z);

x_range = -ND/2:ND/2-1;
y_range = -ND/2:ND/2-1;
z_range = -ND/2:ND/2-1;

[ Ind_Y, Ind_X ,Ind_Z] = meshgrid(y_range,-x_range,z_range); 

Ind_X = single(Ind_X); 
Ind_Y = single(Ind_Y); 
Ind_Z = single(Ind_Z);

Ind_X_vec = Ind_X(:); 
Ind_Y_vec = Ind_Y(:); 
Ind_Z_vec = Ind_Z(:);

% Induction vector components B = (U, V, W);
% U = single(zeros(ND, ND, ND));
% V = single(zeros(ND, ND, ND));
W = single(zeros(ND, ND, ND));

X = single(X);
Y = single(Y);
Z = single(Z);

% The curve is discretized in Nl points, we iterate on the Nl-1
% segments. Each segment is discretized with a "ds" length step
% to evaluate a "dB" increment of the induction "B".

Lx = [];
Ly = [];
Lz = [];
Nl = numel(L)/3; % Number of points of the curve
Npi = zeros(Nl-1,1);


for pCurv = 1:Nl-1
% for pCurv = 1
    % Length of the curve element
    len = norm(L_Final(pCurv,:) - L_Final(pCurv+1,:));
    % Number of points for the curve-element discretization                
    Npi(pCurv) = ceil(len/Res_CurElem);

    if Npi(pCurv) < 3
%      close(Wait);
%       error('Integration step is too big, Reduce Res_CurElem !!')
        Npi(pCurv) = ceil(len/(Res_CurElem/2));
    end
    
    % Curve-element discretization
    Lx = cat(2,Lx,linspace(L_Final(pCurv,1), L_Final(pCurv+1,1), Npi(pCurv)));
    Ly = cat(2,Ly,linspace(L_Final(pCurv,2), L_Final(pCurv+1,2), Npi(pCurv)));
    Lz = cat(2,Lz,linspace(L_Final(pCurv,3), L_Final(pCurv+1,3), Npi(pCurv)));

end

dLx_v = diff(Lx);
dLy_v = diff(Ly);
dLz_v = diff(Lz);
  
Npi_Total = length(Lx);
cNpi = cumsum(Npi);
Factor = -1*Current*gamma/4/pi;
                
toc                            
for ind = 1 : length(Ind_X_vec)
    
     i = Ind_X_vec(ind) + ND/2+1;
     j = Ind_Y_vec(ind) + ND/2+1;
     k = Ind_Z_vec(ind) + ND/2+1;
         
    if i > ND || j > ND ||  k > ND
        continue
    else
%      waitbar(ind/length(Ind_X_vec), Wait)
       % Ptest is the point of the field where we calculate induction             
       %%%%% X Grid should be up-down, Y Grid should be left-right, Z Grid
       %%%%% should be Foot-Head
         pTest = [X(i,j,k) Y(i,j,k) Z(i,j,k)]; 

         % Integration
            for s = 1:Npi_Total-1
                                                     
                % Vector connecting the infinitesimal curve-element 
                % point and field point "pTest"
%                 
                Rx = Lx(s) - pTest(1);
                Ry = Ly(s) - pTest(2);
                Rz = Lz(s) - pTest(3);
                                            
                % Infinitesimal curve-element components
%                 dLx = Lx(s+1) - Lx(s);
%                 dLy = Ly(s+1) - Ly(s);
%                 dLz = Lz(s+1) - Lz(s);

                % Modules
%                dL = sqrt(dLx^2 + dLy^2 + dLz^2);
                dL = sqrt(dLx_v(s)^2 + dLy_v(s)^2 + dLz_v(s)^2);
                R = sqrt(Rx^2 + Ry^2 + Rz^2);
                                               
                % Biot-Savart
%               dU = -1*Current*gamma/4/pi*(dLy*Rz - dLz*Ry)/R/R/R;
%               dV = -1*Current*gamma/4/pi*(dLz*Rx - dLx*Rz)/R/R/R;
%               dW = -1*Current*gamma/4/pi*(dLx*Ry - dLy*Rx)/R/R/R;
%                  
%                dW = Factor*(dLx*Ry - dLy*Rx)/R/R/R;
                 dW = Factor*(dLx_v(s)*Ry - dLy_v(s)*Rx)/R/R/R;
                 
                 
                 %%%% Split Spear Loop
                 if leads == 1
                      %%% Halfing the field for the split spear coil tip
                      if spearloop_flag == 1
                       if s > cNpi(1) && s <= cNpi(end-1)
                            dW = dW /2;
                       end
                      elseif spearloop_flag == 2
                           dW = dW /2;
                           
                      elseif  speardoublewireloop_flag == 1  %%%% Split Spear Double Loop
                          dW = dW /2;
                      end
                 else
                       %%% Halfing the field for the split spear coil tip
                      if spearloop_flag == 1 || spearloop_flag == 2             
                            dW = dW /2;      
                      end
                 end
                 
                 
                 
                      
                % Add increment to the main field
%                U(i,j,k) = U(i,j,k) + dU;
%                V(i,j,k) = V(i,j,k) + dV;
                W(i,j,k) = W(i,j,k) + dW;
                
            end              
    end
           
end

%close(Wait);

%%% Convert to Hz : 1T = 42.57*1e6 Hz for proton
% U = U * 42.57*1e6;  
% V = V * 42.57*1e6;  
W = W * 42.57*1e6;  

toc


   