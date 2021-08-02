%%
clear;clc;
% -------------------------------------------------------------------------
% user input
% -------------------------------------------------------------------------
%load horizons
[~,top_Sele] = mhdrload('Top_sele');
[~,top_Hod] = mhdrload('Top_hod');
[~,top_Kim] = mhdrload('Top_kim');
[~,Zone3] = mhdrload('zone3');
[~,Zone6] = mhdrload('zone6');
[~,Zone9] = mhdrload('zone9');
[~,top_Pentland] = mhdrload('Top_pentland');
[~,top_Triassic] = mhdrload('Top_triassic');
[~,top_Salt] = mhdrload('Top_salt');
% -------------------------------------------------------------------------
% load strains (top to bottom)
%layer 1 Top Sele
% [~,exx1] = mhdrload('Avg_Exx_OB');
% [~,eyy1] = mhdrload('Avg_Eyy_OB');
[~,ezz1] = mhdrload('Ezz_Top_sele');

%layer 2 Top Hod
% [~,exx2] = mhdrload('Avg_Exx_Res');
% [~,eyy2] = mhdrload('Avg_Eyy_Res');
[~,ezz2] = mhdrload('Ezz_Top_Hod');

%layer 3 Top Kim
% [~,exx3] = mhdrload('Avg_Exx_UB');
% [~,eyy3] = mhdrload('Avg_Eyy_UB');
[~,ezz3] = mhdrload('Ezz_Top_Kim');

%layer 4 Zone3
% [~,exx4] = mhdrload('Avg_Exx_UB');
% [~,eyy4] = mhdrload('Avg_Eyy_UB');
[~,ezz4] = mhdrload('Ezz_Zone3');

%layer 5 Zone6
% [~,exx5] = mhdrload('Avg_Exx_UB');
% [~,eyy5] = mhdrload('Avg_Eyy_UB');
[~,ezz5] = mhdrload('Ezz_Zone6');

%layer 6 Zone9
% [~,exx6] = mhdrload('Avg_Exx_UB');
% [~,eyy6] = mhdrload('Avg_Eyy_UB');
[~,ezz6] = mhdrload('Ezz_Zone9');

%layer 7 Top Kim
% [~,exx7] = mhdrload('Avg_Exx_UB');
% [~,eyy7] = mhdrload('Avg_Eyy_UB');
[~,ezz7] = mhdrload('Ezz_Top_Pentland');

%layer 8 Top Kim
% [~,exx8] = mhdrload('Avg_Exx_UB');
% [~,eyy8] = mhdrload('Avg_Eyy_UB');
[~,ezz8] = mhdrload('Ezz_Top_Triassic');

%layer 9 Top Salt
% [~,exx9] = mhdrload('Avg_Exx_UB');
% [~,eyy9] = mhdrload('Avg_Eyy_UB');
[~,ezz9] = mhdrload('Ezz_Top_Salt');


% load intersection
coord_A = [616940 6316040];
coord_B = [616940 6323040];

lat_step  = 50; %m
vert_step = 15; %ft
td = 32000;
%% ------------------------------------------------------------------------
% extract 2D section from the maps
% -------------------------------------------------------------------------
dist_section = sqrt(sum(coord_A-coord_B).^2);
dx = lat_step * (coord_B(1)-coord_A(1)) / sqrt(sum(coord_A-coord_B).^2);
dy = lat_step * (coord_B(2)-coord_A(2)) / sqrt(sum(coord_A-coord_B).^2);

null_2d = meshgrid(1:length(0:vert_step:td),1:length( coord_A(2):dy:coord_B(2) ))';

if dx==0
    section_grid = [coord_A(1)*ones(1,length( coord_A(2):dy:coord_B(2) )); coord_A(2):dy:coord_B(2)]';
elseif dy==0
    section_grid = [coord_A(1):dx:coord_B(1); coord_A(2)*ones(1,length(coord_A(1):dx:coord_B(1)))]';
else
    section_grid = [coord_A(1):dx:coord_B(1);coord_A(2):dy:coord_B(2)]';
end

%
top_Sele_ind = repmat(ceil(griddata(top_Sele(:,1),top_Sele(:,2),top_Sele(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
top_Hod_depth_ind = repmat(ceil(griddata(top_Hod(:,1),top_Hod(:,2),top_Hod(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
top_Kim_ind = repmat(ceil(griddata(top_Kim(:,1),top_Kim(:,2),top_Kim(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
Zone3_ind = repmat(ceil(griddata(Zone3(:,1),Zone3(:,2),Zone3(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
Zone6_ind = repmat(ceil(griddata(Zone6(:,1),Zone6(:,2),Zone6(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
Zone9_ind = repmat(ceil(griddata(Zone9(:,1),Zone9(:,2),Zone9(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
top_Pentland_ind = repmat(ceil(griddata(top_Pentland(:,1),top_Pentland(:,2),top_Pentland(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
top_Triassic_ind = repmat(ceil(griddata(top_Triassic(:,1),top_Triassic(:,2),top_Triassic(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
top_Salt_ind = repmat(ceil(griddata(top_Salt(:,1),top_Salt(:,2),top_Salt(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);

% layer props
%exx1_sec = repmat(griddata(exx1(:,1),exx1(:,2),exx1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy1_sec = repmat(griddata(eyy1(:,1),eyy1(:,2),eyy1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz1_sec = repmat(griddata(ezz1(:,1),ezz1(:,2),ezz1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx2_sec = repmat(griddata(exx2(:,1),exx1(:,2),exx2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy2_sec = repmat(griddata(eyy2(:,1),eyy1(:,2),eyy2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz2_sec = repmat(griddata(ezz2(:,1),ezz2(:,2),ezz2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz3_sec = repmat(griddata(ezz3(:,1),ezz3(:,2),ezz3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz4_sec = repmat(griddata(ezz4(:,1),ezz4(:,2),ezz4(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz5_sec = repmat(griddata(ezz5(:,1),ezz5(:,2),ezz5(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz6_sec = repmat(griddata(ezz6(:,1),ezz6(:,2),ezz6(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz7_sec = repmat(griddata(ezz7(:,1),ezz7(:,2),ezz7(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz8_sec = repmat(griddata(ezz8(:,1),ezz8(:,2),ezz8(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

%exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
%eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz9_sec = repmat(griddata(ezz9(:,1),ezz9(:,2),ezz9(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';


%% ------------------------------------------------------------------------
% make 2D gridz
% -------------------------------------------------------------------------
% figure
% exx_2d = null_2d*0;
% exx_2d(null_2d<top_depth_ind) = exx1_sec(null_2d<top_depth_ind) ; imagesc(exx_2d)
% exx_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = exx2_sec(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) ; imagesc(exx_2d)
% exx_2d(null_2d>bottom_depth_ind) = exx3_sec(null_2d>bottom_depth_ind) ; imagesc(exx_2d)

% figure
% eyy_2d = null_2d*0;
% eyy_2d(null_2d<top_depth_ind) = eyy1_sec(null_2d<top_depth_ind) ; imagesc(eyy_2d)
% eyy_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = eyy2_sec(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) ; imagesc(eyy_2d)
% eyy_2d(null_2d>bottom_depth_ind) = eyy3_sec(null_2d>bottom_depth_ind) ; imagesc(eyy_2d)
%% ------------------------------------------------------------------------
figure
ezz_2d = null_2d*0;
ezz_2d(null_2d<top_Sele_ind) = ezz1_sec(null_2d<top_Sele_ind) ; 
ezz_2d(null_2d>=top_Sele_ind & null_2d<=top_Hod_depth_ind) = ezz2_sec(null_2d>=top_Sele_ind & null_2d<=top_Hod_depth_ind) ;
ezz_2d(null_2d>=top_Hod_depth_ind & null_2d<=top_Kim_ind) = ezz3_sec(null_2d>=top_Hod_depth_ind & null_2d<=top_Kim_ind) ;
ezz_2d(null_2d>=top_Kim_ind & null_2d<=Zone3_ind) = ezz4_sec(null_2d>=top_Kim_ind & null_2d<=Zone3_ind) ;
ezz_2d(null_2d>=Zone3_ind & null_2d<=Zone6_ind) = ezz5_sec(null_2d>=Zone3_ind & null_2d<=Zone6_ind) ;
ezz_2d(null_2d>=Zone6_ind & null_2d<=Zone9_ind) = ezz6_sec(null_2d>=Zone6_ind & null_2d<=Zone9_ind) ;
ezz_2d(null_2d>=Zone9_ind & null_2d<=top_Pentland_ind) = ezz7_sec(null_2d>=Zone9_ind & null_2d<=top_Pentland_ind) ;
ezz_2d(null_2d>=top_Pentland_ind & null_2d<=top_Triassic_ind) = ezz8_sec(null_2d>=top_Pentland_ind & null_2d<=top_Triassic_ind) ;
ezz_2d(null_2d>=top_Triassic_ind & null_2d<=top_Salt_ind) = ezz9_sec(null_2d>=top_Triassic_ind & null_2d<=top_Salt_ind) ;imagesc(ezz_2d)
%ezz_2d(null_2d>bottom_depth_ind) = ezz3_sec(null_2d>bottom_depth_ind) ;imagesc(ezz_2d)
velocity_2d = null_2d*0;
velocity_2d(null_2d<Zone3_ind) = 3133;
velocity_2d(null_2d>=Zone3_ind & null_2d<=Zone9_ind) = 2770;
velocity_2d(null_2d>Zone9_ind) = 3133;
imagesc(velocity_2d);