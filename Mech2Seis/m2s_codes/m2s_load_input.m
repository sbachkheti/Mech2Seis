%%
clear;clc;
% -------------------------------------------------------------------------
% user input
% -------------------------------------------------------------------------
%load horizons
[~,top_Sele] = mhdrload('Top_Sele');
[~,top_Hod] = mhdrload('Top_Hod');
[~,top_Kim] = mhdrload('Top_Kim');
[~,Zone3] = mhdrload('Zone3');
[~,Zone6] = mhdrload('Zone6');
[~,Zone9] = mhdrload('Zone9');
[~,top_Pentland] = mhdrload('Top_Pentland');
[~,top_Triassic] = mhdrload('Top_Triassic');
[~,top_Salt] = mhdrload('Top_salt');

% -------------------------------------------------------------------------
% load strains (top to nbottom)
%layer 1
[~,exx1] = mhdrload('Avg_Exx_OB');
[~,eyy1] = mhdrload('Avg_Eyy_OB');
[~,ezz1] = mhdrload('Avg_Ezz_OB');

%layer 2
[~,exx2] = mhdrload('Avg_Exx_Res');
[~,eyy2] = mhdrload('Avg_Eyy_Res');
[~,ezz2] = mhdrload('Avg_Ezz_Res');

%layer 3
[~,exx3] = mhdrload('Avg_Exx_UB');
[~,eyy3] = mhdrload('Avg_Eyy_UB');
[~,ezz3] = mhdrload('Avg_Ezz_UB');

% -------------------------------------------------------------------------
% load intersection
coord_A = [616940 6316040];
coord_B = [616940 6323040];

% coord_A = [616326 6316027];
% coord_B = [608534 6323477];

%
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
top_depth_ind = repmat(ceil(griddata(top_depth(:,1),top_depth(:,2),top_depth(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);
bottom_depth_ind = repmat(ceil(griddata(bottom_depth(:,1),bottom_depth(:,2),bottom_depth(:,3),section_grid(:,1),section_grid(:,2))/vert_step)',[size(null_2d,1),1]);

% layer props
exx1_sec = repmat(griddata(exx1(:,1),exx1(:,2),exx1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
eyy1_sec = repmat(griddata(eyy1(:,1),eyy1(:,2),eyy1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz1_sec = repmat(griddata(ezz1(:,1),ezz1(:,2),ezz1(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

exx2_sec = repmat(griddata(exx2(:,1),exx1(:,2),exx2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
eyy2_sec = repmat(griddata(eyy2(:,1),eyy1(:,2),eyy2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz2_sec = repmat(griddata(ezz2(:,1),ezz1(:,2),ezz2(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

exx3_sec = repmat(griddata(exx3(:,1),exx1(:,2),exx3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
eyy3_sec = repmat(griddata(eyy3(:,1),eyy1(:,2),eyy3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';
ezz3_sec = repmat(griddata(ezz3(:,1),ezz1(:,2),ezz3(:,3),section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';



%% ------------------------------------------------------------------------
% make 2D gridz
% -------------------------------------------------------------------------
figure
exx_2d = null_2d*0;
exx_2d(null_2d<top_depth_ind) = exx1_sec(null_2d<top_depth_ind) ; imagesc(exx_2d)
exx_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = exx2_sec(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) ; imagesc(exx_2d)
exx_2d(null_2d>bottom_depth_ind) = exx3_sec(null_2d>bottom_depth_ind) ; imagesc(exx_2d)

figure
eyy_2d = null_2d*0;
eyy_2d(null_2d<top_depth_ind) = eyy1_sec(null_2d<top_depth_ind) ; imagesc(eyy_2d)
eyy_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = eyy2_sec(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) ; imagesc(eyy_2d)
eyy_2d(null_2d>bottom_depth_ind) = eyy3_sec(null_2d>bottom_depth_ind) ; imagesc(eyy_2d)
%% ------------------------------------------------------------------------
figure
ezz_2d = null_2d*0;
ezz_2d(null_2d<top_depth_ind) = ezz1_sec(null_2d<top_depth_ind) ; imagesc(ezz_2d)
ezz_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = ezz2_sec(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) ;imagesc(ezz_2d)
ezz_2d(null_2d>bottom_depth_ind) = ezz3_sec(null_2d>bottom_depth_ind) ;imagesc(ezz_2d)

R1 = 5;
R2 = 2 ;
R3 = 5;
R_2d = null_2d*0;
R_2d(null_2d<top_depth_ind) = R1 ; imagesc(R_2d)
R_2d(null_2d>=top_depth_ind & null_2d<=bottom_depth_ind) = R2; imagesc(R_2d)
R_2d(null_2d>bottom_depth_ind) = R3 ; imagesc(R_2d)

% plot(top_depth_ind,'DisplayName','bottom_depth_sec');hold on;plot(bottom_depth_ind,'DisplayName','top_depth_sec');hold off;
% plot the line
% plot(section_grid(:,1),section_grid(:,2),'.')
%% ------------------------------------------------------------------------
% rock phsyics
% -------------------------------------------------------------------------
R_ob = 5;
R1_res = 2;
Rp = 0.5*(R_ob*ezz1_sec - (1+ R1_res)*ezz2_sec);
Rp_sec = repmat(griddata(0,Rp,0,section_grid(:,1),section_grid(:,2)),[1,size(null_2d,1),1])';

% Rp0 = 0.5*(R_2d(1:end-1,:).*ezz_2d(1:end-1,:) - (R_2d(2:end,:)).*ezz_2d(2:end,:));
% figure;plot(Rp0(:,140),[1:2133])
%%-------------------------------------------------------------------------
% Plot the reflectivity
figure
Rp_2d = null_2d*0;
Rp_2d(null_2d<top_depth_ind) = 0 ;
Rp_2d(null_2d==top_depth_ind) = Rp(null_2d==top_depth_ind);
Rp_2d(null_2d>bottom_depth_ind) = 0;






