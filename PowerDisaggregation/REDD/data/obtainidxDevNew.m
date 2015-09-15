idxDevNew=cell(6,19);  % index 1 runs over houses; index 2 runs over devices
% Convention for index 2:
% 1-lightning
% 2-refrigerator
% 3-dishwasher
% 4-microwave
% 5-oven
% 6-kitchen_outlets
% 7-washer_dryer
% 8-bathroom_gfi
% 9-stove
% 10-electric_heat
% 11-disposal
% 12-outlets_unknown
% 13-electronics
% 14-furance
% 15-smoke_alarms
% 16-air_conditioning
% 17-miscellaeneous
% 18-subpanel
% 19-outdoor_outlets

%% H1
% 1 mains
% 2 mains
% 3 oven
% 4 oven
% 5 refrigerator
% 6 dishwaser
% 7 kitchen_outlets
% 8 kitchen_outlets
% 9 lighting
% 10 washer_dryer
% 11 microwave
% 12 bathroom_gfi
% 13 electric_heat
% 14 stove
% 15 kitchen_outlets
% 16 kitchen_outlets
% 17 lighting
% 18 lighting
% 19 washer_dryer
% 20 washer_dryer
idxDevNew{1,1}=[9,17,18];
idxDevNew{1,2}=5;
idxDevNew{1,3}=6;
idxDevNew{1,4}=11;
idxDevNew{1,5}=[3,4];
idxDevNew{1,6}=[7,8,15,16];
idxDevNew{1,7}=[10,19,20];
idxDevNew{1,8}=12;
idxDevNew{1,9}=14;
idxDevNew{1,10}=13;
idxDevNew{1,11}=[];
idxDevNew{1,12}=[];
idxDevNew{1,13}=[];
idxDevNew{1,14}=[];
idxDevNew{1,15}=[];
idxDevNew{1,16}=[];
idxDevNew{1,17}=[];
idxDevNew{1,18}=[];
idxDevNew{1,19}=[];

%% H2
% 1 mains
% 2 mains
% 3 kitchen_outlets
% 4 lighting
% 5 stove
% 6 microwave
% 7 washer_dryer
% 8 kitchen_outlets
% 9 refrigerator
% 10 dishwaser
% 11 disposal
idxDevNew{2,1}=4;
idxDevNew{2,2}=9;
idxDevNew{2,3}=10;
idxDevNew{2,4}=6;
idxDevNew{2,5}=[];
idxDevNew{2,6}=[3,8];
idxDevNew{2,7}=7;
idxDevNew{2,8}=[];
idxDevNew{2,9}=5;
idxDevNew{2,10}=[];
idxDevNew{2,11}=11;
idxDevNew{2,12}=[];
idxDevNew{2,13}=[];
idxDevNew{2,14}=[];
idxDevNew{2,15}=[];
idxDevNew{2,16}=[];
idxDevNew{2,17}=[];
idxDevNew{2,18}=[];
idxDevNew{2,19}=[];

%% H3
% 1 mains
% 2 mains
% 3 outlets_unknown
% 4 outlets_unknown
% 5 lighting
% 6 electronics
% 7 refrigerator
% 8 disposal
% 9 dishwaser
% 10 furance
% 11 lighting
% 12 outlets_unknown
% 13 washer_dryer
% 14 washer_dryer
% 15 lighting
% 16 microwave
% 17 lighting
% 18 smoke_alarms
% 19 lighting
% 20 bathroom_gfi
% 21 kitchen_outlets
% 22 kitchen_outlets
idxDevNew{3,1}=[5,11,15,17,19];
idxDevNew{3,2}=7;
idxDevNew{3,3}=9;
idxDevNew{3,4}=16;
idxDevNew{3,5}=[];
idxDevNew{3,6}=[21 22];
idxDevNew{3,7}=[13 14];
idxDevNew{3,8}=20;
idxDevNew{3,9}=[];
idxDevNew{3,10}=[];
idxDevNew{3,11}=8;
idxDevNew{3,12}=[3,4,12];
idxDevNew{3,13}=6;
idxDevNew{3,14}=10;
idxDevNew{3,15}=18;
idxDevNew{3,16}=[];
idxDevNew{3,17}=[];
idxDevNew{3,18}=[];
idxDevNew{3,19}=[];

%% H4
% 1 mains
% 2 mains
% 3 lighting
% 4 furance
% 5 kitchen_outlets
% 6 outlets_unknown
% 7 washer_dryer
% 8 stove
% 9 air_conditioning
% 10 air_conditioning
% 11 miscellaeneous
% 12 smoke_alarms
% 13 lighting
% 14 kitchen_outlets
% 15 dishwaser
% 16 bathroom_gfi
% 17 bathroom_gfi
% 18 lighting
% 19 lighting
% 20 air_conditioning
idxDevNew{4,1}=[3,13,18,19];
idxDevNew{4,2}=[];
idxDevNew{4,3}=15;
idxDevNew{4,4}=[];
idxDevNew{4,5}=[];
idxDevNew{4,6}=[5,14];
idxDevNew{4,7}=7;
idxDevNew{4,8}=[16,17];
idxDevNew{4,9}=8;
idxDevNew{4,10}=[];
idxDevNew{4,11}=[];
idxDevNew{4,12}=6;
idxDevNew{4,13}=[];
idxDevNew{4,14}=4;
idxDevNew{4,15}=12;
idxDevNew{4,16}=[9 10 20];
idxDevNew{4,17}=11;
idxDevNew{4,18}=[];
idxDevNew{4,19}=[];

%% H5
% 1 mains
% 2 mains
% 3 microwave
% 4 lighting
% 5 outlets_unknown
% 6 furance
% 7 outlets_unknown
% 8 washer_dryer
% 9 washer_dryer
% 10 subpanel
% 11 subpanel
% 12 electric_heat
% 13 electric_heat
% 14 lighting
% 15 outlets_unknown
% 16 bathroom_gfi
% 17 lighting
% 18 refrigerator
% 19 lighting
% 20 dishwaser
% 21 disposal
% 22 electronics
% 23 lighting
% 24 kitchen_outlets
% 25 kitchen_outlets
% 26 outdoor_outlets
idxDevNew{5,1}=[4,14,17,19,23];
idxDevNew{5,2}=18;
idxDevNew{5,3}=20;
idxDevNew{5,4}=3;
idxDevNew{5,5}=[];
idxDevNew{5,6}=[24,25];
idxDevNew{5,7}=[8,9];
idxDevNew{5,8}=16;
idxDevNew{5,9}=[];
idxDevNew{5,10}=[12,13];
idxDevNew{5,11}=21;
idxDevNew{5,12}=[5,7,15];
idxDevNew{5,13}=22;
idxDevNew{5,14}=6;
idxDevNew{5,15}=[];
idxDevNew{5,16}=[];
idxDevNew{5,17}=[];
idxDevNew{5,18}=[10,11];
idxDevNew{5,19}=26;

%% H6
% 1 mains
% 2 mains
% 3 kitchen_outlets
% 4 washer_dryer
% 5 stove
% 6 electronics
% 7 bathroom_gfi
% 8 refrigerator
% 9 dishwaser
% 10 outlets_unknown
% 11 outlets_unknown
% 12 electric_heat
% 13 kitchen_outlets
% 14 lighting
% 15 air_conditioning
% 16 air_conditioning
% 17 air_conditioning
idxDevNew{6,1}=14;
idxDevNew{6,2}=8;
idxDevNew{6,3}=9;
idxDevNew{6,4}=[];
idxDevNew{6,5}=[];
idxDevNew{6,6}=[3,13];
idxDevNew{6,7}=4;
idxDevNew{6,8}=7;
idxDevNew{6,9}=5;
idxDevNew{6,10}=12;
idxDevNew{6,11}=[];
idxDevNew{6,12}=[10,11];
idxDevNew{6,13}=6;
idxDevNew{6,14}=[];
idxDevNew{6,15}=[];
idxDevNew{6,16}=[15,16,17];
idxDevNew{6,17}=[];
idxDevNew{6,18}=[];
idxDevNew{6,19}=[];

%% Sanity check
nH = [20 11 22 20 26 17];
for h=1:6
    aux = [];
    for i=1:19
        aux = [aux idxDevNew{h,i}];
    end
    if(length(unique(aux))~=nH(h)-2)
        disp(['Warning: something is wrong with House ' num2str(h)]);
        disp(['   Missing indexes: ' num2str(setdiff(1:nH(h),unique(aux)))]);
    end
end

%% Save to file
save('idxDevNew.mat','idxDevNew');
