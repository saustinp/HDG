fprintf('  --> Initializing HDG ... v1\n')
d0=fileparts([pwd,filesep]);
addpath([d0,'/util']);
addpath([d0,'/mesh/mkmesh']);
addpath([d0,'/mesh/cmesh']);
addpath([d0,'/mesh/foilmesh']);
addpath([d0,'/mesh/airfoilTools']);
addpath([d0,'/mesh/airfoilTools/geometries']);
addpath([d0,'/mesh']);
addpath([d0,'/kernel']);
addpath([d0,'/master']);
addpath([d0,'/plot']);
addpath([d0,'/shockcapturing']);
addpath([d0,'/problem/ionicwind2/mesh']);
addpath([d0,'/problem/ionicwind2/swarm']);
fprintf(' \n Done. \n\n');
clear d0
