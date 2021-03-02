% clear all
% 
% DELTA   = 0.08;  %depreciation rate
% ALFA    = 0.32;  %capital share
% BETTA   = 0.98; %discount rate
% G       = 1.014;
% SIGM    = 0.9;
% % LAMBDAP = 0.95;
% LAMBDAZ = 0.92;
% sigma_Z = 0.017;
% sigma_P = 0.03;
% FRISCHELAS = 0.5;
% STEADYSTATEL = 0.3;
% MU = 0.086;
% T=200;
% shock = "none";
% k0_mult=1;
% startopposite = 0;
% regimechanges = 0;
% regime_change_frequency = 50;
% order = 4;
% randomseq = 2;
% 
% MultiplicativeU = 1;
% 
% addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
% import casadi.*
% 
% if MultiplicativeU
%     utility_function_str = "Multiplicative";
% else
%     utility_function_str = "Additive";
% end
% 
% LAMBDAPs = [0.8,0.8,0.95,0.95];
% useopposites = [0,1,0,1];
% 
% 
% 
%     
% N = length(LAMBDAPs);
% output_vars=9;
% output = NaN([T output_vars N]);
% if shock == "historical"
%     output = NaN([33 output_vars N]);
% end
% tic
% parfor i = 1:N
%     LAMBDAP = LAMBDAPs(i);
%     useopposite = useopposites(i);
%     [output(:,:,i)] = model2_testing(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,regimechanges,regime_change_frequency,randomseq,order,useopposite);
% end
% toc
% filename = join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/Tests.xlsx"],"");
% sheetloc = 'A5';
% writematrix(output,filename,'Sheet',1,'Range',sheetloc)
% % writematrix(output(:,4,:),filename,'Sheet',1,'Range','A4')
% % writematrix(output(:,8,:),filename,'Sheet',2,'Range','A4')
% % writematrix(output,filename,'Sheet',5,'Range','A5')
% % writematrix(output,filename,'Sheet',6,'Range','B4')


clear all

DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
% LAMBDAP = 0.95;
LAMBDAZ = 0.92;
sigma_Z = 0.017;
sigma_P = 0.03;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU = 0.086;
T=110;
k0_mult=1;
startopposite = 0;
regimechanges = 0;
regime_change_frequency = 50;
order = 4;
randomseq = 2;
useopposite = 0;
MultiplicativeU = 1;
% sigma_P = 0.0433;
shiftstart = "kPZ";

addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
import casadi.*

nstartpoints = 10;
LAMBDAPs = repmat([0.8,0.8,0.95,0.95],1,nstartpoints); %,0.8,0.95,0.8,0.95,0.8,0.95,0.8,0.95,0.8,0.95];
shocks = repmat(["irf_P","irf_none","irf_P","irf_none"],1,nstartpoints);
% shocks = ["irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z","irf_Z",...
%     "irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P","irf_P"]; %,"irf_P","irf_P","irf_K","irf_K"]; %,"irf_P","irf_P","irf_K","irf_K","irf_Z","irf_Z","irf_P","irf_P","irf_K","irf_K"];
% shiftstarts = ["k","k","k","k","Z","Z","Z","Z","P","P","P","P",...
%     "k","k","k","k","Z","Z","Z","Z","P","P","P","P"];
% shiftdirs = ["high","high","low","low","high","high","low","low","high","high","low","low",...
%     "high","high","low","low","high","high","low","low","high","high","low","low"];
startpoints = repelem(1:nstartpoints,1,4);
    
N = length(shocks);
output_vars = 9;
hardcode_irf_T = 100;
output = NaN([hardcode_irf_T output_vars N]);
tic
parfor i = 1:N
    LAMBDAP = LAMBDAPs(i);
    shock = shocks(i);
    startpoint = startpoints(i);
    [output(:,:,i)] = model2_testing(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,...
        sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,...
        startopposite,regimechanges,regime_change_frequency,randomseq,order,useopposite,shiftstart,startpoint);
end
toc
filename = join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/IRFsAwayFromSteadyState.xlsx"],"");
sheetloc = 'A4';
writematrix(output(:,[1,4,8,9],1:2:20)-output(:,[1,4,8,9],2:2:20),filename,'Sheet',1,'Range',sheetloc)
writematrix(output(:,[1,4,8,9],21:2:40)-output(:,[1,4,8,9],22:2:40),filename,'Sheet',2,'Range',sheetloc)
