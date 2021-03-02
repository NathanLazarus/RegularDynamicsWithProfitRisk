function [rgwkcl_mat] = model2_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,LAMBDAphi,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,regimechanges,regime_change_frequency,randomseq,order)

% clear all


hardcode_irf_T = 100;

use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients = true;

defaults = {0.08,0.32,0.98,1.014,...
    0.9,0.95,0.92,0.74,0.017,...
    0.03,0.2,0.086,0.5,0.3,...
    200,"historical",1,0,...
    1,0,50,2,4};

var = ["DELTA","ALFA","BETTA","G",...
    "SIGM","LAMBDAP","LAMBDAZ","LAMBDAphi","sigma_Z",...
    "sigma_P","sigma_phi","MU","FRISCHELAS","STEADYSTATEL",...
    "T","shock","k0_mult","MultiplicativeU",...
    "startopposite","regimechanges","regime_change_frequency","randomseq","order"];



for i = 1:length(defaults)
    if ~exist(var{i},"var")
        if class(defaults{i}) == "string"
            eval(sprintf('%s = "%s";',var(i),defaults{i}))
        else
            eval(sprintf("%s = %g;",var(i),defaults{i}))
        end
    end
end

folder = fileparts(which(mfilename));
addpath(join([folder,'\','MyHelperFunctions'],""))
addpath(join([folder,'\','SchmittGroheUribeHelperFunctions'],""))

if MultiplicativeU
    u = multiplicative_u;
    utility_function_str = "Multiplicative";
else
    u = additive_u;
    utility_function_str = "Additive";
end

max_order_I_can_handle = 4;
nstates_plus_shocks = 10;
if use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients
    SchmittGrohe_decision_func = "SchmittGrohe";
    max_order_I_can_handle = 2;
else
    SchmittGrohe_decision_func = "";
end

decision_func_to_use = str2func(join(["decision_func_","order",max_order_I_can_handle,SchmittGrohe_decision_func],""));

LAMBDAPhigh = 0.95;
LAMBDAPlow = 0.8;
if LAMBDAP == LAMBDAPhigh
    LAMBDAPopposite = LAMBDAPlow;
end
if LAMBDAP == LAMBDAPlow
    LAMBDAPopposite = LAMBDAPhigh;
end

ZSTAR = 1; %steady-state value of technology shock 
PHISTAR = 1;
PSTAR = G^(1/(1-LAMBDAP)); %steady state markup
PSTAR = G^(1/(1-0.8));

sym_labor_supply = laborsupply(u);
intertemporal_euler_ss = dupdcp_over_dudc(u,1);
intertemporal_euler_sym = dupdcp_over_dudc(u,0);

addpath("C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1")
import casadi.*

value_of_P_where_LSTAR_equals_STEADYSTATEL = G^(1/(1-LAMBDAPlow));

[~,~,~,~,~,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,PHISTAR,sym_labor_supply,intertemporal_euler_ss,u);
% GAMA = 28.1677381532213;
[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model2_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTAR,ETA,GAMA,SIGM,ZSTAR,PHISTAR,sym_labor_supply,intertemporal_euler_ss);
stockSTAR = G*((PSTAR-1)/PSTAR)*y_func(KSTAR,LSTAR,ZSTAR,ALFA,PHISTAR)/((1+RSTAR) - G);
% KSTAR = 0.835777225015446;
% LSTAR = 0.273094921033578;
% CSTAR = 0.312063920996547;
% RSTAR = little_r(KSTAR,LSTAR,PSTAR,1,ALFA,DELTA)
fprintf("{%.15g, %.15g, %.15g, %.15g, %.15g, %.15g, %d}\n",KSTAR,CSTAR,LSTAR,stockSTAR,GAMA,ETA,MultiplicativeU)
% [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model2_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,PSTAR,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u)

ssvals = [KSTAR,ZSTAR,ZSTAR,ZSTAR,ZSTAR,PSTAR,PHISTAR];

if startopposite
    PSTARopposite = G^(1/(1-LAMBDAPopposite));
    [KSTARopposite,CSTARopposite,LSTARopposite,WSTARopposite,RSTARopposite]=model2_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTARopposite,ETA,GAMA,SIGM,ZSTAR,PHISTAR,sym_labor_supply,intertemporal_euler_ss);
    stockSTARopposite = G*((PSTARopposite-1)/PSTARopposite)*y_func(KSTARopposite,LSTARopposite,ZSTAR,ALFA,PHISTAR)/((1+RSTARopposite) - G);
    k0_mult = KSTARopposite/KSTAR;
    ssvalsopposite = [KSTARopposite,ZSTAR,ZSTAR,ZSTAR,ZSTAR,PSTARopposite,PHISTAR];
end

if string(shock) == "historical_postwar_trend"
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks_relative_to_postwar_trend.csv";
    shock = "historical";
elseif string(shock) == "historical_endogenous_P_postwar_trend"
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks_relative_to_postwar_trend.csv";
    shock = "historical_endogenous_P";
else
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks.csv";
end

if use_SchmittGrohe_Uribe_Matlab_code_to_find_coefficients

    get_coefs_SchmittGrohe_Uribe
    
else

    for order_to_check_for_coefs=1:max_order_I_can_handle
        if order >= order_to_check_for_coefs
            eval(sprintf('order%1$d = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/MathematicaCoefs/",utility_function_str,"UtilityLambdaP",string(round(LAMBDAP*100)),"coefs",%1$d,".csv"],""));',order_to_check_for_coefs))
            eval(sprintf('order%1$d_opposite = readmatrix(join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model2/MathematicaCoefs/",utility_function_str,"UtilityLambdaP",string(round(LAMBDAPopposite*100)),"coefs%1$d.csv"],""));',order_to_check_for_coefs))
        else
            eval(sprintf('order%1$d = zeros([4 nstates_plus_shocks^%1$d]);',order_to_check_for_coefs))
            eval(sprintf('order%1$d_opposite = zeros([4 nstates_plus_shocks^%1$d]);',order_to_check_for_coefs))
        end
    end

    dec_k=[KSTAR,order1(1,:),1/2*order2(1,:),1/6*order3(1,:),1/24*order4(1,:)];
    dec_stock=[stockSTAR,order1(2,:),1/2*order2(2,:),1/6*order3(2,:),1/24*order4(2,:)];
    dec_l=[LSTAR,order1(3,:),1/2*order2(3,:),1/6*order3(3,:),1/24*order4(3,:)];
    dec_c=[CSTAR,order1(4,:),1/2*order2(4,:),1/6*order3(4,:),1/24*order4(4,:)];

    if startopposite
        dec_k_opposite=[KSTARopposite,order1_opposite(1,:),1/2*order2_opposite(1,:),1/6*order3_opposite(1,:),1/24*order4_opposite(1,:)];
        dec_stock_opposite=[stockSTARopposite,order1_opposite(2,:),1/2*order2_opposite(2,:),1/6*order3_opposite(2,:),1/24*order4_opposite(2,:)];
        dec_l_opposite=[LSTARopposite,order1_opposite(3,:),1/2*order2_opposite(3,:),1/6*order3_opposite(3,:),1/24*order4_opposite(3,:)];
        dec_c_opposite=[CSTARopposite,order1_opposite(4,:),1/2*order2_opposite(4,:),1/6*order3_opposite(4,:),1/24*order4_opposite(4,:)];
    end
    
end

rng(13466910+randomseq,"twister");
rho_Z = normrnd(0,sigma_Z,[1 T]);
if ~(string(shock) == "none"||string(shock) == "historical"|| string(shock) == "historical_endogenous_P")
    rho_Z(1:5)=shock;
end

rng(123142+randomseq,"twister");
rho_phi = normrnd(0,sigma_phi,[1 T]);

rng(123140+randomseq,"twister");
rho_P = normrnd(0,sigma_P,[1 T]);

if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
    T = 37;
    realTFPshocks = readmatrix(tfpfile);
    rho_Z = zeros([1 T]) + realTFPshocks(realTFPshocks(:,1)>1980.5&realTFPshocks(:,1)<2017.5,2)';
    if string(shock) == "historical"
        realProfits = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/ProfitShare.csv");
        True_P_Path = (1./(1-realProfits(realProfits(:,1)>1980.5&realProfits(:,1)<2017.5,2)))';
        historical_Z_path = zeros([1 T]) + ZSTAR;
        rho_P = NaN([1 T]);
        for i=1:T
            historical_Z_path(i)=historical_Z_path(max(i-1,1))^LAMBDAZ*exp(rho_Z(i));
            if i == 1
                if startopposite
                    historical_P_implied_by_Z = P_func(G,PSTARopposite,LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
                else
                    historical_P_implied_by_Z = P_func(G,PSTAR,LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
                end
            else
                historical_P_implied_by_Z = P_func(G,True_P_Path(max(i-1,1)),LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
            end
            rho_P(i) = log(True_P_Path(i)/historical_P_implied_by_Z);
        end
    end
    if string(shock) == "historical_endogenous_P"
        rho_P = zeros([1 T]);
    end
end


shock_character_vector = char(shock);

if string(shock_character_vector(1:3)) == "irf"
    rho_Z = zeros([1 T]);
    rho_P = zeros([1 T]);
    if string(shock_character_vector(5)) == "Z"
        rho_Z(T-(hardcode_irf_T-2)) = sigma_Z;
    end
    if string(shock_character_vector(5)) == "z"
        rho_Z(T-(hardcode_irf_T-2)) = sigma_Z;
        rho_Z(T-(hardcode_irf_T-3)) = log(1/exp(sigma_Z)^LAMBDAZ);
    end
    if string(shock_character_vector(5)) == "P"
        rho_P(T-(hardcode_irf_T-2)) = sigma_P;
    end
    if string(shock_character_vector(5)) == "K"
        k1_mult = 0.5;
    end
end
        
k_sim = zeros([1 T]) + KSTAR*k0_mult;
Z_sim = zeros([1 T]) + ZSTAR;
phi_sim = zeros([1 T]) + PHISTAR;
P_sim = zeros([1 T]) + PSTAR;
if startopposite
    P_sim = zeros([1 T]) + PSTARopposite;
end

state_vars = repmat([k_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),Z_sim(1),P_sim(1),phi_sim(1)],T,1);
c_sim = NaN([1 T]);
l_sim = NaN([1 T]);
stock_sim = NaN([1 T]);

for i=1:T
    Z_sim(i)=Z_sim(max(i-1,1))^LAMBDAZ*exp(rho_Z(i));
    phi_sim(i) = phi_sim(max(i-1,1))^LAMBDAphi * exp(rho_phi(i));
%     in_opposite_regime = mod(floor((i-1)/regime_change_frequency),2) == 1;

    P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
    
    state_vars(i,:) = [k_sim(i),Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),P_sim(i),phi_sim(i)];
        
    l_sim(i)=decision_func_to_use(dec_l_hat,state_vars(i,:),ssvals,sigma_Z,sigma_P);
    
    c_sim(i)=decision_func_to_use(dec_c,[state_vars(i,:),l_sim(i)],[ssvals,LSTAR],sigma_Z,sigma_P);
    stock_sim(i)=decision_func_to_use(dec_stock,[state_vars(i,:),l_sim(i)],[ssvals,LSTAR],sigma_Z,sigma_P);
    if i < T
        k_sim(i+1) = decision_func_to_use(dec_k,[state_vars(i,:),l_sim(i)],[ssvals,LSTAR],sigma_Z,sigma_P);
    end
    if i+1 == 5 && (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
        k_sim(i+1) = KSTAR*k0_mult;
    end
    
end

w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA,phihat_sim);
r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA,phihat_sim);
y_sim = y_func(k_sim,l_sim,Z_sim,ALFA,phi_sim);
g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
profits_sim = ((P_sim-1)./P_sim).*y_sim;

rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',stock_sim',y_sim',phi_sim',phihat_sim',P_sim',Z_sim'];

if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
    rgwkcl_mat = rgwkcl_mat(5:37,:);
end


if string(shock_character_vector(1:3)) == "irf"
    rgwkcl_mat = rgwkcl_mat((T-(hardcode_irf_T-1)):T,:);
end

writematrix(rgwkcl_mat,"C:/Users/Nathan/Downloads/PerturbationMethods/rgwkcl_mat.xlsx",'Sheet',1) %,'Range','B4')