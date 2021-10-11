%script to calculate the global and local minima for QWS-symmetry over the
%whole modulation index range for a spefific pulsenumber

% clear all;
% clc;

M = 8; % total number of switching events within one HWS (twice of QWS), eg: 2, 4, 6, 8 etc
m = 1.24;% modulation index

%Machine data
R = 0.14;
W = 2*100*pi;
Ld = 0.0091;
Lq = 0.0229;
Psi_f = 1.0696;
np=M;
%minimum allowed angle difference

delta_t_sw_min = 10e-6;% minimum time for one switch on switch off process
n_max = 6000; % maximum speed
p = 2;
%worst case estimation for maximum difference angle
delta_alpha_min_deg = ceil(n_max/60*p*2*pi*delta_t_sw_min*2*pi/360);
delta_alpha_min = delta_alpha_min_deg*2*pi/360;
%% calculate optimum switching angles[rad]
%multirun over all modulation indexes --< yes: 1, no: 0
multirun = 1;
 
%%
% % loop over modulation index

if multirun == 1
    
    %Define modulation index and theta_dq array resolution
    m_array= [0.01:0.01:1.27];
    theta_D=[90 : 1 : 180];
    theta_R=[(90/180)*pi : (pi/2) : pi ];
    x_over_m = zeros(np, length(m_array));
    x0_over_m = zeros(np, length(m_array));
    WTHD_over_m=zeros(np, length(m_array));
    %local_min = zeros(1, length(m_array));
    
    count = 0;
    th=0;
    
    for th = 1:length(theta_R)
        for count = 1:length(m_array)
            
            disp(count);
            
            
            %find initial value
            if count == 1
                
                if np == 1
                    
                    x0= rand *pi;
                    
                else
                    
                  %  x0 = init_angles_OPP_half_IPMSM(m_array(count),theta_R(th) , np, R, Ld,Lq,W,Psi_f);
                  % Can be extended for other values of switching
                  % events,eg: 4, 10, 12 etc
                  
                  if np==6
                  x0(1)=x_over_m_QWS_1(1,1);     %data from QWS for 3 switching events per quarter and extended for HWS. Similarly can be done for other switching as well
                   x0(2)=x_over_m_QWS_2(1,1);
                   x0(3)=x_over_m_QWS_3(1,1);
                   x0(4)=pi-x_over_m_QWS_3(1,1);
                   x0(5)=pi-x_over_m_QWS_2(1,1);
                   x0(6)=pi-x_over_m_QWS_1(1,1);
                  elseif np==8
                   x0(1)=x_over_m_QWS_1(1,1);    
                   x0(2)=x_over_m_QWS_2(1,1);
                   x0(3)=x_over_m_QWS_3(1,1);
                   x0(4)=x_over_m_QWS_4(1,1);
                   x0(5)=pi-x_over_m_QWS_4(1,1);
                   x0(6)=pi-x_over_m_QWS_3(1,1);
                   x0(7)=pi-x_over_m_QWS_2(1,1);
                   x0(8)=pi-x_over_m_QWS_1(1,1);
                  end
                  end
            else
                
                x0 =  x_over_m(:,count-1); % take previous optimum angle as initial value for next optimum calculation
                
            end
            
            [problem] = createProblem(m_array(count), theta_R(th), np,x0,delta_alpha_min);
            
            ms = MultiStart;
            ms.UseParallel = true;%use multiple processors
            ms.MaxTime = 10800; %maximum time to calculate: 3 h
            %calculate optimum switching angles
            [x,f,flag,outpt,allmins] = run(ms,problem,2000);%use 1000 different start points
            %result data
            local_min{count} = allmins;% save all local minima
            x_over_m(:,count) = x;
            x0_over_m(:,count) = x0;
            WTHD_over_m(:,count) = f;
            clear('x','x0','f');
            
            
        end
        
        save(['results_OPP_HW_np' num2str(np) '_theta_' num2str(theta_D(th)) '_localmin_min_diff_angle'],'m_array','x0_over_m','x_over_m','WTHD_over_m','local_min');
        
    end
    
    
elseif multirun == 0
    x0 = init_angles_OPP_quarter_IPMSM_halfwave_nina(m_array, np, R, Ld, Lq, W,Psi_f);
    [problem,x0] = createProblem(m, np,x0,delta_alpha_min);
    
    ms = MultiStart;
    ms.UseParallel = true;%use more multiple processors
    ms.MaxTime = 10800; %maximum time to calculate: 3 h
    %calculate optimum switching angles
    [x,f,flag,outpt,allmins] = run(ms,problem,100);%use 50 different start points
    
end

%% objective function to minimize
function [problem,x0] = createProblem(m,theta1,np,x0,delta_alpha_min)

nonlcon = @nonlconstraints;
%x(1)<x(2)<x(3)<x(np) in A*x<b -form
A = zeros(np-1, np);
b = -delta_alpha_min*ones(np-1,1);
for i= 1: (np -1)
    
    A(i,i) = 1;
    A(i,i+1) = -1;
    
end
Aeq = [];
beq = [];
%lower and upper bound for the calculated switching angles:0<x<pi
lb = zeros(1,np);
ub = pi * ones(1,np);


opts = optimoptions(@fmincon,'Algorithm','sqp');

problem = createOptimProblem('fmincon',...
    'x0',x0,'objective',@minfunction, ...
    'Aineq',A,'bineq',b,'lb' , lb, 'ub',ub,'nonlcon',nonlcon,'options', opts);


    function [Distortion] = minfunction(x)
        
        R = 0.14;
        W = 2*100*pi;
        Ld = 0.0091;
        Lq = 0.0229;
        Psi_f = 1.0696;
        
        u_k = 0;
        
        u_1 = 0; a_1=0; b_1=0;
        
        for j = 1: np
            
            a_1 = a_1 +((-1)^(j))*sin(x(j));
            
            b_1 = b_1 +((-1)^(1 + j))*cos(x(j));
        end
        
        angle_1=0;
        
        angle_1= atan2(b_1,a_1);
        
        the = (theta1)  ;
        
        u_1 = 2/(pi)*1500*(sqrt(a_1^2 + b_1^2));
        
        I_q1 = ( R*( u_1*sin(the) - (W*Psi_f)) - (W*Ld*u_1*cos(the)))/((R^2)  +  ((W^2)*(Ld*Lq)));
        
        I_d1 = ( (u_1 *cos(the)) + (W*Lq*I_q1) ) / R;
        
        %Fundamnetal current
        I1 = sqrt((I_d1^2) + (I_q1^2));
        
        
        %uneven and nontriplen harmonics up to k = 31
        k = [];
        k_n = [];
        k_p= [];
        
        for j = 1:5
            
            k_n(j)= 6*j-1;
            k_p(j) =6*j+1;
            
        end
        k = sort([k_n k_p]);
        
        Acc_Ud = [];
        Acc_Ud7 = [];
        Acc_Ud5 = [];
        Acc_Ih = [];
        Acc_Ih1 = [];
        Acc_Ih2 = [];
        Acc_init_Angle = [];
        Acc_Signal = [];
        Add_Signal = 0;
        
        Acc_7_cos7 = [];
        Acc_7_sin7 = [];
        Acc_5_cos7 = [];
        Acc_5_sin7 = [];
        
        
        for kk = 1:length(k_p) % for 6n+1 harmonic
            
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;a_k=0; b_k=0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j))*sin(x(j)*k_p(kk))));
                b_k = b_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j+1))*cos(x(j)*k_p(kk))));
                
            end
            
            angle_p=0;
            
            angle_p =  - atan2(b_k,a_k) + k_p(kk)*(angle_1 +  theta1);
            
            U_sk =(sqrt(a_k^2 +b_k^2));
            
            Ih1 = (1/2)*(   (U_sk)/(W*((6*kk)+1))  )*((1/Ld) -(1/Lq));
            Ih2 = (1/2)*(   (U_sk)/(W*((6*kk)+1))  )*((1/Ld) +(1/Lq));
            
            Acc_7_cos7 = [Acc_7_cos7 (Ih2*cos(angle_p))];
            Acc_7_sin7 = [Acc_7_sin7 -(Ih2*sin(angle_p))];
            Acc_5_cos7 = [Acc_5_cos7 (Ih1*cos(angle_p))];
            Acc_5_sin7 = [Acc_5_sin7 -(Ih1*sin(angle_p))];
            
            Acc_Ih2 = [Acc_Ih2 Ih1 Ih2];
            Acc_Ud7 = [Acc_Ud7 U_sk];                          
        end
        
        
        Acc_7_cos5 = [];
        Acc_7_sin5 = [];
        Acc_5_cos5 = [];
        Acc_5_sin5 = [];
        
        for kk =1:length(k_n) % for 6n-1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;a_k=0; b_k=0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j))*sin(x(j)*k_n(kk))));
                b_k = b_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(x(j)*k_n(kk))));
            end
            
            angle_n=0;
            
            angle_n =  - atan2(b_k,a_k) +  k_n(kk)*(angle_1 +  theta1 );
            
            U_sk =(sqrt(a_k^2 +b_k^2));
            
            Ih1 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) +(1/Lq));
            Ih2 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) -(1/Lq));
            
            Acc_7_cos5 = [Acc_7_cos5 (Ih2*cos(angle_n))];
            Acc_7_sin5 = [Acc_7_sin5 -(Ih2*sin(angle_n))];
            Acc_5_cos5 = [Acc_5_cos5 (Ih1*cos(angle_n))];
            Acc_5_sin5 = [Acc_5_sin5 -(Ih1*sin(angle_n))];
            
            Acc_Ud5 = [Acc_Ud5 U_sk];
            Acc_Ih1 = [Acc_Ih1 Ih1 Ih2];
          
        end
        
        %6n+1 current component
        Acc_7_cos7 = Acc_7_cos7 + Acc_7_cos5;
        Acc_7_sin7 = Acc_7_sin7 + Acc_7_sin5;
        
        %6n-1 current component
        Acc_5_cos5 = Acc_5_cos7 + Acc_5_cos5;
        Acc_5_sin5 = Acc_5_sin7 + Acc_5_sin5;
        
        Acc_7net = sqrt((Acc_7_cos7).^2+(Acc_7_sin7).^2);
        Acc_5net = sqrt((Acc_5_cos5).^2 + (Acc_5_sin5).^2);
        
        Acc_Ih = (Acc_7_cos7).^2+(Acc_7_sin7).^2 + (Acc_5_cos5).^2 + (Acc_5_sin5).^2;
        
        d= sqrt(sum(Acc_Ih)); %distortion current
        
        Distortion = (d);
        
    end

%% function to set the nonlinear constraint of the modulation index
    function [ c, ceq ] = nonlconstraints(x)
        
        %sum over np pulses
        u_1 = 0; a_1 = 0; b_1 =0;
        for j = 1: np
            
             a_1 = a_1 +((-1)^(j))*sin(x(j));
        
            b_1 = b_1 +((-1)^(1 + j))*cos(x(j));
        
        end
        
        u_1 =sqrt(a_1^2 + b_1^2);
        
        u_1 = (2/(pi))*(u_1);
        
        %equality constraint c(x) = 0
        ceq = u_1- m;
        
        %inequality constraint c(x) <= 0
        c = [];
        
    end

end