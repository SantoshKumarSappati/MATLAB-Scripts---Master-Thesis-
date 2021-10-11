%% Define modulation index and theta_dq
m_array= [0.01:0.01:1.27];                        % define m_array
theta = [(90/180)*pi:(1/180)*pi:(180/180)*pi];    % define theta_dq

har_curr= zeros(length(theta),length(m_array));   % define total current array
har_d=zeros(length(theta),length(m_array));       % define distortion current array
har_I=zeros(length(theta),length(m_array));       % define fundamental current array

M = 4  ;                                         % define total switching events per quarter wave
Acc_u_111=[];

%define machine parameters

R = 0.14;
W = 2*(3000/30)*pi;       % Distortion current and total curent calculated for speed=3000rpm
Ld = 0.0091;
Lq = 0.0229;
Psi_f = 1.0696;
pa=[];
np=M;

for dis=1:length(theta)
    
    for i=1:length(m_array)
        
        m_ref1=m_array;
        
        alpha(1,1) = x_over_m_QWS_1(dis,i);
        alpha(2,1) = x_over_m_QWS_2(dis,i);
        alpha(3,1) = x_over_m_QWS_3(dis,i);
        alpha(4,1) = x_over_m_QWS_4(dis,i);
        
        theta_ref1 = theta(dis);
        
        u_1=0;
        
        for j = 1: np
            
            u_1 = u_1 + (4/(pi))*(-1)^(j+1)*cos(alpha(j,1));
        end
        
        u_1 = 1500*(u_1);
        
        I_q1 = ( R*( u_1*sin(theta_ref1) - (W*Psi_f)) - (W*Ld*u_1*cos(theta_ref1)))/((R^2)  +  ((W^2)*(Ld*Lq)));
        
        I_d1 = ( (u_1 *cos(theta_ref1)) + (W*Lq*I_q1) ) / R;
        
        % Fundamnetal current
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
        
        c =0;
        d = 0;
        d1=0;
        d_har=[];
        B1=0;
        B2=0;
        
        Acc_Ud = [];
        Acc_Ud7 = [];
        Acc_Ud5 = [];
        Acc_Ih = [];
        Acc_Ih1 = [];
        Acc_Ih2 = [];
        
        
        Acc_7_cos7 = [];
        Acc_7_sin7 = [];
        Acc_5_cos7 = [];
        Acc_5_sin7 = [];
        
        
        for kk = 1:length(k_p) % for 6n+1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;
            
            for j = 1: np
                U_sk = U_sk + (1500*(4/(pi))*(1/k_p(kk))*(((-1)^(j+1))*cos(alpha(j,1)*k_p(kk))));
            end
            
            angle_p =(theta_ref1+0.5*pi)*k_p(kk);
            
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
        
        for kk = 1:length(k_n) % for 6n-1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;
            
            for j = 1: np
                U_sk = U_sk + (1500*(4/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(alpha(j,1)*k_n(kk))));
            end
            
            angle_n =(theta_ref1+0.5*pi)*k_n(kk);
            
            Ih1 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) +(1/Lq));
            Ih2 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) -(1/Lq));
            
            Acc_7_cos5 = [Acc_7_cos5 (Ih2*cos(angle_n))];
            Acc_7_sin5 = [Acc_7_sin5 -(Ih2*sin(angle_n))];
            Acc_5_cos5 = [Acc_5_cos5 (Ih1*cos(angle_n))];
            Acc_5_sin5 = [Acc_5_sin5 -(Ih1*sin(angle_n))];
            
            Acc_Ud5 = [Acc_Ud5 U_sk];
            Acc_Ih1 = [Acc_Ih1 Ih1 Ih2];
            
        end
        
        Acc_7_c7 = [];
        Acc_7_s7 = [];
        Acc_5_c5 = [];
        Acc_5_s5 = [];
        
        % 6n+1 current components
        Acc_7_c7 = Acc_7_cos7 + Acc_7_cos5;
        Acc_7_s7 = Acc_7_sin7 + Acc_7_sin5;
        
        net_7=  sqrt((Acc_7_c7).^2+(Acc_7_s7).^2);
        
        % 6n-1 current components
        Acc_5_c5 = Acc_5_cos7 + Acc_5_cos5;
        Acc_5_s5 = Acc_5_sin7 + Acc_5_sin5;
        
        net_5=  sqrt((Acc_5_c5).^2+(Acc_5_s5).^2);
        
        Acc_Ih = (Acc_7_c7).^2+(Acc_7_s7).^2 + (Acc_5_c5).^2 + (Acc_5_s5).^2;
        
        d= sqrt(sum(Acc_Ih)); % Distortion current
        
        har_d(dis,i) = (d);  % Distortion current array over modulation index and theta_dq
        
        har_I(dis,i) = I1;   % Fundamental current array over modulation index and theta_dq
        
        har_curr(dis,i) = (100*d)/I1;
        
    end
end
