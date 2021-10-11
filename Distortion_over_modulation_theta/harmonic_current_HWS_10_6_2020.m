%% Define modulation index and theta_dq
m_array1= [0.01:0.01:1.27];                         % define m_array
theta = [(90/180)*pi:(3/180)*pi:(180/180)*pi];      % define theta_dq

har_curr_HWS= zeros(length(theta),length(m_array1));% define total current array
har_d_HWS=zeros(length(theta),length(m_array1));    % define distortion current array
har_I_HWS=zeros(length(theta),length(m_array1));    % define fundamental current array

M=6;                                                % define total switching events per half wave

Acc_the = [];
Acc_u_123 = [];

% Define machine parameters
R = 0.14;
W = 2*(3000/30)*pi;                                 % Distortion current and total curent calculated for speed=3000rpm
Ld = 0.0091;
Lq = 0.0229;
Psi_f = 1.0696;
pa=[];
np=M;

for dis=1:length(theta)
    for i=1:length(m_array1)
        
        m_ref1=m_array1;
        theta_ref1=theta(dis);
        
        
        alpha(1,1) = x_over_m_HWS_1(dis,i);
        alpha(2,1) = x_over_m_HWS_2(dis,i);
        alpha(3,1) = x_over_m_HWS_3(dis,i);
        
        alpha(4,1) = x_over_m_HWS_4(dis,i);
        alpha(5,1) = x_over_m_HWS_5(dis,i);
        alpha(6,1) = x_over_m_HWS_6(dis,i);
        
        u_1=0;
        I_q1=0;
        I_d1=0;
        a_1=0;
        b_1=0;
        the=0;
        
        
        for j = 1: np
            
            a_1 = a_1 +((2/(pi))*((-1)^(j)))*sin(alpha(j,1));
            
            b_1 = b_1 +((2/(pi))*((-1)^(1 + j)))*cos(alpha(j,1));
            
        end
        
        angle_1=0;
        
        angle_1= atan2(b_1,a_1);
        
        the = (theta(dis))  ;
        
        u_1 = 1500*((a_1^2 +b_1^2)^0.5);
        
        Acc_ud(dis,i)=(u_1*cos(the));
        Acc_uq(dis,i)=(u_1*sin(the));
        
        Acc_the(dis,i)=angle_1;
        Acc_the_degree(dis,i)=angle_1*(180/pi);
        
        I_q1 = ( R*( u_1*sin(the) - (W*Psi_f)) - (W*Ld*u_1*cos(the)))/((R^2)  +  ((W^2)*(Ld*Lq)));
        
        I_d1 = ( (u_1 *cos(the)) + (W*Lq*I_q1) ) / R;
        
        %fundamental current
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
        ni=1;
        wei1=0;
        
        Acc_time_Signal = [];
        Add_time_Signal = [];
        
        
        Acc_Ud = [];
        Acc_Ud7_HWS = [];
        Acc_Ud5_HWS = [];
        
        Acc_Ih = [];
        Acc_Ih1 = [];
        Acc_Ih2 = [];
        
        Acc_7_cos7 = [];
        Acc_7_sin7 = [];
        Acc_5_cos7 = [];
        Acc_5_sin7 = [];
        
        
        for kk = 1:length(k_p) % for 6n+1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;angle1=0;
            a_k=0;
            b_k=0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j))*sin(alpha(j,1)*k_p(kk))));
                b_k = b_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j+1))*cos(alpha(j,1)*k_p(kk))));
            end
            
            
            angle_p=0;
            
            angle_p =  - atan2(b_k,a_k) + k_p(kk)*(angle_1 +  theta(dis));
            
            U_sk =(sqrt(a_k^2 +b_k^2));
            
            Ih1 = (1/2)*(   (U_sk)/(W*((6*kk)+1))  )*((1/Ld) -(1/Lq));
            Ih2 = (1/2)*(   (U_sk)/(W*((6*kk)+1))  )*((1/Ld) +(1/Lq));
            
            Acc_7_cos7 = [Acc_7_cos7 (Ih2*cos(angle_p))];
            Acc_7_sin7 = [Acc_7_sin7 -(Ih2*sin(angle_p))];
            Acc_5_cos7 = [Acc_5_cos7 (Ih1*cos(angle_p))];
            Acc_5_sin7 = [Acc_5_sin7 -(Ih1*sin(angle_p))];
            
            
            Acc_Ih2 = [Acc_Ih2 Ih1 Ih2];
            Acc_Ud7_HWS = [Acc_Ud7_HWS U_sk];
           
        end
        
        Acc_7_cos5 = [];
        Acc_7_sin5 = [];
        Acc_5_cos5 = [];
        Acc_5_sin5 = [];
        
        for kk = 1:length(k_n) % for 6n-1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;a_k=0;b_k=0; angle2 = 0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j))*sin(alpha(j,1)*k_n(kk))));
                b_k = b_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(alpha(j,1)*k_n(kk))));
            end
            
            angle_n=0;
            angle_n =  - atan2(b_k,a_k) +  k_n(kk)*(angle_1 +  theta(dis) );
            
            U_sk =(sqrt(a_k^2 +b_k^2));
            
            Ih1 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) +(1/Lq));
            Ih2 = (1/2)*(   (U_sk)/(W*((6*kk)-1))  )*((1/Ld) -(1/Lq));
            
            Acc_7_cos5 = [Acc_7_cos5 (Ih2*cos(angle_n))];
            Acc_7_sin5 = [Acc_7_sin5 -(Ih2*sin(angle_n))];
            Acc_5_cos5 = [Acc_5_cos5 (Ih1*cos(angle_n))];
            Acc_5_sin5 = [Acc_5_sin5 -(Ih1*sin(angle_n))];
            
            Acc_Ud5_HWS = [Acc_Ud5_HWS U_sk];
            Acc_Ih1 = [Acc_Ih1 Ih1 Ih2];
            
        end
        
        % 6n+1 current components
        Acc_7_cos7 = Acc_7_cos7 + Acc_7_cos5;
        Acc_7_sin7 = Acc_7_sin7 + Acc_7_sin5;
        
        % 6n-1 current components
        Acc_5_cos5 = Acc_5_cos7 + Acc_5_cos5;
        Acc_5_sin5 = Acc_5_sin7 + Acc_5_sin5;
        
        
        Acc_7net = sqrt((Acc_7_cos7).^2+(Acc_7_sin7).^2);
        Acc_5net = sqrt((Acc_5_cos5).^2 + (Acc_5_sin5).^2);
        
        Acc_Ih = (Acc_7_cos7).^2+(Acc_7_sin7).^2 + (Acc_5_cos5).^2 + (Acc_5_sin5).^2;
        
        c= sqrt(sum(Acc_Ih));               % Distortion current
        
        har_d_HWS(dis,i) = (c);             % Distortion current array over modulation index and theta_dq
        har_I_HWS(dis,i) = I1;              % Fundamental current array over modulation index and theta_dq
        har_curr_HWS(dis,i) = (100*c)/I1;   % Total current array over modulation index and theta_dq
        
    end
end