%find good start angles x0 for the angles which satisfy the constraints
function [x0 ]= init_angles_OPP_quarter_IPMSM(m,the, np, R, Ld, Lq, W,Psi_f)


x0_array = cell(3000, 1);

WTHD = ones(3000,1)*10000;

tic
for i = 1:3000
    
    x0 = zeros(np,1);
    %create random values for angles between 0 and pi/2
    for j = 2:np
        
        x0(j)= rand *pi/2;
        
    end
    
    x0 = sort(x0);
    
    u_1 = 0;
    
    
    for j = 2:np
        u_1 = u_1 +(-1)^(j+1)*cos(x0(j));
    end
    
    %based on created angles, the first angle is calculated to fulfill the
    %constraint
    
    x0(1) = acos(((m*pi)/4) - u_1);
    
    
    if imag(x0(1))==0 && x0(1) > 0 && x0(1)<x0(2)%check if calculated angle is valid
        
        
        
        u_k = 0;
        
        u_1 = 0;
        
        
        for j = 1: np
            
            u_1 = u_1 + (-1)^(j+1)*cos(x0(j));
        end
        
        u_1 = 4/(pi)*1500*(u_1);
        
        
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
        
        
        for kk = 1:length(k_p) % for 6n-1 harmonic
            
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;
            
            for j = 1: np
                U_sk = U_sk + (1500*(4/(pi))*(1/k_p(kk))*(((-1)^(j+1))*cos(x0(j)*k_p(kk))));
            end
            
            angle_p =(the+0.5*pi)*k_p(kk);
            
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
        
        
        
        for kk =1:length(k_n) % for 6n+1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;
            
            for j = 1: np
                U_sk = U_sk + (1500*(4/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(x0(j)*k_n(kk))));
            end
            
            angle_n =(the+0.5*pi)*k_n(kk);
            
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
        
        Acc_7_c7 = Acc_7_cos7 + Acc_7_cos5;
        Acc_7_s7 = Acc_7_sin7 + Acc_7_sin5;
        
        Acc_5_c5 = Acc_5_cos7 + Acc_5_cos5;
        Acc_5_s5 = Acc_5_sin7 + Acc_5_sin5;
        
        Acc_Ih = (Acc_7_c7).^2+(Acc_7_s7).^2 + (Acc_5_c5).^2 + (Acc_5_s5).^2;
        
        d= sqrt(sum(Acc_Ih));
        
        WTHD(i) = (d);
        x0_array{i}= x0;
        
    end
    
end

% find best set of angles by looking for the minimum WTHD
minidx = find (WTHD == min(WTHD));
x0 = x0_array{minidx};

loop_time =toc;


end