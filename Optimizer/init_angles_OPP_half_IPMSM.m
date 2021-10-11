%find good start angles x0 for the angles which satisfy the constraints
function [x0 ]= init_angles_OPP_half_IPMSM(m, theta1,np, R, Ld, Lq, W,Psi_f)


x0_array = cell(3000, 1);

WTHD = ones(3000,1)*10000;

tic
for i = 1:3000
    
    x0 = zeros(np,1);
    %create random values for angles between 0 and pi
    for j = 2:np
        
        x0(j)= rand *pi;
        
    end
    
    x0 = sort(x0);
    
    a_1 = 0;
    b_1 = 0;
    
    for j = 2:np
        a_1 = a_1 +((-1)^(j))*sin(x0(j));
        
        b_1 = b_1 +((-1)^(1 + j))*cos(x0(j));
    end
    
    angle_11 = atan2(a_1,-b_1);
    
    
    x0(1) = acos(((1/2)*(1 + a_1^2 + b_1^2 -  (((m*pi)/2)^2)))    /((a_1^2 +b_1^2 )^0.5)) - angle_11;
    
    
    if imag(x0(1))==0 && x0(1) > 0 && x0(1)<x0(2) && x0(1)< pi%check if calculated angle is valid
        
        a_1 = 0;
        b_1 = 0;
        
        for j = 1: np
            
            a_1 = a_1 +((2/(pi))*1500*((-1)^(j))*sin(x0(j)));
            
            b_1 = b_1 +((2/(pi))*1500*((-1)^(1 + j))*cos(x0(j)));
        end
        
        angle_1=0;
        
        angle_1= atan2(b_1,a_1);
        
        the = (theta1);
        
        u_1 = 2/(pi)*1500*(sqrt(a_1^2 + b_1^2));
        
        I_q1 = ( R*( u_1*sin(the) - (W*Psi_f)) - (W*Ld*u_1*cos(the)))/((R^2)  +  ((W^2)*(Ld*Lq)));
        I_d1 = ( (u_1 *cos(the)) + (W*Lq*I_q1) ) / R;
        
        %fundamental current
        I1 = sqrt((I_d1^2) + (I_q1^2));
        
        %uneven and nontriplen harmonics up to k = 49
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
        
        for kk = 1:length(k_p) % for 6n-1 harmonic
            
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0; a_k =0; b_k=0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j))*sin(x0(j)*k_p(kk))));
                b_k = b_k + ((3000/pi)*(1/k_p(kk))*(((-1)^(j+1))*cos(x0(j)*k_p(kk))));
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
        
        
        for kk =1:length(k_n) % for 6n+1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0; a_k=0; b_k=0;
            
            for j = 1: np
                a_k = a_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j))*sin(x0(j)*k_n(kk))));
                b_k = b_k + ((3000/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(x0(j)*k_n(kk))));
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
        
        Acc_7_cos7 = Acc_7_cos7 + Acc_7_cos5;
        Acc_7_sin7 = Acc_7_sin7 + Acc_7_sin5;
        
        Acc_5_cos5 = Acc_5_cos7 + Acc_5_cos5;
        Acc_5_sin5 = Acc_5_sin7 + Acc_5_sin5;
        
        Acc_7net = sqrt((Acc_7_cos7).^2+(Acc_7_sin7).^2);
        Acc_5net = sqrt((Acc_5_cos5).^2 + (Acc_5_sin5).^2);
        
        Acc_Ih = (Acc_7_cos7).^2+(Acc_7_sin7).^2 + (Acc_5_cos5).^2 + (Acc_5_sin5).^2;
        
        c= sqrt(sum(Acc_Ih));
        
        
        
        WTHD(i) = (c);
        
        
        x0_array{i}= x0;
        
    end
    
end

% find best set of angles by looking for the minimum WTHD
minidx = find (WTHD == min(WTHD));
x0 = x0_array{minidx};

loop_time =toc


end