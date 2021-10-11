
theta=[90: 3 :180];             % Define theta_dq array resolution for which the simulation has been carried out
m_array= [0.01:0.01:1.27];     % Define modulation index array resolution
M=6;                           % Define the number of switching events per quarter wave eg:3,4 etc

x_over_m_HWS_1=zeros;
x_over_m_HWS_2=zeros;
x_over_m_HWS_3=zeros;
x_over_m_HWS_4=zeros;
x_over_m_HWS_5=zeros; 
x_over_m_HWS_6=zeros;

np=M;

for j =1:length(theta)
    
    % load(['results_NPC_OPP_type_A_QW_np3_' num2str(theta(j)) '_localmin_min_diff_angle1'])
    load(['results_OPP_type_A_HW_np' num2str(np) '_theta_' num2str(theta(j)) '_localmin_min_diff_angle'])
    
    for i=1:length(m_array)
        
            x_over_m_HWS_1(j,i)=x_over_m(1,i);
            x_over_m_HWS_2(j,i)=x_over_m(2,i);
            x_over_m_HWS_3(j,i)=x_over_m(3,i);
            x_over_m_HWS_4(j,i)=x_over_m(4,i);
            x_over_m_HWS_5(j,i)=x_over_m(5,i);
            x_over_m_HWS_6(j,i)=x_over_m(6,i);
        
    end
end


%%note: done for 1 case where, 6 switching event per half wave and follows
%%half wave symmetry (lack of compuation time within thesis time frame)
