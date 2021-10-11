
theta=[90: 1 :180];             % Define theta_dq array resolution for which the simulation has been carried out
m_array= [0.01:0.01:1.27];     % Define modulation index array resolution
M=4;                           % Define the number of switching events per quarter wave eg:3,4 etc

x_over_m_QWS_1=zeros;
x_over_m_QWS_2=zeros;
x_over_m_QWS_3=zeros;
x_over_m_QWS_4=zeros;
 
np=M;

for j =1:length(theta)
    
    % load(['results_NPC_OPP_type_A_QW_np3_' num2str(theta(j)) '_localmin_min_diff_angle1'])
    load(['results_OPP_type_A_QW_np' num2str(np) '_theta_' num2str(theta(j)) '_localmin_min_diff_angle'])
    
    for i=1:length(m_array)
        
        if np==3
            x_over_m_QWS_1(j,i)=x_over_m(1,i);
            x_over_m_QWS_2(j,i)=x_over_m(2,i);
            x_over_m_QWS_3(j,i)=x_over_m(3,i);
            
        elseif np==4
            x_over_m_QWS_1(j,i)=x_over_m(1,i);
            x_over_m_QWS_2(j,i)=x_over_m(2,i);
            x_over_m_QWS_3(j,i)=x_over_m(3,i);
            x_over_m_QWS_4(j,i)=x_over_m(4,i);
        end
        
    end
end

