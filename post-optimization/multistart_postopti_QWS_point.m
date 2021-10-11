%script to calculate the post-optimized minima for QWS-symmetry over the
% dfined modulation index range for a spefific pulse number

% clear all;
% clc;

np = 3; % number of switching angles within one QWS
m = 1.24;% modulation index
the=pi;

R = 0.14;
W = 2*100*pi;
Ld = 0.0091;
Lq = 0.0229;
Psi_f = 1.0696;

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
%% loop over modulation index

if multirun == 1
    
    % Define modulation index array for post-optimization
    m_array = [0.6:0.01:1.2];
    column=38;
    
    % Dividing the defined array in half starting from the center
    m_array1=[0.9:0.01:1.2];
    m_array2=[0.9:-0.01:0.6];
    
    % Define array to save post-optimized reults
    x_over_m = zeros(np, length(m_array));
    x0_over_m = zeros(np, length(m_array));
    WTHD_over_m = zeros(np, length(m_array));
    
    % Define theta_dq resolution
    th_a = [(135/180)*pi:(3/180)*pi:pi];      % 2nd half
    th_b = [(135/180)*pi:-(3/180)*pi:0.5*pi]; % 1st half
    
    theta = [0.5*pi:(1/60)*pi:pi];            % Complete theta_dq array
    
    x_al =0.09;  % allowed change in switching angles between neighbouring points
    
    % define the array for saving the post-optimized switching angles
    alpha_1_over_m_theta = zeros(length(theta),length(m_array));
    alpha_2_over_m_theta = zeros(length(theta),length(m_array));
    alpha_3_over_m_theta = zeros(length(theta),length(m_array));
    
    % Provide the data from the local minima list saved from the
    % discontinious set of data. provide the local minima for the point modulation index '0.9' and theta_dq 135 degree or '0.75*pi'. 
    
    
   % The starting point is provided to the centre of the surface (16,31)
   % for the chosen area. 
    alpha_1_over_m_theta(16,31) = 0.4934;    
    alpha_2_over_m_theta(16,31) = 0.6419;
    alpha_3_over_m_theta(16,31) = 0.8930;
    
    alpha_1a_theta=alpha_1_over_m_theta(16,31);
    alpha_2a_theta=alpha_2_over_m_theta(16,31);
    alpha_3a_theta=alpha_3_over_m_theta(16,31);
    
    alpha_1b_theta=alpha_1_over_m_theta(16,31);
    alpha_2b_theta=alpha_2_over_m_theta(16,31);
    alpha_3b_theta=alpha_3_over_m_theta(16,31);
    
    count = 0;
    
    
    %  Calculates all the post optimized switching angles for constant modulation index with varying theta_dq. first 2 loops calculate the entire switching anlges in the first line.
    % Loop 1:
    
    for th = 2 : length(th_a)
        
        lb_array_alpha1 = [alpha_1a_theta-x_al];
        lb_array_alpha2 = [alpha_2a_theta-x_al];
        lb_array_alpha3 = [alpha_3a_theta-x_al];
        
        ub_array_alpha1 = [alpha_1a_theta+x_al];
        ub_array_alpha2 = [alpha_2a_theta+x_al];
        ub_array_alpha3 = [alpha_3a_theta+x_al];
        
        lb = [max(lb_array_alpha1);...
            max(lb_array_alpha2);...
            max(lb_array_alpha3)];
        
        ub = [min(ub_array_alpha1);...
            min(ub_array_alpha2);...
            min(ub_array_alpha3)];
        
        lb(lb<0) =0;
        ub(ub>(pi/2))=pi/2;
        
        x0=(lb+ub)/2;
        
        % constant value of modulation index=0.9 for the given example
        [problem] = createProblem(0.9,th_a(th),np,x0,delta_alpha_min,lb,ub);
        
        ms = MultiStart;
        ms.UseParallel = true;%use multiple processors
        ms.MaxTime = 10800; %maximum time to calculate: 3 h
        %calculate optimum switching angles
        [x,f,flag,outpt,allmins] = run(ms,problem,1000);%use 1000 different start points
%         %result data
%         local_min{th} = allmins;% save all local minima
%         x_over_m(:,th) = x;
%         x0_over_m(:,th) = x0;
%         WTHD_over_m(:,th) = f;
        alpha_1_over_m_theta(15+th,31) = x(1,1);
        alpha_2_over_m_theta(15+th,31) = x(2,1);
        alpha_3_over_m_theta(15+th,31) = x(3,1);
        
        alpha_1a_theta=x(1,1);
        alpha_2a_theta=x(2,1);
        alpha_3a_theta=x(3,1);
        
        clear('x','x0','f');
        
    end
    
    % Loop 2:
    
    for th = 2 : length(th_b)
        
        lb_array_alpha1 = [alpha_1b_theta-x_al];
        lb_array_alpha2 = [alpha_2b_theta-x_al];
        lb_array_alpha3 = [alpha_3b_theta-x_al];
        
        ub_array_alpha1 = [alpha_1b_theta+x_al];
        ub_array_alpha2 = [alpha_2b_theta+x_al];
        ub_array_alpha3 = [alpha_3b_theta+x_al];
        
        lb = [max(lb_array_alpha1);...
            max(lb_array_alpha2);...
            max(lb_array_alpha3)];
        
        ub = [min(ub_array_alpha1);...
            min(ub_array_alpha2);...
            min(ub_array_alpha3)];
        
        lb(lb<0) =0;
        ub(ub>(pi/2))=pi/2;
        
        x0=(lb+ub)/2;
        
        % constant value of modulation index=0.9 for the given example
        [problem] = createProblem(0.9,th_b(th),np,x0,delta_alpha_min,lb,ub);
        
        ms = MultiStart;
        ms.UseParallel = true;%use multiple processors
        ms.MaxTime = 10800; %maximum time to calculate: 3 h
        %calculate optimum switching angles
        [x,f,flag,outpt,allmins] = run(ms,problem,1000);%use 1000 different start points
%         %result data
%         local_min{th} = allmins;% save all local minima
%         x_over_m(:,th) = x;
%         x0_over_m(:,th) = x0;
%         WTHD_over_m(:,th) = f;
        alpha_1_over_m_theta(17-th,31) = x(1,1);
        alpha_2_over_m_theta(17-th,31) = x(2,1);
        alpha_3_over_m_theta(17-th,31) = x(3,1);
        
        alpha_1b_theta=x(1,1);
        alpha_2b_theta=x(2,1);
        alpha_3b_theta=x(3,1);
        
        clear('x','x0','f');
        
    end
    
    alpha_1a_over_m_theta=zeros(length(theta),length(m_array1));
    alpha_2a_over_m_theta=zeros(length(theta),length(m_array1));
    alpha_3a_over_m_theta=zeros(length(theta),length(m_array1));
    
    alpha_1b_over_m_theta=zeros(length(theta),length(m_array1));
    alpha_2b_over_m_theta=zeros(length(theta),length(m_array1));
    alpha_3b_over_m_theta=zeros(length(theta),length(m_array1));
    
    
    alpha_1a_over_m_theta(:,1) = alpha_1_over_m_theta(:,31);
    alpha_2a_over_m_theta(:,1) = alpha_2_over_m_theta(:,31);
    alpha_3a_over_m_theta(:,1) = alpha_3_over_m_theta(:,31);
    
    alpha_1b_over_m_theta(:,1) = alpha_1_over_m_theta(:,31);
    alpha_2b_over_m_theta(:,1) = alpha_2_over_m_theta(:,31);
    alpha_3b_over_m_theta(:,1) = alpha_3_over_m_theta(:,31);
    
    % The next loop calculates all the points that are in the right half of the
    % pláne, for the given example all the points right to '0.9' modulation index value..
    for count = 2:length(m_array1)
        for th=1:length(theta)
            
            
            disp(count);
            disp(th);
            
            
            % dependency on 2 neighbouring points
            if th==1
                
                lb_array_alpha1 = [alpha_1a_over_m_theta(th,count-1)-x_al;...
                    alpha_1a_over_m_theta(th+1,count-1)-x_al];
                
                
                lb_array_alpha2 = [alpha_2a_over_m_theta(th,count-1)-x_al;...
                    alpha_2a_over_m_theta(th+1,count-1)-x_al];
                
                
                lb_array_alpha3 = [alpha_3a_over_m_theta(th,count-1)-x_al;...
                    alpha_3a_over_m_theta(th+1,count-1)-x_al];
                
                ub_array_alpha1 = [alpha_1a_over_m_theta(th,count-1)+x_al;...
                    alpha_1a_over_m_theta(th+1,count-1)+x_al];
                
                
                ub_array_alpha2 = [alpha_2a_over_m_theta(th,count-1)+x_al;...
                    alpha_2a_over_m_theta(th+1,count-1)+x_al];
                
                
                ub_array_alpha3 = [alpha_3a_over_m_theta(th,count-1)+x_al;...
                    alpha_3a_over_m_theta(th+1,count-1)+x_al];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
                
                % dependency on 3 neighbouring points
            elseif th== length(theta)
                
                lb_array_alpha1 = [alpha_1a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_1a_over_m_theta(th,count-1)-x_al;...
                    alpha_1a_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha2 = [alpha_2a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_2a_over_m_theta(th,count-1)-x_al;...
                    alpha_2a_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha3 = [alpha_3a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_3a_over_m_theta(th,count-1)-x_al;...
                    alpha_3a_over_m_theta(th-1,count)-x_al;];
                
                ub_array_alpha1 = [alpha_1a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_1a_over_m_theta(th,count-1)+x_al;...
                    alpha_1a_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha2 = [alpha_2a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_2a_over_m_theta(th,count-1)+x_al;...
                    alpha_2a_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha3 = [alpha_3a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_3a_over_m_theta(th,count-1)+x_al;...
                    alpha_3a_over_m_theta(th-1,count)+x_al;];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
                
                
            else    % dependency on 4 neighbouring points
                lb_array_alpha1 = [alpha_1a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_1a_over_m_theta(th,count-1)-x_al;...
                    alpha_1a_over_m_theta(th+1,count-1)-x_al;...
                    alpha_1a_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha2 = [alpha_2a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_2a_over_m_theta(th,count-1)-x_al;...
                    alpha_2a_over_m_theta(th+1,count-1)-x_al;...
                    alpha_2a_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha3 = [alpha_3a_over_m_theta(th-1,count-1)-x_al;...
                    alpha_3a_over_m_theta(th,count-1)-x_al;...
                    alpha_3a_over_m_theta(th+1,count-1)-x_al;...
                    alpha_3a_over_m_theta(th-1,count)-x_al;];
                
                ub_array_alpha1 = [alpha_1a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_1a_over_m_theta(th,count-1)+x_al;...
                    alpha_1a_over_m_theta(th+1,count-1)+x_al;...
                    alpha_1a_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha2 = [alpha_2a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_2a_over_m_theta(th,count-1)+x_al;...
                    alpha_2a_over_m_theta(th+1,count-1)+x_al;...
                    alpha_2a_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha3 = [alpha_3a_over_m_theta(th-1,count-1)+x_al;...
                    alpha_3a_over_m_theta(th,count-1)+x_al;...
                    alpha_3a_over_m_theta(th+1,count-1)+x_al;...
                    alpha_3a_over_m_theta(th-1,count)+x_al;];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
            end
            
            x0=(lb+ub)/2;
            
            [problem] = createProblem(m_array1(count),theta(th),np,x0,delta_alpha_min,lb,ub);
            
            ms = MultiStart;
            ms.UseParallel = true;%use multiple processors
            ms.MaxTime = 10800; %maximum time to calculate: 3 h
            %calculate optimum switching angles
            [x,f,flag,outpt,allmins] = run(ms,problem,1000);%use 1000 different start points
            %result data
            local_min{count} = allmins;% save all local minima
            x_over_m(:,count) = x;
            x0_over_m(:,count) = x0;
            WTHD_over_m(:,count) = f;
            alpha_1_over_m_theta(th,30+count) = x(1,1);
            alpha_2_over_m_theta(th,30+count) = x(2,1);
            alpha_3_over_m_theta(th,30+count) = x(3,1);
            
            
            alpha_1a_over_m_theta(th,count) = x(1,1);
            alpha_2a_over_m_theta(th,count) = x(2,1);
            alpha_3a_over_m_theta(th,count) = x(3,1);
            
            clear('x','x0','f');
        end
        
    end
    
    % The next loop calculates all the points that are in the left half of the
    % pláne, for the given example all the points left to '0.9' modulation index value..
    for count = 2:length(m_array2)
        for th=1:length(theta)
            
            disp(count);
            disp(th);
            if th==1     % dependency on 2 neighbouring points
                
                lb_array_alpha1 = [alpha_1b_over_m_theta(th,count-1)-x_al;...
                    alpha_1b_over_m_theta(th+1,count-1)-x_al];
                
                
                lb_array_alpha2 = [alpha_2b_over_m_theta(th,count-1)-x_al;...
                    alpha_2b_over_m_theta(th+1,count-1)-x_al];
                
                
                lb_array_alpha3 = [alpha_3b_over_m_theta(th,count-1)-x_al;...
                    alpha_3b_over_m_theta(th+1,count-1)-x_al];
                
                ub_array_alpha1 = [alpha_1b_over_m_theta(th,count-1)+x_al;...
                    alpha_1b_over_m_theta(th+1,count-1)+x_al];
                
                
                ub_array_alpha2 = [alpha_2b_over_m_theta(th,count-1)+x_al;...
                    alpha_2b_over_m_theta(th+1,count-1)+x_al];
                
                
                ub_array_alpha3 = [alpha_3b_over_m_theta(th,count-1)+x_al;...
                    alpha_3b_over_m_theta(th+1,count-1)+x_al];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
                
            elseif th== length(theta)       % dependency on 3 neighbouring points
                
                lb_array_alpha1 = [alpha_1b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_1b_over_m_theta(th,count-1)-x_al;...
                    alpha_1b_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha2 = [alpha_2b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_2b_over_m_theta(th,count-1)-x_al;...
                    alpha_2b_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha3 = [alpha_3b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_3b_over_m_theta(th,count-1)-x_al;...
                    alpha_3b_over_m_theta(th-1,count)-x_al;];
                
                ub_array_alpha1 = [alpha_1b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_1b_over_m_theta(th,count-1)+x_al;...
                    alpha_1b_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha2 = [alpha_2b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_2b_over_m_theta(th,count-1)+x_al;...
                    alpha_2b_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha3 = [alpha_3b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_3b_over_m_theta(th,count-1)+x_al;...
                    alpha_3b_over_m_theta(th-1,count)+x_al;];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
                
                
            else        % dependency on 4 neighbouring points
                lb_array_alpha1 = [alpha_1b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_1b_over_m_theta(th,count-1)-x_al;...
                    alpha_1b_over_m_theta(th+1,count-1)-x_al;...
                    alpha_1b_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha2 = [alpha_2b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_2b_over_m_theta(th,count-1)-x_al;...
                    alpha_2b_over_m_theta(th+1,count-1)-x_al;...
                    alpha_2b_over_m_theta(th-1,count)-x_al;];
                
                lb_array_alpha3 = [alpha_3b_over_m_theta(th-1,count-1)-x_al;...
                    alpha_3b_over_m_theta(th,count-1)-x_al;...
                    alpha_3b_over_m_theta(th+1,count-1)-x_al;...
                    alpha_3b_over_m_theta(th-1,count)-x_al;];
                
                ub_array_alpha1 = [alpha_1b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_1b_over_m_theta(th,count-1)+x_al;...
                    alpha_1b_over_m_theta(th+1,count-1)+x_al;...
                    alpha_1b_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha2 = [alpha_2b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_2b_over_m_theta(th,count-1)+x_al;...
                    alpha_2b_over_m_theta(th+1,count-1)+x_al;...
                    alpha_2b_over_m_theta(th-1,count)+x_al;];
                
                ub_array_alpha3 = [alpha_3b_over_m_theta(th-1,count-1)+x_al;...
                    alpha_3b_over_m_theta(th,count-1)+x_al;...
                    alpha_3b_over_m_theta(th+1,count-1)+x_al;...
                    alpha_3b_over_m_theta(th-1,count)+x_al;];
                
                lb = [max(lb_array_alpha1);...
                    max(lb_array_alpha2);...
                    max(lb_array_alpha3)];
                
                ub = [min(ub_array_alpha1);...
                    min(ub_array_alpha2);...
                    min(ub_array_alpha3)];
                
                lb(lb<0) =0;
                ub(ub>(pi/2))=pi/2;
            end
            
            x0=(lb+ub)/2;
            
            [problem] = createProblem(m_array2(count),theta(th),np,x0,delta_alpha_min,lb,ub);
            
            ms = MultiStart;
            ms.UseParallel = true;%use multiple processors
            ms.MaxTime = 10800; %maximum time to calculate: 3 h
            %calculate optimum switching angles
            [x,f,flag,outpt,allmins] = run(ms,problem,1000);%use 1000 different start points
            %result data
            local_min{count} = allmins;% save all local minima
            x_over_m(:,count) = x;
            x0_over_m(:,count) = x0;
            WTHD_over_m(:,count) = f;
            alpha_1_over_m_theta(th,count) = x(1,1);
            %         alpha_2_over_m_theta(th,count) = x(2,1);
            %         alpha_3_over_m_theta(th,count) = x(3,1);
            %         %alpha_4_over_m_theta(th,count) = x(4,1);
            %         clear('x','x0','f');
            %
            alpha_1_over_m_theta(th,32-count) = x(1,1);
            alpha_2_over_m_theta(th,32-count) = x(2,1);
            alpha_3_over_m_theta(th,32-count) = x(3,1);
            
            
            alpha_1b_over_m_theta(th,count) = x(1,1);
            alpha_2b_over_m_theta(th,count) = x(2,1);
            alpha_3b_over_m_theta(th,count) = x(3,1);
            clear('x','x0','f');
        end
        
    end
    
    
elseif multirun == 0
    x0 = init_angles_OPP_quarter_IPMSM_QWS_theta_m(m_array,theta, np, R, Ld, Lq, W,Psi_f);
    [problem,x0] = createProblem(m, np,x0,delta_alpha_min);
    
    ms = MultiStart;
    ms.UseParallel = true;%use more multiple processors
    ms.MaxTime = 10800; %maximum time to calculate: 3 h
    %calculate optimum switching angles
    [x,f,flag,outpt,allmins] = run(ms,problem,100);%use 50 different start points
    
end




%% objective function to minimize
function [problem,x0] = createProblem(m,the,np,x0,delta_alpha_min,lb,ub)



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

opts = optimoptions(@fmincon,'Algorithm','sqp');

problem = createOptimProblem('fmincon',...
    'x0',x0,'objective',@minfunction, ...
    'Aineq',A,'bineq',b,'lb' , lb, 'ub',ub,'nonlcon',nonlcon,'options', opts);


    function [WTHD] = minfunction(x)
        
        R = 0.14;
        W = 2*100*pi;
        Ld = 0.0091;
        Lq = 0.0229;
        Psi_f = 1.0696;
        
        u_k = 0;
        
        u_1 = 0;
        
        for j = 1: np
            
            u_1 = u_1 + (-1)^(j+1)*cos(x(j));
        end
        
        u_1 = 4/(pi)*1500*(u_1);
        
        I_q1 = ( R*( u_1*sin(the) - (W*Psi_f)) - (W*Ld*u_1*cos(the)))/((R^2)  +  ((W^2)*(Ld*Lq)));
        
        I_d1 = ( (u_1 *cos(the)) + (W*Lq*I_q1) ) / R;
        
        % Fundamenatal current
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
        
        d = 0;
        B1=0;
        B2=0;
        
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
                U_sk = U_sk + (1500*(4/(pi))*(1/k_p(kk))*(((-1)^(j+1))*cos(x(j)*k_p(kk))));
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
        
        for kk =1:length(k_n) % for 6n-1 harmonic
            U_sk = 0; B1 = 0; B2 = 0; Id = 0; Iq = 0;
            
            for j = 1: np
                U_sk = U_sk + (1500*(4/pi)*(1/k_n(kk))*(((-1)^(j+1))*cos(x(j)*k_n(kk))));
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
        
        d= sqrt(sum(Acc_Ih));   % Distortion current
        
        WTHD = (d);
        
    end

%% function to set the nonlinear constraint of the modulation index
    function [ c, ceq ] = nonlconstraints(x)
        
        %sum over np pulses
        u_1 = 0;
        for j = 1: np
            
            u_1 = u_1 + ((-1)^(1+j))*cos(x(j));
            
        end
        
        u_1 = (4/(pi))*(u_1);
        %equality constraint c(x) = 0
        ceq = u_1- m;
        
        %inequality constraint c(x) <= 0
        c = [];
        
    end

end