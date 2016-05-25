% This code is about Delta A/C flight condition (2) dynamics in longitudinal and lateral modes with all degrees of freedom of approximation 
% And use pole placement (Aker) and LQR methods to controll some states as follow in the code
% You must fill INPUT DATA section to get solution
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
%% INPUT DATA
        %% % Genral parameter
        S=576; % sing area (m^2)
        AR=7.75; % aspect ratio
        c=9.17; % chord (m)
        T_tot=730; % total related thurst (kN)
        CG=0.3*c; % C.G.
        I_X_P=25;
        I_Z_P=2.5;
        W=3e5; % weight (kg)
        g=9.81; % m/s^2
        %% % Inertias (kg m^2)
        I_xx=3.77e7;
        I_yy=4.31e7;
        I_zz=7.62e7;
        I_xz=3.35e6;
        %% Parameter
        H=6100; % height (m)
        M=0.6; % mach number
        u_0=190; % initial velocity (m/s)
        q_bar=11730; % dynamic prassure (Nm^2)
        alpha_0=2.2; % degree
        gama_0=0; % degree
        theta_0=alpha_0+gama_0; % degree
        %% longitudinal stability derivative
        X_u=-0.003;
        X_w=0.04;
        X_delta_e=0.26;
        X_delta_T=0.15e-4;
        Z_u=-0.08;
        Z_w=-0.618;
        Z_wdot=0; % assumed
        Z_q=0; % assumed
        Z_delta_e=-6.83;
        Z_delta_T=0.05e-5;
        M_u=3.28e-4;
        M_w=-0.007;
        M_wdot=-0.001;
        M_q=-0.77;
        M_delta_e=-1.25;
        M_delta_T=1.42e-5;
        %% lateral stability derivative
        Y_v=-0.11;
        Y_p=0;% assumed
        Y_r=0;%assumed
        Y_delta_a_star=-0.29e-4;
        Y_delta_r_star=0.0055;
        L_beta_dash=-1.33;
        L_p_dash=-1;
        L_r_dash=0.28;
        L_delta_a_dash=0.43;
        L_delta_r_dash=0.187;
        N_beta_dash=0.432;
        N_p_dash=-0.09;
        N_r_dash=-0.2;
        N_delta_a_dash=0.03;
        N_delta_r_dash=-0.52;
        %% Pole placement (Aker)
                %% Longitudinal Mode
                        %% State (u) with delta_e
                        requird_poles_u=[-1,-5,-2,-3];                        
        %% LQR paramater
                %% Longitudinal Mode
                        %% State (u) with delta_T
                        Q_u=diag([1000 0 0 0]);
                        R_u=1/1000^2;
                        %% State (theta) with delta_e
                        Q_theta=diag([0 0 0 100]);
                        R_theta=1/(20*pi/180);
                        %% State (w) with delta_e
                        Q_w=diag([1 100 0 0]);
                        R_w=1/(20*pi/180)^2;
                %% Latral Mode
                        %% State (psi) with delta_r
                        Q_psi=diag([0 0 0 0 1]);
                        R_psi=1/(20*pi/180);
        %% Plotting Control
        % Plot = 0 --> No Plotting data    Plot = 1 --> Plotting data
        Plot=1;  
        %% Save figures control
        % Save = 0 --> No Saving figures    Save = 1 --> Saving figures
        % The figures will saved in running folder directory
        Save=0;                         
%% -------------------------------------------------------------------------------------------%%
%% Longitudinal Mode
        %% exact linearized longitudinal dynamics
        % SS model
        % x = [u w q theta]'
        % u=[delta_e, delta_T]'
            A1_long=[1,                                     0,     0,      0;...
                             0,                  1/(1-Z_wdot),     0,      0;...
                             0,       M_wdot/(1-Z_wdot),    1,      0;...
                             0,                                      0,   0,      1];
            A2_long=[X_u,      X_w,                           0,        -g*cosd(theta_0);...
                             Z_u,      Z_w,               u_0+Z_q,        -g*sind(theta_0);...
                             M_u,    M_w,                      M_q,                                  0;...
                             0,              0,                            1,                                  0];
            B_long=[X_delta_e,     X_delta_T;...
                           Z_delta_e,     Z_delta_T;...
                           M_delta_e,    M_delta_T;...
                           0,                      0];
            A_long=A1_long*A2_long;
            B_long=A1_long*B_long;
            C_long=diag([1,1,1,1]);
            D_long=0;
        %% eigen vectors and eigen values
        [Eig_Vec_long,Eig_Val_long]=eig(A_long);
        %% eigen vectors norm
        for n=1:4
            for m=1:4
                Eig_Vec_long_Norm(n,m)=norm(Eig_Vec_long(n,m));
            end
        end
        %% TFs
        SS_long = ss(A_long,B_long,C_long,D_long);
        TF_long=tf(SS_long);
        [Sys_long_wn,Sys_long_zeta,Sys_long_poles]=damp(TF_long);
%% -------------------------------------------------------------------------------------------%%
%% linearized long period longitudinal dynamics
        %% SS x=[u theta],  u=[delta_e delta_T]
        A_long_period=[A_long(1,1),A_long(1,4);-A_long(2,1)/u_0,-A_long(2,4)/u_0];
        B_long_period=[B_long(1,:);-B_long(2,:)/u_0];
        C_long_period=eye(2);  % assume the output is al the states
        D_long_period=0;
        %% TF
        SS_long_period = ss(A_long_period,B_long_period,C_long_period,D_long_period);
        TF_long_period=tf(SS_long_period);
        [Sys_long_period_wn,Sys_long_period_zeta,Sys_long_period_poles]=damp(TF_long_period);
%% -------------------------------------------------------------------------------------------%%
%% linearized short period longitudinal dynamics
        %% SS x=[w q],  u=[delta_e delta_T]
        A_short_period=A_long(2:3,2:3);
        B_short_period=[B_long(2,:);B_long(3,:)];
        C_short_period=eye(2);  % assume the output is al the states
        D_short_period=0;
        %% TF
        SS_short_period = ss(A_short_period,B_short_period,C_short_period,D_short_period);
        TF_short_period =tf(SS_short_period);
        [Sys_short_period_wn,Sys_short_period_zeta,Sys_short_period_poles] =damp(TF_short_period );
%% -------------------------------------------------------------------------------------------%%
%% u/de  TF
        %% exact
        u_de_TF_exact=tf([TF_long.num{1}],[TF_long.den{1}]);
        %% long period
        u_de_TF_long_period=tf([TF_long_period.num{1}],[TF_long_period.den{1}]);
%% -------------------------------------------------------------------------------------------%%
%% w/de  TF
        %% exact
        w_de_TF_exact=tf([TF_long.num{2}],[TF_long.den{2}]);
        %% short period
        w_de_TF_short_period=tf([TF_short_period.num{1}],[TF_short_period.den{1}]);
%% -------------------------------------------------------------------------------------------%%
%% q/de  TF
        %% exact
        q_de_TF_exact=tf([TF_long.num{3}],[TF_long.den{3}]);
        %% short period
        q_de_TF_short_period=tf([TF_short_period.num{2}],[TF_short_period.den{2}]);
%% -------------------------------------------------------------------------------------------%%
%% Theta/de  TF
        %% exact
        theta_de_TF_exact=tf([TF_long.num{4}],[TF_long.den{4}]);
        %% long period
        theta_de_TF_long_period=tf([TF_long_period.num{2}],[TF_long_period.den{2}]);
%% -------------------------------------------------------------------------------------------%%
%% u/dT  TF
        %% exact
        u_dT_TF_exact=tf([TF_long.num{5}],[TF_long.den{5}]);
        %% long period
        u_dT_TF_long_period=tf([TF_long_period.num{3}],[TF_long_period.den{3}]);
%% -------------------------------------------------------------------------------------------%%
%% w/dT  TF
        %% exact
        w_dT_TF_exact=tf([TF_long.num{6}],[TF_long.den{6}]);
        %% short period
        w_dT_TF_short_period=tf([TF_short_period.num{3}],[TF_short_period.den{3}]);
%% -------------------------------------------------------------------------------------------%%
%% q/dT  TF
        %% exact
        q_dT_TF_exact=tf([TF_long.num{7}],[TF_long.den{7}]);
        %% short period
        q_dT_TF_short_period=tf([TF_short_period.num{4}],[TF_short_period.den{4}]);
%% -------------------------------------------------------------------------------------------%%
%% Theta/dT  TF
        %% exact
        theta_dT_TF_exact=tf([TF_long.num{8}],[TF_long.den{8}]);
        %% long period
        theta_dT_TF_long_period=tf([TF_long_period.num{4}],[TF_long_period.den{4}]);
%% -------------------------------------------------------------------------------------------%%
%% Lateral Mode        
        %% exact linearized lateral dynamics
        % SS model
        % x = [beta, p ,r ,phi, psi]
        % u = [delta_a, delta_r]
        A_latral = [Y_v,                       Y_p/u_0,                    -(1-Y_r/u_0),                   -g*cosd(theta_0)/u_0,               0;...
                          L_beta_dash,       L_p_dash,                         L_r_dash,                                                    0,               0;...
                          N_beta_dash,     N_p_dash,                        N_r_dash,                                                    0,               0;...
                          0,                                       1,                                     0,                                                    0,               0;...
                          0,                                       0,                                     1,                                                    0,               0];
        B_latral = [Y_delta_a_star,        Y_delta_r_star;...
                          L_delta_a_dash,      L_delta_r_dash;...
                          N_delta_a_dash,    N_delta_r_dash;...
                          0,                                                    0;...
                          0,                                                    0];
        C_latral = eye(5,5); %Outputs are All States
        D_latral = 0;
        %% eigen vectors and eigen values
        [Eig_Vec_latral,Eig_Val_latral]=eig(A_latral);
        %% eigen vectors norm
        for n=1:4
            for m=1:4
                Eig_Vec_latral_Norm(n,m)=norm(Eig_Vec_latral(n,m));
            end
        end
        %% TFs
        SS_latral = ss(A_latral,B_latral,C_latral,D_latral);
        TF_latral=tf(SS_latral);
        [Sys_latral_wn,Sys_latral_zeta,Sys_latral_poles]=damp(TF_latral);
%% -------------------------------------------------------------------------------------------%%
%% 3DOF Approx. of lateral dynamics (Dutch and Roll Mode)
        %% SS x=[beta, p, r], u=[delta_a, delta_r]
        A_3dof_dutch=A_latral(1:3,1:3);
        B_3dof_dutch=B_latral(1:3,1:2);
        C_3dof_dutch=eye(3);
        D_3dof_dutch=0;
        %% TF
        SS_3dof_dutch = ss(A_3dof_dutch,B_3dof_dutch,C_3dof_dutch,D_3dof_dutch);
        TF_3dof_dutch=tf(SS_3dof_dutch);
        [Sys_3dof_dutch_wn,Sys_3dof_dutch_zeta,Sys_3dof_dutch_poles]=damp(TF_3dof_dutch);
%% -------------------------------------------------------------------------------------------%%
%% 3DOF Approx. of lateral dynamics (Spiral and Roll Mode)        
        %% SS x=[p, r, phi], u=[delta_a, delta_r]
        A_3dof_spiral=A_latral(2:4,2:4);
        B_3dof_spiral=B_latral(2:4,1:2);
        C_3dof_spiral=eye(3);
        D_3dof_spiral=0;
        %% TF
        SS_3dof_spiral = ss(A_3dof_spiral,B_3dof_spiral,C_3dof_spiral,D_3dof_spiral);
        TF_3dof_spiral=tf(SS_3dof_spiral);
        [Sys_3dof_spiral_wn,Sys_3dof_spiral_zeta,Sys_3dof_spiral_poles]=damp(TF_3dof_spiral);
%% -------------------------------------------------------------------------------------------%%
%% 2DOF Approx. of lateral dynamics (Dutch and Roll Mode)                
        %% SS x=[beta, r], u=[delta_a, delta_r]
        A_2dof_dutch=[A_latral(1,1), A_latral(1,3); A_latral(3,1), A_latral(3,3)];
        B_2dof_dutch=[B_latral(1,1:2); B_latral(3,1:2)];
        C_2dof_dutch=eye(2);
        D_2dof_dutch=0;
        %% TF
        SS_2dof_dutch = ss(A_2dof_dutch,B_2dof_dutch,C_2dof_dutch,D_2dof_dutch);
        TF_2dof_dutch=tf(SS_2dof_dutch);
        [Sys_2dof_dutch_wn,Sys_2dof_dutch_zeta,Sys_2dof_dutch_poles]=damp(TF_2dof_dutch);
%% -------------------------------------------------------------------------------------------%%
%% 2DOF Approx. of lateral dynamics (Roll Mode)                
        %% SS x=[p, phi], u=[delta_a, delta_r]
        A_2dof_roll=[A_latral(2,2), A_latral(2,4); A_latral(4,2), A_latral(4,4)];
        B_2dof_roll=[B_latral(2,1:2); B_latral(4,1:2)];
        C_2dof_roll=eye(2);
        D_2dof_roll=0;
        %% TF
        SS_2dof_roll = ss(A_2dof_roll,B_2dof_roll,C_2dof_roll,D_2dof_roll);
        TF_2dof_roll=tf(SS_2dof_roll);
        [Sys_2dof_roll_wn,Sys_2dof_roll_zeta,Sys_2dof_roll_poles]=damp(TF_2dof_roll);
%% -------------------------------------------------------------------------------------------%%
%% 1DOF Approx. of lateral dynamics   
        %% SS x=[p], u=[delta_a, delta_r]
        A_1dof=A_latral(2,2);
        B_1dof=[B_latral(2,1:2)];
        C_1dof=eye(1);
        D_1dof=0;
        %% TF
        SS_1dof = ss(A_1dof,B_1dof,C_1dof,D_1dof);
        TF_1dof=tf(SS_1dof);
        [Sys_1dof_wn,Sys_1dof_zeta,Sys_1dof_poles]=damp(TF_1dof);
%% -------------------------------------------------------------------------------------------%%
%% beta/de TF
        %% exact
        beta_de_TF_exact=tf([TF_latral.num{1}],[TF_latral.den{1}]);
        %% 3DOF(Dutch and Roll Mode)
        beta_de_TF_3dof_dutch=tf([TF_3dof_dutch.num{1}],[TF_3dof_dutch.den{1}]);
        %% 2DOF(Dutch and Roll Mode)
        beta_de_TF_2dof_dutch=tf([TF_2dof_dutch.num{1}],[TF_2dof_dutch.den{1}]);
%% -------------------------------------------------------------------------------------------%%
%% p/de TF
        %% exact
        p_de_TF_exact=tf([TF_latral.num{2}],[TF_latral.den{2}]);
        %% 3DOF (Dutch and Roll Mode)
        p_de_TF_3dof_dutch=tf([TF_3dof_dutch.num{2}],[TF_3dof_dutch.den{2}]);
        %% 3DOF (Spiral and Roll Mode)
        p_de_TF_3dof_spiral=tf([TF_3dof_spiral.num{1}],[TF_3dof_spiral.den{1}]);
        %% 2DOF (Roll Mode)
        p_de_TF_2dof_roll=tf([TF_2dof_roll.num{1}],[TF_2dof_roll.den{1}]);
        %% 1DOF
        p_de_TF_1dof=tf([TF_1dof.num{1}],[TF_1dof.den{1}]);        
%% -------------------------------------------------------------------------------------------%%
%% r/de TF
        %% exact
        r_de_TF_exact=tf([TF_latral.num{3}],[TF_latral.den{3}]);
        %% 3DOF (Dutch and Roll Mode)
        r_de_TF_3dof_dutch=tf([TF_3dof_dutch.num{3}],[TF_3dof_dutch.den{3}]);
        %% 3DOF (Spiral and Roll Mode)        
        r_de_TF_3dof_spiral=tf([TF_3dof_spiral.num{2}],[TF_3dof_spiral.den{2}]);
        %% 2DOF (Dutch and Roll Mode) 
        r_de_TF_2dof_dutch=tf([TF_2dof_dutch.num{2}],[TF_2dof_dutch.den{2}]);
%% -------------------------------------------------------------------------------------------%%
%% phi/de TF
        %% exact
        phi_de_TF_exact=tf([TF_latral.num{4}],[TF_latral.den{4}]);
        %% 3DOF (Spiral and Roll Mode)        
        phi_de_TF_3dof_spiral=tf([TF_3dof_spiral.num{3}],[TF_3dof_spiral.den{3}]);
        %% 2DOF (Roll Mode)
        phi_de_TF_2dof_roll=tf([TF_2dof_roll.num{2}],[TF_2dof_roll.den{2}]);
%% -------------------------------------------------------------------------------------------%%
%% psi/de TF
        %% exact
        psi_de_TF_exact=tf([TF_latral.num{5}],[TF_latral.den{5}]);       
%% -------------------------------------------------------------------------------------------%%
%% beta/dT TF
        %% exact
        beta_dT_TF_exact=tf([TF_latral.num{6}],[TF_latral.den{6}]);
        %% 3DOF(Dutch and Roll Mode)
        beta_dT_TF_3dof_dutch=tf([TF_3dof_dutch.num{4}],[TF_3dof_dutch.den{4}]);
        %% 2DOF(Dutch and Roll Mode)
        beta_dT_TF_2dof_dutch=tf([TF_2dof_dutch.num{3}],[TF_2dof_dutch.den{3}]);
%% -------------------------------------------------------------------------------------------%%
%% p/dT TF
        %% exact
        p_dT_TF_exact=tf([TF_latral.num{7}],[TF_latral.den{7}]);
        %% 3DOF (Dutch and Roll Mode)
        p_dT_TF_3dof_dutch=tf([TF_3dof_dutch.num{5}],[TF_3dof_dutch.den{5}]);
        %% 3DOF (Spiral and Roll Mode)
        p_dT_TF_3dof_spiral=tf([TF_3dof_spiral.num{4}],[TF_3dof_spiral.den{4}]);
        %% 2DOF (Roll Mode)
        p_dT_TF_2dof_roll=tf([TF_2dof_roll.num{3}],[TF_2dof_roll.den{3}]);
        %% 1DOF
        p_dT_TF_1dof=tf([TF_1dof.num{2}],[TF_1dof.den{2}]);        
%% -------------------------------------------------------------------------------------------%%
%% r/dT TF
        %% exact
        r_dT_TF_exact=tf([TF_latral.num{8}],[TF_latral.den{8}]);
        %% 3DOF (Dutch and Roll Mode)
        r_dT_TF_3dof_dutch=tf([TF_3dof_dutch.num{6}],[TF_3dof_dutch.den{6}]);
        %% 3DOF (Spiral and Roll Mode)        
        r_dT_TF_3dof_spiral=tf([TF_3dof_spiral.num{4}],[TF_3dof_spiral.den{4}]);
        %% 2DOF (Dutch and Roll Mode) 
        r_dT_TF_2dof_dutch=tf([TF_2dof_dutch.num{4}],[TF_2dof_dutch.den{4}]);        
%% -------------------------------------------------------------------------------------------%%
%% phi/dT TF
        %% exact
        phi_dT_TF_exact=tf([TF_latral.num{9}],[TF_latral.den{9}]);
        %% 3DOF (Spiral and Roll Mode)        
        phi_dT_TF_3dof_spiral=tf([TF_3dof_spiral.num{6}],[TF_3dof_spiral.den{6}]);
        %% 2DOF (Roll Mode)
        phi_dT_TF_2dof_roll=tf([TF_2dof_roll.num{4}],[TF_2dof_roll.den{4}]);
%% -------------------------------------------------------------------------------------------%%
%% psi/dT TF
        %% exact
        psi_dT_TF_exact=tf([TF_latral.num{10}],[TF_latral.den{10}]);
%% -------------------------------------------------------------------------------------------%%
%% Pole Placement (Aker)
        %% Longitudinal dynamics
                %% State (u) with delta_T
                KC_u = acker(A_long,B_long(:,2),requird_poles_u);
                AC_long_u=A_long-B_long(:,2)*KC_u;
                BC_long_u=B_long(:,2);
                CC_long_u=C_long;
                DC_long_u=D_long;
                X_u=ss(AC_long_u,BC_long_u,CC_long_u,DC_long_u);
                XX_u=tf(X_u);
                [Y_u,T_u]=step(XX_u(1));
                DC_gain_u=1/Y_u(end);
                [y_u,t_u]=step(DC_gain_u*XX_u);
%% -------------------------------------------------------------------------------------------%%
%% LQR
        %% Longitudinal dynamics
                %% State (u) with delta_T
                [kC_u_opt, P_u_opt, poles_u_opt] = lqr(A_long,B_long(:,2),Q_u,R_u);
                AC_long_u_opt = [(A_long-B_long(:,2)*kC_u_opt)];
                BC_long_u_opt = [B_long(:,2)];
                CC_long_u_opt = [C_long];
                DC_long_u_opt = [D_long];
                X_u_opt = ss(AC_long_u_opt,-BC_long_u_opt,CC_long_u_opt,DC_long_u_opt);
                XX_u_opt=tf(X_u_opt);
                [Y_u_opt,T_u_opt]=step(XX_u_opt(1));
                DC_gain_u_opt=1/Y_u_opt(end);
                [y_u_opt,t_u_opt]=step(DC_gain_u_opt*XX_u_opt);
                %% State (theta) with delta_e
                [kC_theta_opt, P_theta_opt, poles_theta_opt] = lqr(A_long,B_long(:,1),Q_theta,R_theta);
                AC_long_theta_opt = [(A_long-B_long(:,1)*kC_theta_opt)];
                BC_long_theta_opt = [B_long(:,1)];
                CC_long_theta_opt = [C_long];
                DC_long_theta_opt = [D_long];
                X_theta_opt = ss(AC_long_theta_opt,-BC_long_theta_opt,CC_long_theta_opt,DC_long_theta_opt);
                XX_theta_opt=tf(X_theta_opt);
                [Y_theta_opt,T_theta_opt]=step(XX_theta_opt(4));
                DC_gain_theta_opt=1/Y_theta_opt(end);
                [y_theta_opt,t_theta_opt]=step(DC_gain_theta_opt*XX_theta_opt);
                %% State (w) with delta_e
                [kC_w_opt, P_w_opt, poles_w_opt] = lqr(A_long,B_long(:,1),Q_w,R_w);
                AC_long_w_opt = [(A_long-B_long(:,1)*kC_w_opt)];
                BC_long_w_opt = [B_long(:,1)];
                CC_long_w_opt = [C_long];
                DC_long_w_opt = [D_long];
                X_w_opt = ss(AC_long_w_opt,-BC_long_w_opt,CC_long_w_opt,DC_long_w_opt);
                XX_w_opt=tf(X_w_opt);
                [Y_w_opt,T_w_opt]=step(XX_w_opt(2));
                DC_gain_w_opt=1/Y_w_opt(end);
                [y_w_opt,t_w_opt]=step(DC_gain_w_opt*XX_w_opt);
        %% Latral Mode
                %% State (psi) with delta_r
                [KC_psi, P_psi, poles_psi] = lqr(A_latral,B_latral(:,2),Q_psi,R_psi);
                AC_psi = [(A_latral-B_latral(:,2)*KC_psi)];
                BC_psi = [B_latral(:,2)];
                CC_psi = [C_latral];
                DC_psi = [D_latral];
                X_psi = ss(AC_psi,-BC_psi,CC_psi,DC_psi);
                XX_psi=tf(X_psi);
                [Y_psi,T_psi]=step(XX_psi(5));
                DC_gain_psi=1/Y_psi(end);
                [y_psi,t_psi]=step(DC_gain_psi*XX_psi,T_psi(end));
%% -------------------------------------------------------------------------------------------%%
%% Plotting 
if Plot ==1
        %% Longitudinal dynamics
                %% u/de
                set(0,'defaultfigureposition',[0 50 1700 630]);
                figure('Name','Open Loop Response of State u','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(u_de_TF_exact);
                legend('Exact')
                title('u/\delta_e Step Response')
                subplot(2,1,2)
                step(u_de_TF_long_period);
                legend('Long Period')
                title('u/\delta_e Step Response')
                %% u/dT
                figure('Name','Open Loop Response of State u','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(u_dT_TF_exact);
                legend('Exact')
                title('u/\delta_T Step Response')
                subplot(2,1,2)
                step(u_dT_TF_long_period);
                legend('Long Period')
                title('u/\delta_T Step Response')
                %% w/de
                figure('Name','Open Loop Response of State w','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(w_de_TF_exact);
                legend('Exact')
                title('w/\delta_e Step Response')
                subplot(2,1,2)
                step(w_de_TF_short_period);
                legend('Short Period')
                title('w/\delta_e Step Response')
                %% w/dT
                figure('Name','Open Loop Response of State w','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(w_dT_TF_exact);
                legend('Exact')
                title('w/\delta_T Step Response')
                subplot(2,1,2)
                step(w_dT_TF_short_period);
                legend('Short Period')
                title('w/\delta_T Step Response')
                %% q/de
                figure('Name','Open Loop Response of State q','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(q_de_TF_exact);
                legend('Exact')
                title('q/\delta_e Step Response')
                subplot(2,1,2)
                step(q_de_TF_short_period);
                legend('Short Period')
                title('q/\delta_e Step Response')
                %% q/dT
                figure('Name','Open Loop Response of State q','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(q_dT_TF_exact);
                legend('Exact')
                title('q/\delta_T Step Response')
                subplot(2,1,2)
                step(q_dT_TF_short_period);
                legend('Short Period')
                title('q/\delta_T Step Response')
                %% theta/de
                figure('Name','Open Loop Response of State Theta','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(theta_de_TF_exact);
                legend('Exact')
                title('\theta/\delta_e Step Response')
                subplot(2,1,2)
                step(theta_de_TF_long_period);
                legend('Long Period')
                title('\theta/\delta_e Step Response')
                %% theta/dT
                figure('Name','Open Loop Response of State Theta','NumberTitle','off')
                set(gcf,'color','w');
                subplot(2,1,1)
                step(theta_dT_TF_exact);
                legend('Exact')
                title('\theta/\delta_T Step Response')
                subplot(2,1,2)
                step(theta_dT_TF_long_period);
                legend('Long Period')
                title('\theta/\delta_T Step Response')
        %% Lateral dynamics
                %% beta/de
                figure('Name','Open Loop Response of State Beta','NumberTitle','off')
                set(gcf,'color','w');
                subplot(3,1,1)
                step(beta_de_TF_exact);
                legend('Exact')
                title('\beta/\delta_e Step Response')
                subplot(3,1,2)
                step(beta_de_TF_3dof_dutch);
                legend('3DOF Dutch')
                title('\beta/\delta_e Step Response')
                subplot(3,1,3)
                step(beta_de_TF_2dof_dutch);
                legend('2DOF Dutch')
                title('\beta/\delta_e Step Response')
                %% beta/dT
                figure('Name','Open Loop Response of State Beta','NumberTitle','off')
                set(gcf,'color','w');
                subplot(3,1,1)
                step(beta_dT_TF_exact);
                legend('Exact')
                title('\beta/\delta_T Step Response')
                subplot(3,1,2)
                step(beta_dT_TF_3dof_dutch);
                legend('3DOF Dutch')
                title('\beta/\delta_T Step Response')
                subplot(3,1,3)
                step(beta_dT_TF_2dof_dutch);
                legend('2DOF Dutch')
                title('\beta/\delta_T Step Response')
                %% p/de
                figure('Name','Open Loop Response of State p','NumberTitle','off')
                set(gcf,'color','w');
                subplot(5,1,1)
                step(p_de_TF_exact);
                legend('Exact')
                title('p/\delta_e Step Response')
                subplot(5,1,2)
                step(p_de_TF_3dof_dutch);
                title(' ')
                legend('3DOF Dutch')
                subplot(5,1,3)
                step(p_de_TF_3dof_spiral);
                title(' ')
                legend('3DOF Spiral')
                subplot(5,1,4)
                step(p_de_TF_2dof_roll);
                title(' ')
                legend('2DOF Roll')
                subplot(5,1,5)
                step(p_de_TF_1dof);
                title(' ')
                legend('1DOF')
                %% p/dT
                figure('Name','Open Loop Response of State p','NumberTitle','off')
                set(gcf,'color','w');
                subplot(5,1,1)
                step(p_dT_TF_exact);
                legend('Exact')
                title('p/\delta_T Step Response')
                subplot(5,1,2)
                step(p_dT_TF_3dof_dutch);
                legend('3DOF Dutch')
                title(' ')
                subplot(5,1,3)
                step(p_dT_TF_3dof_spiral);
                legend('3DOF Spiral')
                title(' ')
                subplot(5,1,4)
                step(p_dT_TF_2dof_roll);
                legend('2DOF Roll')
                title(' ')
                subplot(5,1,5)
                step(p_dT_TF_1dof);
                legend('1DOF')
                title(' ')
                %% r/de
                figure('Name','Open Loop Response of State r','NumberTitle','off')
                set(gcf,'color','w');
                subplot(4,1,1)
                step(r_de_TF_exact);
                legend('Exact')
                title('r/\delta_e Step Response')
                subplot(4,1,2)
                step(r_de_TF_3dof_dutch);
                legend('3DOF Dutch')
                title(' ')
                subplot(4,1,3)
                step(r_de_TF_3dof_spiral);
                legend('3DOF Spiral')
                title(' ')
                subplot(4,1,4)
                step(r_de_TF_2dof_dutch);
                legend('2DOF Dutch')
                title(' ')
                %% r/dT
                figure('Name','Open Loop Response of State r','NumberTitle','off')
                set(gcf,'color','w');
                subplot(4,1,1)
                step(r_dT_TF_exact);
                legend('Exact')
                title('r/\delta_T Step Response')
                subplot(4,1,2)
                step(r_dT_TF_3dof_dutch);
                legend('3DOF Dutch')
                title(' ')
                subplot(4,1,3)
                step(r_dT_TF_3dof_spiral);
                legend('3DOF Spiral')
                title(' ')
                subplot(4,1,4)
                step(r_dT_TF_2dof_dutch);
                legend('2DOF Dutch')
                title(' ')
                %% phi/de
                figure('Name','Open Loop Response of State Phi','NumberTitle','off')
                set(gcf,'color','w');
                subplot(3,1,1)
                step(phi_de_TF_exact);
                legend('Exact')
                title('\phi/\delta_e Step Response')
                subplot(3,1,2)
                step(phi_de_TF_3dof_spiral);
                legend('3DOF Spiral')
                title('\phi/\delta_e Step Response')
                subplot(3,1,3)
                step(phi_de_TF_2dof_roll);
                legend('2DOF Roll')
                title('\phi/\delta_e Step Response')                
                %% phi/dT
                figure('Name','Open Loop Response of State Phi','NumberTitle','off')
                set(gcf,'color','w');
                subplot(3,1,1)
                step(phi_dT_TF_exact);
                legend('Exact')
                title('\phi/\delta_T Step Response')
                subplot(3,1,2)
                step(phi_dT_TF_3dof_spiral);
                legend('3DOF Spiral')
                title('\phi/\delta_T Step Response')
                subplot(3,1,3)
                step(phi_dT_TF_2dof_roll);
                legend('2DOF Roll')
                title('\phi/\delta_T Step Response')                
                %% psi/de
                figure('Name','Open Loop Response of State Psi','NumberTitle','off')
                hold all;
                set(gcf,'color','w');
                step(psi_de_TF_exact);
                legend('Exact')
                title('\psi/\delta_e Step Response')                             
                %% psi/dT
                figure('Name','Open Loop Response of State Psi','NumberTitle','off')
                hold all;
                set(gcf,'color','w');
                step(psi_dT_TF_exact);
                legend('Exact')
                title('\psi/\delta_T Step Response')
        %% Pole Placement (Aker)
                %% Longitudinal dynamics
                        %% State (u) with delta_T 
                        figure('Name','Control on State u Using Pole Placement (Aker)','NumberTitle','off')
                        subplot(4,2,1)
                        plot(t_u,y_u(:,1),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('u(t)','fontsize',18)
                        title('Step Response of u','fontsize',18)
                        subplot(4,2,2)
                        plot(t_u,-KC_u(1)*y_u(:,1),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        title('Control Action','fontsize',18)
                        subplot(4,2,3)
                        plot(t_u,y_u(:,2),'linewidth',2)
                        grid on;
                        ylabel('w(t)','fontsize',18)
                        subplot(4,2,4)
                        plot(t_u,-KC_u(2)*y_u(:,2),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        subplot(4,2,5)
                        plot(t_u,y_u(:,3),'linewidth',2)
                        grid on;
                        ylabel('q(t)','fontsize',18)
                        subplot(4,2,6)
                        plot(t_u,-KC_u(3)*y_u(:,3),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        subplot(4,2,7)
                        plot(t_u,y_u(:,4),'linewidth',2)
                        grid on;
                        ylabel('\theta(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        subplot(4,2,8)
                        plot(t_u,-KC_u(4)*y_u(:,4),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
        %% LQR
                %% Longitudinal dynamics
                        %% State (u) with delta_T
                        figure('Name','Control on State u Using LQR','NumberTitle','off')
                        subplot(4,2,1)
                        plot(t_u_opt,y_u_opt(:,1),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('u(t)','fontsize',18)
                        title('Step Response of u','fontsize',18)
                        subplot(4,2,2)
                        plot(t_u_opt,-kC_u_opt(1)*y_u_opt(:,1),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        title('Control Action','fontsize',18)
                        subplot(4,2,3)
                        plot(t_u_opt,y_u_opt(:,2),'linewidth',2)
                        grid on;
                        ylabel('w(t)','fontsize',18)
                        subplot(4,2,4)
                        plot(t_u_opt,-kC_u_opt(2)*y_u_opt(:,2),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        subplot(4,2,5)
                        plot(t_u_opt,y_u_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('q(t)','fontsize',18)
                        subplot(4,2,6)
                        plot(t_u_opt,-kC_u_opt(3)*y_u_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        subplot(4,2,7)
                        plot(t_u_opt,y_u_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\theta(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        subplot(4,2,8)
                        plot(t_u_opt,-kC_u_opt(4)*y_u_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\delta_T(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        %% State (theta) with delta_e
                        figure('Name','Control on State theta Using LQR','NumberTitle','off')
                        subplot(4,2,1)
                        plot(t_theta_opt,y_theta_opt(:,1),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('u(t)','fontsize',18)
                        title('Step Response of \theta','fontsize',18)
                        subplot(4,2,2)
                        plot(t_theta_opt,-kC_theta_opt(1)*y_theta_opt(:,1),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        title('Control Action','fontsize',18)
                        subplot(4,2,3)
                        plot(t_theta_opt,y_theta_opt(:,2),'linewidth',2)
                        grid on;
                        ylabel('w(t)','fontsize',18)
                        subplot(4,2,4)
                        plot(t_theta_opt,-kC_theta_opt(2)*y_theta_opt(:,2),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        subplot(4,2,5)
                        plot(t_theta_opt,y_theta_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('q(t)','fontsize',18)
                        subplot(4,2,6)
                        plot(t_theta_opt,-kC_theta_opt(3)*y_theta_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        subplot(4,2,7)
                        plot(t_theta_opt,y_theta_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\theta(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        subplot(4,2,8)
                        plot(t_theta_opt,-kC_theta_opt(4)*y_theta_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        %% State (w) with delta_e
                        figure('Name','Control on State w Using LQR','NumberTitle','off')
                        subplot(4,2,1)
                        plot(t_w_opt,y_w_opt(:,1),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('u(t)','fontsize',18)
                        title('Step Response of w','fontsize',18)
                        subplot(4,2,2)
                        plot(t_w_opt,-kC_w_opt(1)*y_w_opt(:,1),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        title('Control Action','fontsize',18)
                        subplot(4,2,3)
                        plot(t_w_opt,y_w_opt(:,2),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('w(t)','fontsize',18)
                        subplot(4,2,4)
                        plot(t_w_opt,-kC_w_opt(2)*y_w_opt(:,2),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        subplot(4,2,5)
                        plot(t_w_opt,y_w_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('q(t)','fontsize',18)
                        subplot(4,2,6)
                        plot(t_w_opt,-kC_w_opt(3)*y_w_opt(:,3),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        subplot(4,2,7)
                        plot(t_w_opt,y_w_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\theta(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        subplot(4,2,8)
                        plot(t_w_opt,-kC_w_opt(4)*y_w_opt(:,4),'linewidth',2)
                        grid on;
                        ylabel('\delta_e(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                %% Latral Mode
                        %% State (psi) with delta_r
                        figure('Name','Control on State psi Using LQR','NumberTitle','off')
                        subplot(5,2,1)
                        plot(t_psi,y_psi(:,1),'linewidth',2)
                        set(gcf,'color','w');
                        grid on;
                        ylabel('\beta(t)','fontsize',18)
                        title('Step Response of \psi','fontsize',18)
                        subplot(5,2,2)
                        plot(t_psi,-KC_psi(1)*y_psi(:,1),'linewidth',2)
                        grid on;
                        ylabel('\delta_r(t)','fontsize',18)
                        title('Control Action','fontsize',18)
                        subplot(5,2,3)
                        plot(t_psi,y_psi(:,2),'linewidth',2)
                        grid on;
                        ylabel('p(t)','fontsize',18)
                        subplot(5,2,4)
                        plot(t_psi,-KC_psi(2)*y_psi(:,2),'linewidth',2)
                        grid on;
                        ylabel('\delta_r(t)','fontsize',18)
                        subplot(5,2,5)
                        plot(t_psi,y_psi(:,3),'linewidth',2)
                        grid on;
                        ylabel('r(t)','fontsize',18)
                        subplot(5,2,6)
                        plot(t_psi,-KC_psi(3)*y_psi(:,3),'linewidth',2)
                        grid on;
                        ylabel('\delta_r(t)','fontsize',18)
                        subplot(5,2,7)
                        plot(t_psi,y_psi(:,4),'linewidth',2)
                        grid on;
                        ylabel('\phi(t)','fontsize',18)
                        subplot(5,2,8)
                        plot(t_psi,-KC_psi(4)*y_psi(:,4),'linewidth',2)
                        grid on;
                        ylabel('\delta_r(t)','fontsize',18)                        
                        subplot(5,2,9)
                        plot(t_psi,y_psi(:,5),'linewidth',2)
                        grid on;
                        ylabel('\psi(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
                        subplot(5,2,10)
                        plot(t_psi,-KC_psi(5)*y_psi(:,5),'linewidth',2)
                        grid on;
                        ylabel('\delta_r(t)','fontsize',18)
                        xlabel('Time (sec)','fontsize',18)
end                        
%% Save figures
if Save ==1 && Plot ==1
    for S=1:23
        figure(S);
        saveas(gcf, [num2str(S) '.png']);
    end
end