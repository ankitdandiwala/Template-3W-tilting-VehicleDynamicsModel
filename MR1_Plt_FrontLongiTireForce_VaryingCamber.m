
%This code is for Longitudinal Force (With Both: PureSlip & CombineSlip), From:PACEJKA, 2006
clear;

%Constants
Fz0=4000;%Newton

%Constants for PureSlip
pDx1=1.210;
pDx2=-0.037;
pCx1=1.685;
pEx1=0.344;
pEx2=0.095;
pEx3=-0.02;
pEx4=0;
qsy1=0.01;
%qsx1=0;
pKx1=21.51;
pKx2=-0.163;
pKx3=0.245;
pVx1=0;
pVx2=0;
lamda_Fz0=1;
lambda_ux=1;
lambda_muv=0;
%R0=;
lambda_cx=1;
lambda_Ex=1;
lambda_My=1;
lambda_Kxk=1;
epsilon_vx=1;
lambda_vx=1;
A_mu=10;
epsilon_x=0;
term_Vcx=1;

%Constant for Combined Slip
c1=8;
c2=1.33;
c5=1;
rBx1=12.35;
rBx2=-10.77;
rBx3=0;
rCx1=1.092;
rHx1=0.007;
Fzf0c=3000; %Newton

%Vehicle's output: Input to model
%VariableName ending with 'A' means: Force vary, camber constant (zero)
gamma1000A=0;
gamma2000A=0;
gamma3000A=0;
gamma4000A=0;
Fzf1000A=1000; %Newton
Fzf2000A=2000; %Newton
Fzf3000A=3000; %Newton
Fzf4000A=4000; %Newton
%VariableName ending with 'B' means: Camber vary, Force constant (3000 N)
gamma1000B=deg2rad(-5);
gamma2000B=deg2rad(0);
gamma3000B=deg2rad(5);
gamma4000B=deg2rad(10);
Fzf1000B=3000; %Newton
Fzf2000B=3000; %Newton
Fzf3000B=3000; %Newton
Fzf4000B=3000; %Newton

%Creating Numeric Matrix for plot
NormalizedFxf01000=[];
NormalizedFxf02000=[];
NormalizedFxf03000=[];
NormalizedFxf04000=[];
Fxf01000=[];
Fxf02000=[];
Fxf03000=[];
Fxf04000=[];
k1000=[];
k2000=[];
k3000=[];
k4000=[];
Fxf1000=[];
Fxf2000=[];
Fxf3000=[];
Fxf4000=[];

for k=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    alpha=0.5;
    %Fzf=Fzf1000A;
    gamma=gamma1000B;
    Fzf=Fzf1000B;
    %gamma=gamma1000B;
    
    %Model (Pure Slip):
    dfzf1000=(Fzf-Fz0)/Fz0;
    Kxk=Fzf*(pKx1+pKx2*dfzf1000)*2.71^(pKx3*dfzf1000);
    Svx=Fzf*(pVx1+pVx2*dfzf1000)*term_Vcx;
    Shx=-(qsy1*Fzf+Svx)/Kxk;
    mu_x=(pDx1+pDx2*dfzf1000);
    Ex=(pEx1+pEx2*dfzf1000+pEx3*dfzf1000^2);
    Cx=pCx1;
    Dx=mu_x+Fzf;
    kx=k+Shx;
    Bx=Kxk/(Cx*Dx+epsilon_x);
    Fxf0 = Dx*sin(Cx*atan(Bx*kx-Ex*(Bx*kx-atan(Bx*kx))))+Svx;
    Fxf01000=[Fxf01000;Fxf0];
    k1000=[k1000;k];
    
    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0c*sin(2*atan(Fzf/c2/Fzf0c));
    Cfgamma=c5*Fzf;
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Shxalpha=rHx1;
    Cxalpha=rCx1;
    Bxalpha=(rBx1+rBx3*gamma^2)*cos(atan(rBx2*k));
    alphas=alpha_star*+Shxalpha;
    Gxalpha0=cos(Cxalpha*atan(Bxalpha*Shxalpha));
    Gxalpha=cos(Cxalpha*atan(Bxalpha*alphas))/Gxalpha0;
    Fxf=Gxalpha*Fxf0;
    Fxf1000=[Fxf1000;Fxf];
end


for k=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    alpha=0.5;
    %Fzf=Fzf2000A;
    gamma=gamma2000B;
    Fzf=Fzf2000B;
    %gamma=gamma2000B;
    
    %Model (Pure Slip):
    dfzf2000=(Fzf-Fz0)/Fz0;
    Kxk=Fzf*(pKx1+pKx2*dfzf2000)*2.71^(pKx3*dfzf2000);
    Svx=Fzf*(pVx1+pVx2*dfzf2000)*term_Vcx;
    Shx=-(qsy1*Fzf+Svx)/Kxk;
    mu_x=(pDx1+pDx2*dfzf2000);
    Ex=(pEx1+pEx2*dfzf2000+pEx3*dfzf2000^2);
    Cx=pCx1;
    Dx=mu_x+Fzf;
    kx=k+Shx;
    Bx=Kxk/(Cx*Dx+epsilon_x);
    Fxf0 = Dx*sin(Cx*atan(Bx*kx-Ex*(Bx*kx-atan(Bx*kx))))+Svx;
    Fxf02000=[Fxf02000;Fxf0];
    k2000=[k2000;k];
    
    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0c*sin(2*atan(Fzf/c2/Fzf0c));
    Cfgamma=c5*Fzf;
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Shxalpha=rHx1;
    Cxalpha=rCx1;
    Bxalpha=(rBx1+rBx3*gamma^2)*cos(atan(rBx2*k));
    alphas=alpha_star*+Shxalpha;
    Gxalpha0=cos(Cxalpha*atan(Bxalpha*Shxalpha));
    Gxalpha=cos(Cxalpha*atan(Bxalpha*alphas))/Gxalpha0;
    Fxf=Gxalpha*Fxf0;
    Fxf2000=[Fxf2000;Fxf];
end


for k=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    alpha=0.5;
    %Fzf=Fzf3000A;
    gamma=gamma3000B;
    Fzf=Fzf3000B;
    %gamma=gamma3000B;
    
    %Model (Pure Slip):
    dfzf3000=(Fzf-Fz0)/Fz0;
    Kxk=Fzf*(pKx1+pKx2*dfzf3000)*2.71^(pKx3*dfzf3000);
    Svx=Fzf*(pVx1+pVx2*dfzf3000)*term_Vcx;
    Shx=-(qsy1*Fzf+Svx)/Kxk;
    mu_x=(pDx1+pDx2*dfzf3000);
    Ex=(pEx1+pEx2*dfzf3000+pEx3*dfzf3000^2);
    Cx=pCx1;
    Dx=mu_x+Fzf;
    kx=k+Shx;
    Bx=Kxk/(Cx*Dx+epsilon_x);
    Fxf0 = Dx*sin(Cx*atan(Bx*kx-Ex*(Bx*kx-atan(Bx*kx))))+Svx;
    Fxf03000=[Fxf03000;Fxf0];
    k3000=[k3000;k];
    
    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0c*sin(2*atan(Fzf/c2/Fzf0c));
    Cfgamma=c5*Fzf;
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Shxalpha=rHx1;
    Cxalpha=rCx1;
    Bxalpha=(rBx1+rBx3*gamma^2)*cos(atan(rBx2*k));
    alphas=alpha_star*+Shxalpha;
    Gxalpha0=cos(Cxalpha*atan(Bxalpha*Shxalpha));
    Gxalpha=cos(Cxalpha*atan(Bxalpha*alphas))/Gxalpha0;
    Fxf=Gxalpha*Fxf0;
    Fxf3000=[Fxf3000;Fxf];
end

for k=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    alpha=0.5;
    %Fzf=Fzf4000A;
    gamma=gamma4000B;
    Fzf=Fzf4000B;
    %gamma=gamma4000B;
    
    %Model (Pure Slip):
    dfzf4000=(Fzf-Fz0)/Fz0;
    Kxk=Fzf*(pKx1+pKx2*dfzf4000)*2.71^(pKx3*dfzf4000);
    Svx=Fzf*(pVx1+pVx2*dfzf4000)*term_Vcx;
    Shx=-(qsy1*Fzf+Svx)/Kxk;
    mu_x=(pDx1+pDx2*dfzf4000);
    Ex=(pEx1+pEx2*dfzf4000+pEx3*dfzf4000^2);
    Cx=pCx1;
    Dx=mu_x+Fzf;
    kx=k+Shx;
    Bx=Kxk/(Cx*Dx+epsilon_x);
    Fxf0 = Dx*sin(Cx*atan(Bx*kx-Ex*(Bx*kx-atan(Bx*kx))))+Svx;
    Fxf04000=[Fxf04000;Fxf0];
    k4000=[k4000;k];
    
    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0c*sin(2*atan(Fzf/c2/Fzf0c));
    Cfgamma=c5*Fzf;
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Shxalpha=rHx1;
    Cxalpha=rCx1;
    Bxalpha=(rBx1+rBx3*gamma^2)*cos(atan(rBx2*k));
    alphas=alpha_star*+Shxalpha;
    Gxalpha0=cos(Cxalpha*atan(Bxalpha*Shxalpha));
    Gxalpha=cos(Cxalpha*atan(Bxalpha*alphas))/Gxalpha0;
    Fxf=Gxalpha*Fxf0;
    Fxf4000=[Fxf4000;Fxf];
end

%%Plot for Pure Slip
% hold off;
% plot(k1000, Fxf01000,'LineWidth',1,'LineStyle','-','Color','k');
% hold on
% plot(k2000, Fxf02000,'LineWidth',1,'LineStyle','--','Color','k');
% hold on
% plot(k3000, Fxf03000,'LineWidth',1,'LineStyle','-.','Color','k');
% hold on
% plot(k4000, Fxf04000,'LineWidth',1,'LineStyle',':','Color','k');
% h=legend('\gamma=-5\circ', '\gamma=0\circ','\gamma=5\circ','\gamma=10\circ');
% % %axis equal;
% 
% % set(gcf, 'PaperPositionMode', 'auto');
% % print trajectory_all.jpg '-bestfit'
% % close;
% 
% %Title  = title ('Lateral tire force');
% XLabel = xlabel('k'         );
% YLabel = ylabel('F_{xf_0} (N)'         );
% set(gca,'Box'         , 'off'     , ...                                                                                                                 
%     'TickDir'     , 'out'     , ...
%     'TickLength'  , [0.020 0.020] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'off'      ,...
%   'XGrid'       , 'off'  );
% % yticks([-50 0 50 ])


%Plot for Combined Slip
hold off;
plot(k1000, Fxf1000,'LineWidth',1,'LineStyle','-','Color','k');
hold on
plot(k2000, Fxf2000,'LineWidth',1,'LineStyle','--','Color','k');
hold on
plot(k3000, Fxf3000,'LineWidth',1,'LineStyle','-.','Color','k');
hold on
plot(k4000, Fxf4000,'LineWidth',1,'LineStyle',':','Color','k');
h=legend('\gamma= -5\circ', '\gamma= 0\circ','\gamma= 5\circ','\gamma= 10\circ');
%h=legend('\gamma=-5\circ', '\gamma=0\circ','\gamma=5\circ','\gamma=10\circ');
% %axis equal;

% set(gcf, 'PaperPositionMode', 'auto');
% print trajectory_all.jpg '-bestfit'
% close;

%Title  = title ('Lateral tire force');
XLabel = xlabel('k');
YLabel = ylabel('F_{xf} (N)');
set(gca,'Box','off', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [0.020 0.020] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      ,...
  'XGrid'       , 'off'  );
% yticks([-50 0 50 ])






