
%This code is for Longitudinal Force (With Both: PureSlip & CombineSlip), From:PACEJKA, 2006
clear;

%Constants
Fz0=4000;%Newton
d4=1.2;
d6=0.1;
d7=0.15;
d8=1.6;
C=d8;

%Constant for Combined Slip
rBy1=6.461;
rBy2=4.196;
rBy3=-0.015;
rBy4=0;
rCy1=1.081;
rHy1=0.009;
rVy1=0.053;
rVy2=-0.073;
rVy3=0.517;
rVy4=35.44;
rVy5=1.9;
rVy6=-10.71;
Mu_y=0.8;
Fzf0=4000;
c1=8;
c2=1.33;
c5=1;


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

%slip ratio
k=0.5;

%Creating Numeric Matrix for plot
NormalizedFyf01000=[];
NormalizedFyf02000=[];
NormalizedFyf03000=[];
NormalizedFyf04000=[];
Fyf01000=[];
Fyf02000=[];
Fyf03000=[];
Fyf04000=[];
alpha1000=[];
alpha2000=[];
alpha3000=[];
alpha4000=[];
Fyf1000=[];
Fyf2000=[];
Fyf3000=[];
Fyf4000=[];

for alpha=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    Fzf=Fzf1000A;
    gamma=gamma1000A;
    %Fzf=Fzf1000B;
    %gamma=gamma1000B;
    dfzf1000=(Fzf-Fzf0)/Fzf0;
    
    %Model (Pure Slip):
    Calpha=9.74*Fzf;
    Cgamma=0.86*Fzf;
    D=d4*Fzf/(1+d7*gamma^2);
    B=Calpha/C/D;
    Shf=Cgamma*gamma/Calpha;
    Sv=d6*Fzf*gamma;
    Sh=Shf-(Sv/Calpha);
    Fyf0 = D*sin(C*atan(B*(alpha+Sh)))+Sv;
    Fyf01000=[Fyf01000;Fyf0];
    alpha1000=[alpha1000;alpha];

    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0*sin(2*atan(Fzf/c2/Fzf0));
    Cfgamma=c5*Fzf;
    disp(["Cfalpha",Cfalpha,"Cfgamma",Cfgamma]);
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Dvyk=Mu_y*Fzf*(rVy1+rVy2*dfzf1000+rVy3*gamma)*cos(atan(rVy4*alpha_star));
    Shyk=rHy1;
    Cyk=rCy1;
    Byk=(rBy1+rBy4*gamma^2)*cos(atan(rBy2*(alpha_star-rBy3)));
    ks=k*+Shyk;
    Svyk=Dvyk*sin(rVy5*atan(rVy6*k));
    Gyk0=cos(Cyk*atan(Byk*Shyk));
    Gyk=cos(Cyk*atan(Byk*ks))/Gyk0;
    disp(["Gyk",Gyk,"Svyk",Svyk]);
    Fyf=Gyk*Fyf0+Svyk;
    Fyf1000=[Fyf1000;Fyf];
end

for alpha=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    Fzf=Fzf2000A;
    gamma=gamma2000A;
    %Fzf=Fzf1000B;
    %gamma=gamma1000B;
    dfzf2000=(Fzf-Fzf0)/Fzf0;
    
    %Model (Pure Slip):
    Calpha=9.74*Fzf;
    Cgamma=0.86*Fzf;
    D=d4*Fzf/(1+d7*gamma^2);
    B=Calpha/C/D;
    Shf=Cgamma*gamma/Calpha;
    Sv=d6*Fzf*gamma;
    Sh=Shf-(Sv/Calpha);
    Fyf0 = D*sin(C*atan(B*(alpha+Sh)))+Sv;
    Fyf02000=[Fyf02000;Fyf0];
    alpha2000=[alpha2000;alpha];

    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0*sin(2*atan(Fzf/c2/Fzf0));
    Cfgamma=c5*Fzf;
    disp(["Cfalpha",Cfalpha,"Cfgamma",Cfgamma]);
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Dvyk=Mu_y*Fzf*(rVy1+rVy2*dfzf2000+rVy3*gamma)*cos(atan(rVy4*alpha_star));
    Shyk=rHy1;
    Cyk=rCy1;
    Byk=(rBy1+rBy4*gamma^2)*cos(atan(rBy2*(alpha_star-rBy3)));
    ks=k*+Shyk;
    Svyk=Dvyk*sin(rVy5*atan(rVy6*k));
    Gyk0=cos(Cyk*atan(Byk*Shyk));
    Gyk=cos(Cyk*atan(Byk*ks))/Gyk0;
    disp(["Gyk",Gyk,"Svyk",Svyk]);
    Fyf=Gyk*Fyf0+Svyk;
    Fyf2000=[Fyf2000;Fyf];
end

for alpha=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    Fzf=Fzf3000A;
    gamma=gamma3000A;
    %Fzf=Fzf1000B;
    %gamma=gamma1000B;
    dfzf3000=(Fzf-Fzf0)/Fzf0;
    
    %Model (Pure Slip):
    Calpha=9.74*Fzf;
    Cgamma=0.86*Fzf;
    D=d4*Fzf/(1+d7*gamma^2);
    B=Calpha/C/D;
    Shf=Cgamma*gamma/Calpha;
    Sv=d6*Fzf*gamma;
    Sh=Shf-(Sv/Calpha);
    Fyf0 = D*sin(C*atan(B*(alpha+Sh)))+Sv;
    Fyf03000=[Fyf03000;Fyf0];
    alpha3000=[alpha3000;alpha];

    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0*sin(2*atan(Fzf/c2/Fzf0));
    Cfgamma=c5*Fzf;
    disp(["Cfalpha",Cfalpha,"Cfgamma",Cfgamma]);
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Dvyk=Mu_y*Fzf*(rVy1+rVy2*dfzf3000+rVy3*gamma)*cos(atan(rVy4*alpha_star));
    Shyk=rHy1;
    Cyk=rCy1;
    Byk=(rBy1+rBy4*gamma^2)*cos(atan(rBy2*(alpha_star-rBy3)));
    ks=k*+Shyk;
    Svyk=Dvyk*sin(rVy5*atan(rVy6*k));
    Gyk0=cos(Cyk*atan(Byk*Shyk));
    Gyk=cos(Cyk*atan(Byk*ks))/Gyk0;
    disp(["Gyk",Gyk,"Svyk",Svyk]);
    Fyf=Gyk*Fyf0+Svyk;
    Fyf3000=[Fyf3000;Fyf];
end

for alpha=-0.5:0.001:0.5
    %Vehicle's output: Input to the Model
    Fzf=Fzf4000A;
    gamma=gamma4000A;
    %Fzf=Fzf1000B;
    %gamma=gamma1000B;
    dfzf4000=(Fzf-Fzf0)/Fzf0;
    
    %Model (Pure Slip):
    Calpha=9.74*Fzf;
    Cgamma=0.86*Fzf;
    D=d4*Fzf/(1+d7*gamma^2);
    B=Calpha/C/D;
    Shf=Cgamma*gamma/Calpha;
    Sv=d6*Fzf*gamma;
    Sh=Shf-(Sv/Calpha);
    Fyf0 = D*sin(C*atan(B*(alpha+Sh)))+Sv;
    Fyf04000=[Fyf04000;Fyf0];
    alpha4000=[alpha4000;alpha];

    %Model (Combined Slip):
    Cfalpha=c1*c2*Fzf0*sin(2*atan(Fzf/c2/Fzf0));
    Cfgamma=c5*Fzf;
    disp(["Cfalpha",Cfalpha,"Cfgamma",Cfgamma]);
    alpha_star=alpha+(Cfgamma/Cfalpha)*gamma;
    Dvyk=Mu_y*Fzf*(rVy1+rVy2*dfzf4000+rVy3*gamma)*cos(atan(rVy4*alpha_star));
    Shyk=rHy1;
    Cyk=rCy1;
    Byk=(rBy1+rBy4*gamma^2)*cos(atan(rBy2*(alpha_star-rBy3)));
    ks=k*+Shyk;
    Svyk=Dvyk*sin(rVy5*atan(rVy6*k));
    Gyk0=cos(Cyk*atan(Byk*Shyk));
    Gyk=cos(Cyk*atan(Byk*ks))/Gyk0;
    disp(["Gyk",Gyk,"Svyk",Svyk]);
    Fyf=Gyk*Fyf0+Svyk;
    Fyf4000=[Fyf4000;Fyf];
end

% %%Plots for Pure Slip
% hold off;
% plot(alpha1000, Fyf01000,'LineWidth',1,'LineStyle','-','Color','k');
% hold on
% plot(alpha3000, Fyf02000,'LineWidth',1,'LineStyle','--','Color','k');
% hold on
% plot(alpha3000, Fyf03000,'LineWidth',1,'LineStyle','-.','Color','k');
% hold on
% plot(alpha4000, Fyf04000,'LineWidth',1,'LineStyle',':','Color','k');
% XLabel = xlabel('\alpha_f');

% Plots for Pure Slip
hold off;
plot(alpha1000, Fyf1000,'LineWidth',1,'LineStyle','-','Color','k');
hold on
plot(alpha3000, Fyf2000,'LineWidth',1,'LineStyle','--','Color','k');
hold on
plot(alpha3000, Fyf3000,'LineWidth',1,'LineStyle','-.','Color','k');
hold on
plot(alpha4000, Fyf4000,'LineWidth',1,'LineStyle',':','Color','k');
XLabel = xlabel('\alpha_f');


h=legend('F_{zf}=1000N', 'F_{zf}=2000N','F_{zf}=3000N','F_{zf}=4000N');
YLabel = ylabel('F_{yf} (N)');
set(gca,'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [0.020 0.020] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      ,...
  'XGrid'       , 'off'  );
% yticks([-50 0 50 ])



