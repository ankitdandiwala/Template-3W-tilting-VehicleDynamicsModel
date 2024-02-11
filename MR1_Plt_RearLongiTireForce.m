
%THIS COE IS FOR PURE LONGI (PACEJKA 2006)
clear;

%///////////A free-trajectory, Veneri (2020) Veh Sys Dyn
Fz0=3500;%Newton
pKx1=30.5;
pKx3=0.2766;
pDx1=1.8757;
pDx2=-0.127;
lambda_ux=0.8 %0.93;
pCx1=1.6935;
pEx1=0.07708;
alpha_r=0.1;

%%///////////The effect of Veneri (2021) Veh Sys Dyn
% Fz0=809;%Newton
% pKx1=39.06;
% pKx3=-0.23;
% pDx1=-1.20;
% pDx2=0.71;
% lambda_ux=1.00;
% pCx1=2.31;
% pEx1=1.00;

NormalizedFxr1000=[];
NormalizedFxr2000=[];
NormalizedFxr3000=[];
NormalizedFxr4000=[];
Fxr1000=[];
Fxr2000=[];
Fxr3000=[];
Fxr4000=[];
Sigma_x1000=[];
Sigma_x2000=[];
Sigma_x3000=[];
Sigma_x4000=[];

for sigma_x=-0.5:0.001:0.5
    Fzr=1000;
    %Fz0=8/5*Fzr;
    dfzr1000=(Fzr-Fz0)/Fz0;
    Ex=pEx1;
    Cx=pCx1;
    Dx=(pDx1+pDx2*dfzr1000)*lambda_ux;
    Kx=Fzr*pKx1*2.71^(pKx3*dfzr1000);
    Bx=Kx/Cx/Dx/Fzr;
    %sigma_x=k/(1+k);
    k=sigma_x/(1-sigma_x);
    sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fxr = Fzr*(sigma_x/sigma)*Dx*sin(Cx*atan(Bx*sigma-Ex*(Bx*sigma-atan(Bx*sigma))));
    NormFxr=Fxr/Fzr;
    NormalizedFxr1000=[NormalizedFxr1000;NormFxr];
    Fxr1000=[Fxr1000;Fxr];
    Sigma_x1000=[Sigma_x1000;sigma_x];
end

for sigma_x=-0.5:0.001:0.5
    Fzr=2000;
    %Fz0=1.2*8/5*Fzr;
    dfzr2000=(Fzr-Fz0)/Fz0;
    Ex=pEx1;
    Cx=pCx1;
    Dx=(pDx1+pDx2*dfzr2000)*lambda_ux;
    Kx=Fzr*pKx1*2.71^(pKx3*dfzr2000);
    Bx=Kx/Cx/Dx/Fzr;
    %sigma_x=k/(1+k);
    k=sigma_x/(1-sigma_x);
    sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fxr = Fzr*(sigma_x/sigma)*Dx*sin(Cx*atan(Bx*sigma-Ex*(Bx*sigma-atan(Bx*sigma))));
    NormFxr=Fxr/Fzr;
    NormalizedFxr2000=[NormalizedFxr2000;NormFxr];
    Fxr2000=[Fxr2000;Fxr];
    Sigma_x2000=[Sigma_x2000;sigma_x];
end

for sigma_x=-0.5:0.001:0.5
    Fzr=3000;
    %Fz0=1.2*1.2*8/5*Fzr;
    dfzr3000=(Fzr-Fz0)/Fz0;
    Ex=pEx1;
    Cx=pCx1;
    Dx=(pDx1+pDx2*dfzr3000)*lambda_ux;
    Kx=Fzr*pKx1*2.71^(pKx3*dfzr3000);
    Bx=Kx/Cx/Dx/Fzr;
    %sigma_x=k/(1+k);
    k=sigma_x/(1-sigma_x);
    sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fxr = Fzr*(sigma_x/sigma)*Dx*sin(Cx*atan(Bx*sigma-Ex*(Bx*sigma-atan(Bx*sigma))));
    NormFxr=Fxr/Fzr;
    NormalizedFxr3000=[NormalizedFxr3000;NormFxr];
    Fxr3000=[Fxr3000;Fxr];
    Sigma_x3000=[Sigma_x3000;sigma_x];
end

for sigma_x=-0.5:0.001:0.5
    Fzr=4000;
    %Fz0=x^n*8/5*Fzr;
    dfzr4000=(Fzr-Fz0)/Fz0;
    Ex=pEx1;
    Cx=pCx1;
    Dx=(pDx1+pDx2*dfzr4000)*lambda_ux;
    Kx=Fzr*pKx1*2.71^(pKx3*dfzr4000);
    Bx=Kx/Cx/Dx/Fzr;
    %sigma_x=k/(1+k);
    k=sigma_x/(1-sigma_x);
    sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fxr = Fzr*(sigma_x/sigma)*Dx*sin(Cx*atan(Bx*sigma-Ex*(Bx*sigma-atan(Bx*sigma))));
    NormFxr=Fxr/Fzr;
    NormalizedFxr4000=[NormalizedFxr4000;NormFxr];
    Fxr4000=[Fxr4000;Fxr];
    Sigma_x4000=[Sigma_x4000;sigma_x];
end
% disp(["dfzr1000",dfzr1000]);
% disp(["dfzr2000",dfzr2000]);
% disp(["dfzr3000",dfzr3000]);
% disp(["dfzr4000",dfzr4000]);

% hold off;
% plot(Sigma_x1000, NormalizedFxr1000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'r'         );
% hold on
% plot(Sigma_x2000, NormalizedFxr2000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'g'         );
% hold on
% plot(Sigma_x3000, NormalizedFxr3000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'b'         );
% hold on
% plot(Sigma_x4000, NormalizedFxr4000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'c'         );
% h=legend('F_{zr}=1000N', 'F_{zr}=2000N','F_{zr}=3000N','F_{zr}=4000N');
% % %axis equal;
% 
% % set(gcf, 'PaperPositionMode', 'auto');
% % print trajectory_all.jpg '-bestfit'
% % close;
% 
% %Title  = title ('Lateral tire force');
% XLabel = xlabel('\sigma_x'         );
% YLabel = ylabel('F_{xr} (N)'         );
% set(gca,'Box'         , 'off'     , ...
%     'TickDir'     , 'out'     , ...
%     'TickLength'  , [0.020 0.020] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'on'      ,...
%   'XGrid'       , 'on'  );
% % yticks([-50 0 50 ])



hold off;
plot(Sigma_x1000, Fxr1000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '-'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_x2000, Fxr2000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '--'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_x3000, Fxr3000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '-.'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_x4000, Fxr4000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , ':'        , ...
    'Color'           , 'k'         );
h=legend('F_{zr}=1000N', 'F_{zr}=2000N','F_{zr}=3000N','F_{zr}=4000N');
% %axis equal;

% set(gcf, 'PaperPositionMode', 'auto');
% print trajectory_all.jpg '-bestfit'
% close;

%Title  = title ('Lateral tire force');
XLabel = xlabel('\sigma_x'         );
YLabel = ylabel('F_{xr} (N)'         );
set(gca,'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [0.020 0.020] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      ,...
  'XGrid'       , 'off'  );
% yticks([-50 0 50 ])

