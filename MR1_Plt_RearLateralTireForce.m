clear;
k=0.0; % For Pure Slip
%k=0.1 %For Combineds slip
sigma_x=k/(1+k);
%///////////A free-trajectory, Veneri (2020) Veh Sys Dyn
Fz0=3500;%Newton
pKy1=44.2;
pKy2=2.5977;
pDy1=1.8217;
pDy2=-0.4388;
lambda_uy=0.84;
pCy1=1.733;
pEy1=0.29446;

%%///////////The effect of Veneri (2021) Veh Sys Dyn
% Fz0=809;%Newton
% pKx1=39.06;
% pKx3=-0.23;
% pDx1=-1.20;
% pDx2=0.71;
% lambda_ux=1.00;
% pCx1=2.31;
% pEx1=1.00;

NormalizedFyr1000=[];
NormalizedFyr2000=[];
NormalizedFyr3000=[];
NormalizedFyr4000=[];
Fyr1000=[];
Fyr2000=[];
Fyr3000=[];
Fyr4000=[];
Sigma_y1000=[];
Sigma_y2000=[];
Sigma_y3000=[];
Sigma_y4000=[];

for sigma_y=-0.5:0.001:0.5
    Fzr=1000;
    %Fz0=8/5*Fzr;
    dfzr1000=(Fzr-Fz0)/Fz0;
    Ey=pEy1;
    Cy=pCy1;
    Dy=(pDy1+pDy2*dfzr1000)*lambda_uy;
    Ky=Fz0*pKy1*sin(2*atan(Fzr/pKy2/Fz0));
    By=Ky/Cy/Dy/Fzr;
    
    %k=sigma_x/(1-sigma_x);
    %sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fyr = Fzr*(sigma_y/sigma)*Dy*sin(Cy*atan(By*sigma-Ey*(By*sigma-atan(By*sigma))));
    NormFyr=Fyr/Fzr;
    
    NormalizedFyr1000=[NormalizedFyr1000;NormFyr];
    Fyr1000=[Fyr1000;Fyr];
    Sigma_y1000=[Sigma_y1000;sigma_y];
end

for sigma_y=-0.5:0.001:0.5
    Fzr=2000;
    %Fz0=8/5*Fzr;
    dfzr2000=(Fzr-Fz0)/Fz0;
    Ey=pEy1;
    Cy=pCy1;
    Dy=(pDy1+pDy2*dfzr2000)*lambda_uy;
    Ky=Fz0*pKy1*sin(2*atan(Fzr/pKy2/Fz0));
    By=Ky/Cy/Dy/Fzr;

    %k=sigma_x/(1-sigma_x);
    %sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fyr = Fzr*(sigma_y/sigma)*Dy*sin(Cy*atan(By*sigma-Ey*(By*sigma-atan(By*sigma))));
    NormFyr=Fyr/Fzr;
    
    NormalizedFyr2000=[NormalizedFyr2000;NormFyr];
    Fyr2000=[Fyr2000;Fyr];
    Sigma_y2000=[Sigma_y2000;sigma_y];
end

for sigma_y=-0.5:0.001:0.5
    Fzr=3000;
    %Fz0=8/5*Fzr;
    dfzr3000=(Fzr-Fz0)/Fz0;
    Ey=pEy1;
    Cy=pCy1;
    Dy=(pDy1+pDy2*dfzr3000)*lambda_uy;
    Ky=Fz0*pKy1*sin(2*atan(Fzr/pKy2/Fz0));
    By=Ky/Cy/Dy/Fzr;

    %k=sigma_x/(1-sigma_x);
    %sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fyr = Fzr*(sigma_y/sigma)*Dy*sin(Cy*atan(By*sigma-Ey*(By*sigma-atan(By*sigma))));
    NormFyr=Fyr/Fzr;
    
    NormalizedFyr3000=[NormalizedFyr3000;NormFyr];
    Fyr3000=[Fyr3000;Fyr];
    Sigma_y3000=[Sigma_y3000;sigma_y];
end

for sigma_y=-0.5:0.001:0.5
    Fzr=4000;
    %Fz0=8/5*Fzr;
    dfzr4000=(Fzr-Fz0)/Fz0;
    Ey=pEy1;
    Cy=pCy1;
    Dy=(pDy1+pDy2*dfzr4000)*lambda_uy;
    Ky=Fz0*pKy1*sin(2*atan(Fzr/pKy2/Fz0));
    By=Ky/Cy/Dy/Fzr;

    %k=sigma_x/(1-sigma_x);
    %sigma_y=tan(alpha_r)/(1+k);
    sigma=sqrt(sigma_x^2+sigma_y^2);
    Fyr = Fzr*(sigma_y/sigma)*Dy*sin(Cy*atan(By*sigma-Ey*(By*sigma-atan(By*sigma))));
    NormFyr=Fyr/Fzr;
    
    NormalizedFyr4000=[NormalizedFyr4000;NormFyr];
    Fyr4000=[Fyr4000;Fyr];
    Sigma_y4000=[Sigma_y4000;sigma_y];
end

disp(["dfzr1000",dfzr1000]);
disp(["dfzr2000",dfzr2000]);
disp(["dfzr3000",dfzr3000]);
disp(["dfzr4000",dfzr4000]);

% hold off;
% plot(Sigma_y1000, NormalizedFyr1000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'r'         );
% hold on
% plot(Sigma_y2000, NormalizedFyr2000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'g'         );
% hold on
% plot(Sigma_y3000, NormalizedFyr3000, ...
%     'LineWidth'       , 1           , ...
%     'LineStyle'       , '-'        , ...
%     'Color'           , 'b'         );
% hold on
% plot(Sigma_y4000, NormalizedFyr4000, ...
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
% XLabel = xlabel('\sigma_y'         );
% YLabel = ylabel('F_{yr} (N)'         );
% set(gca,'Box'         , 'off'     , ...
%     'TickDir'     , 'out'     , ...
%     'TickLength'  , [0.020 0.020] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'YGrid'       , 'on'      ,...
%   'XGrid'       , 'on'  );
% % yticks([-50 0 50 ])



hold off;
plot(Sigma_y1000, Fyr1000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '-'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_y2000, Fyr2000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '--'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_y3000, Fyr3000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , '-.'        , ...
    'Color'           , 'k'         );
hold on
plot(Sigma_y4000, Fyr4000, ...
    'LineWidth'       , 1           , ...
    'LineStyle'       , ':'        , ...
    'Color'           , 'k'         );
h=legend('F_{zr}=1000N', 'F_{zr}=2000N','F_{zr}=3000N','F_{zr}=4000N');
% %axis equal;

% set(gcf, 'PaperPositionMode', 'auto');
% print trajectory_all.jpg '-bestfit'
% close;

%Title  = title ('Lateral tire force');
XLabel = xlabel('\sigma_y'         );
YLabel = ylabel('F_{yr} (N)'         );
set(gca,'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [0.020 0.020] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      ,...
  'XGrid'       , 'off'  );
% yticks([-50 0 50 ])

