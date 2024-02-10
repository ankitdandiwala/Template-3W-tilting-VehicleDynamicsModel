%{The purpose of this code is to create a template for simulating vehicle
% dynamics of three-wheeled tilting vehicle including lateral dynamics, 
% multibody dynamics, and tire model. The mathematical equations are 
% capturing the dynamics of the veicle with many assumtions. However, 
% this template can be used to convert it into a sophisticated system-level
% simulation model for the narrow tilting three-wheeled vehicle, referring
% journal articles using these links: 1) https://doi.org/10.1177/14644193231207496
% or 2) https://doi.org/10.1016/j.eswa.2023.121970.
%}

%Vehicle parameters
T = 0.83; %track of the vehicle in meter
mt = 130; %mass of tilting parts in kg
mnt = 159; %mass of non-tilting parts in kg
lf = 1.0226; %longitudinal distance between front axle and CG
lr = 0.5224; %longitudinal distance between rear axle and CG
Izz = 176.91 % Yaw moment of inertia (kgm^2)
htp_t =	0.380; %Height of tilt axis in a yz-plane revolute joint at rear (m)
htp_CG = 0.2515; %Height of tilt axis in a yz-plane passing through CG of the vehicle (m)
htilt = 0.40587; %Height of the CG of the tiltig parts from the c/s of tit axis in a plane pasing through CG of the veicle
hnt = 0.3117; %Height of the CG of the non-tilting parts (m)
l = lf+lr; %wheelbase
m = mt+mnt; %total mass of the vehicle

%Constants
step = 0.001; %stepsize in second
g = 9.81; %gravitational acceleration in m/s^2
t = 0:step:10;

%Initialization and array creation
delta = zeros(length(t)+1,1);
ay = zeros(length(t)+1,1);
theta = zeros(length(t)+1,1);
Fyf = zeros(length(t)+1,1);
FyL = zeros(length(t)+1,1);
FyR = zeros(length(t)+1,1);
Fy = zeros(length(t)+1,1);
Fzf = zeros(length(t)+1,1);
FzL = zeros(length(t)+1,1);
FzR = zeros(length(t)+1,1);
Fz = zeros(length(t)+1,1);
Vy = zeros(length(t)+1,1);


for i = 1:length(t)     
    
    %Steering Input
    if t(i) >= 2 && t(i) <= 8
        delta(i) = 0.2;
    else
        delta(i) = 0;
    end
    %Velocity Input
    Vx = 10;
    
    %Tilt Angle
    ht(i+1) = htp_CG + htilt*cos(theta(i))
    h(i+1) = (mt*ht(i+1)+mnt*hnt)/m % height of the CG (m)
    ycg(i+1) = (mt/m)*htilt*sin(theta(i)); %Lateral displacement of the CG
    theta_CG(i+1) = atan(ycg(i+1)/h(i+1));  %tilt angle of CG (rad)
    % SM(i+1) =1-(2*h(i+1)*(l*abs(ay(i+1))/(T*g*lf);
    % SMthr = -0.00000000006*SM^4 + 0.0000000003*SM^3 - 0.05*SM^2 + 0.515*SM + 0.535
    % DTA = sign(ay(i))*asin((m*T*lf*(SMthr-SM(i)))/(2*mt*l*htilt));
    % if abs(DTA) > 0.4186
    %     SMthr =  0.006*SM^4 - 0.0127*SM^3 + 0.007*SM^2 + 0.9987*SM + 0.2702
    % end
    theta(i+1) = atan(ay(i)/g);


    %Lateral dynamics
    ay(i+1) = Fy(i)/m;
    if i < length(t)-1
        Vy(i+1) = (ay(i+1)-ay(i))/(step);
    end
    %{Include yaw dyamics as given in 1) https://doi.org/10.1177/14644193231207496
    % or 2) https://doi.org/10.1016/j.eswa.2023.121970. Solve equtions of 
    % yaw and lateral dynamics to determin the SLIP ANGLES at front and rear%}

    %Front Tire Model (Pure slip)
    C_alpha(i+1)=9.74*Fzf(i);
    C_gamma(i+1)=0.86*Fzfr(i);
    d4=1.2;
    d6=0.1;
    d7=0.15;
    d8=1.6;
    C=d8;
    D(i+1)=d4*Fzf(i)/(1+d7*theta(i)^2);
    B(i+1)=C_alpha(i+1)/C/D(i+1);
    Shf(i+1)=C_gamma(i+1)*theta(i+1)/C_alpha(i+1);
    Sv(i+1)=d6*Fzfr*theta(i+1);
    Sh(i+1)=Shf(i+1)-(Sv(i+1)/C_alpha(i+1));
    Fyf(i+1) = D(i+1)*sin(C*atan(B*(alpha_f(i+1)+Sh(i+1))))+Sv(i+1);

    %Rear Tire Model (pure slip)
    c1=8;
    c2=1.33;
    C=1.3;
    E=-1;
    muzero=1;
    Fz0=3000; %(Newton)
    C_alpha(i+1)=c1*c2*Fz0*sin(2*atan(Fzr/Fz0));
    D0=muzero*Fz0;
    B0(i+1)=C_alpha(i+1)/C/D0;
    alpha_eq(i+1)=Fz0/Fzr(i)*alpha_r(i+1);
    x(i+1)=tan(alpha_eq(i+1));
    Fy0(i+1)=D0*sin(1.3*atan(B0(i+1)*x(i+1)+(B0(i+1)*x(i+1)-atan(B0(i+1)*x(i+1)))));
    Fyr(i+1) = Fzr(i)/Fz0*Fy0(i+1);

    % Fyf(i+1) = 0.1*Fzf(i);
    % FyL(i+1) = 0.1*FzL(i);
    % FyR(i+1) = 0.1*FzR(i);
    % Fyr(i+1) = FyR(i+1) + FyL(i+1);
    % Fy(i+1) = Fyf(i+1)+FyL(i+1)+FyR(i+1);

    %Multibody Dynamics
    Fzf(i+1) = m*g*(lr/l);
    FzL(i+1) = m*ay(i+1)*(-h/T)+m*g*(lf/(2*l)+h*sin(theta(i+1))/T);
    FzR(i+1) = m*ay(i+1)*(h/T)+m*g*(lf/(2*l)-h*sin(theta(i+1))/T);
    Fz(i+1) = Fzf(i+1)+FzL(i+1)+FzR(i+1);
end
result = table(delta,ay,Fy,Fz,theta,Vy);

