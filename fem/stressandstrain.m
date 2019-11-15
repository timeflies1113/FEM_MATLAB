%% 单元应力、应变
function [stress,strain]=stressandstrain(E,NU,x1,y1,x2,y2,x3,y3,x4,y4,u)
%输入参数：弹性模量、泊松比、四点坐标、位移；输出参数：单元中心点的应力、应变
syms s t;
a = (y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b = (y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c = (x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d = (x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;
B1 = [a*(t-1)/4-b*(s-1)/4 0 ; 0 c*(s-1)/4-d*(t-1)/4 ;
   c*(s-1)/4-d*(t-1)/4 a*(t-1)/4-b*(s-1)/4];
B2 = [a*(1-t)/4-b*(-1-s)/4 0 ; 0 c*(-1-s)/4-d*(1-t)/4 ;
   c*(-1-s)/4-d*(1-t)/4 a*(1-t)/4-b*(-1-s)/4];
B3 = [a*(t+1)/4-b*(s+1)/4 0 ; 0 c*(s+1)/4-d*(t+1)/4 ;
   c*(s+1)/4-d*(t+1)/4 a*(t+1)/4-b*(s+1)/4];
B4 = [a*(-1-t)/4-b*(1-s)/4 0 ; 0 c*(1-s)/4-d*(-1-t)/4 ;
   c*(1-s)/4-d*(-1-t)/4 a*(-1-t)/4-b*(1-s)/4];
Bfirst = [B1 B2 B3 B4];
Jfirst = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ;
   s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8;
B = Bfirst/J;
D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];
strain=B*u;
stress=D*B*u;%应变=B*U，应力=D*应变
% 计算单元中心点的应力应变，代替整个单元的应力应变
strain00= subs(strain, {s,t}, {0,0});
stress00= subs(stress, {s,t}, {0,0});
strain=double(strain00);
stress= double(stress00);
end