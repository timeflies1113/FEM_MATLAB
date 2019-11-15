%% 单元刚度矩阵计算方法一
function k=Stiffness(E,NU,h,x1,y1,x2,y2,x3,y3,x4,y4,p)
%BilinearQuadElementStiffness   This function returns the element
%                               stiffness matrix for a bilinear
%                               quadrilateral element with modulus
%                               of elasticity E, Poisson's ratio
%                               NU, thickness h, coordinates of
%                               node 1 (x1,y1), coordinates
%                               of node 2 (x2,y2), coordinates of
%                               node 3 (x3,y3), and coordinates of
%                               node 4 (x4,y4). Use p = 1 for cases
%                               of plane stress, and p = 2 for
%                               cases of plane strain.
%                               The size of the element
%                               stiffness matrix is 8 x 8.
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
if p == 1
    D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];
elseif p == 2
    D = (E/(1+NU)/(1-2*NU))*[1-NU, NU, 0 ; NU, 1-NU, 0 ; 0, 0, (1-2*NU)/2];
end
BD = J*transpose(B)*D*B;
r = int(int(BD, t, -1, 1), s, -1, 1);
z = h*r;
w = double(z);
end

%% 单元刚度矩阵方法二
function w = BilinearQuadElementStiffness2(E,NU,h,x1,y1,x2,y2,x3,y3,x4,y4,p)
%BilinearQuadElementStiffness   This function returns the element 
%                               stiffness matrix for a bilinear   
%                               quadrilateral element with modulus 
%                               of elasticity E, Poisson's ratio 
%                               NU, thickness h, coordinates of 
%                               node 1 (x1,y1), coordinates 
%                               of node 2 (x2,y2), coordinates of 
%                               node 3 (x3,y3), and coordinates of 
%                               node 4 (x4,y4). Use p = 1 for cases 
%                               of plane stress, and p = 2 for 
%                               cases of plane strain.
%                               The size of the element 
%                               stiffness matrix is 8 x 8.
syms s t;
N1 = (1-s)*(1-t)/4;
N2 = (1+s)*(1-t)/4;
N3 = (1+s)*(1+t)/4;
N4 = (1-s)*(1+t)/4;
x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;
xs = diff(x,s);
xt = diff(x,t);
ys = diff(y,s);
yt = diff(y,t);
J = xs*yt - ys*xt;
N1s = diff(N1,s);
N2s = diff(N2,s);
N3s = diff(N3,s);
N4s = diff(N4,s);
N1t = diff(N1,t);
N2t = diff(N2,t);
N3t = diff(N3,t);
N4t = diff(N4,t);
B11 = yt*N1s - ys*N1t;
B12 = 0;
B13 = yt*N2s - ys*N2t;
B14 = 0;
B15 = yt*N3s - ys*N3t;
B16 = 0;
B17 = yt*N4s - ys*N4t;
B18 = 0;
B21 = 0;
B22 = xs*N1t - xt*N1s;
B23 = 0;
B24 = xs*N2t - xt*N2s;
B25 = 0;
B26 = xs*N3t - xt*N3s;
B27 = 0;
B28 = xs*N4t - xt*N4s;
B31 = xs*N1t - xt*N1s;
B32 = yt*N1s - ys*N1t;
B33 = xs*N2t - xt*N2s;
B34 = yt*N2s - ys*N2t;
B35 = xs*N3t - xt*N3s;
B36 = yt*N3s - ys*N3t;
B37 = xs*N4t - xt*N4s;
B38 = yt*N4s - ys*N4t;
B = [B11 B12 B13 B14 B15 B16 B17 B18 ;
   B21 B22 B23 B24 B25 B26 B27 B28 ;
   B31 B32 B33 B34 B35 B36 B37 B38];
if p == 1 
   D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];
elseif p == 2
   D = (E/(1+NU)/(1-2*NU))*[1-NU, NU, 0 ; NU, 1-NU, 0 ; 0, 0, (1-2*NU)/2];
end
BD = transpose(B)*D*B/J;
r = int(int(BD, t, -1, 1), s, -1, 1);
z = h*r;
w = double(z);
end

%% 刚度矩阵组装
function y = BilinearQuadAssemble(K,k,i,j,m,n)
    %BilinearQuadAssemble   This function assembles the element
    %                       stiffness matrix k of the bilinear
    %                       quadrilateral element with nodes i, j,
    %                       m, and n into the global stiffness
    %                       matrix K.
    %                       This function returns the global stiffness
    %                       matrix K after the element stiffness matrix
    %                       k is assembled.
    K(2*i-1,2*i-1) = K(2*i-1,2*i-1) + k(1,1);
    K(2*i-1,2*i) = K(2*i-1,2*i) + k(1,2);
    K(2*i-1,2*j-1) = K(2*i-1,2*j-1) + k(1,3);
    K(2*i-1,2*j) = K(2*i-1,2*j) + k(1,4);
    K(2*i-1,2*m-1) = K(2*i-1,2*m-1) + k(1,5);
    K(2*i-1,2*m) = K(2*i-1,2*m) + k(1,6);
    K(2*i-1,2*n-1) = K(2*i-1,2*n-1) + k(1,7);
    K(2*i-1,2*n) = K(2*i-1,2*n) + k(1,8);
    K(2*i,2*i-1) = K(2*i,2*i-1) + k(2,1);
    K(2*i,2*i) = K(2*i,2*i) + k(2,2);
    K(2*i,2*j-1) = K(2*i,2*j-1) + k(2,3);
    K(2*i,2*j) = K(2*i,2*j) + k(2,4);
    K(2*i,2*m-1) = K(2*i,2*m-1) + k(2,5);
    K(2*i,2*m) = K(2*i,2*m) + k(2,6);
    K(2*i,2*n-1) = K(2*i,2*n-1) + k(2,7);
    K(2*i,2*n) = K(2*i,2*n) + k(2,8);
    K(2*j-1,2*i-1) = K(2*j-1,2*i-1) + k(3,1);
    K(2*j-1,2*i) = K(2*j-1,2*i) + k(3,2);
    K(2*j-1,2*j-1) = K(2*j-1,2*j-1) + k(3,3);
    K(2*j-1,2*j) = K(2*j-1,2*j) + k(3,4);
    K(2*j-1,2*m-1) = K(2*j-1,2*m-1) + k(3,5);
    K(2*j-1,2*m) = K(2*j-1,2*m) + k(3,6);
    K(2*j-1,2*n-1) = K(2*j-1,2*n-1) + k(3,7);
    K(2*j-1,2*n) = K(2*j-1,2*n) + k(3,8);
    K(2*j,2*i-1) = K(2*j,2*i-1) + k(4,1);
    K(2*j,2*i) = K(2*j,2*i) + k(4,2);
    K(2*j,2*j-1) = K(2*j,2*j-1) + k(4,3);
    K(2*j,2*j) = K(2*j,2*j) + k(4,4);
    K(2*j,2*m-1) = K(2*j,2*m-1) + k(4,5);
    K(2*j,2*m) = K(2*j,2*m) + k(4,6);
    K(2*j,2*n-1) = K(2*j,2*n-1) + k(4,7);
    K(2*j,2*n) = K(2*j,2*n) + k(4,8);
    K(2*m-1,2*i-1) = K(2*m-1,2*i-1) + k(5,1);
    K(2*m-1,2*i) = K(2*m-1,2*i) + k(5,2);
    K(2*m-1,2*j-1) = K(2*m-1,2*j-1) + k(5,3);
    K(2*m-1,2*j) = K(2*m-1,2*j) + k(5,4);
    K(2*m-1,2*m-1) = K(2*m-1,2*m-1) + k(5,5);
    K(2*m-1,2*m) = K(2*m-1,2*m) + k(5,6);
    K(2*m-1,2*n-1) = K(2*m-1,2*n-1) + k(5,7);
    K(2*m-1,2*n) = K(2*m-1,2*n) + k(5,8);
    K(2*m,2*i-1) = K(2*m,2*i-1) + k(6,1);
    K(2*m,2*i) = K(2*m,2*i) + k(6,2);
    K(2*m,2*j-1) = K(2*m,2*j-1) + k(6,3);
    K(2*m,2*j) = K(2*m,2*j) + k(6,4);
    K(2*m,2*m-1) = K(2*m,2*m-1) + k(6,5);
    K(2*m,2*m) = K(2*m,2*m) + k(6,6);
    K(2*m,2*n-1) = K(2*m,2*n-1) + k(6,7);
    K(2*m,2*n) = K(2*m,2*n) + k(6,8);
    K(2*n-1,2*i-1) = K(2*n-1,2*i-1) + k(7,1);
    K(2*n-1,2*i) = K(2*n-1,2*i) + k(7,2);
    K(2*n-1,2*j-1) = K(2*n-1,2*j-1) + k(7,3);
    K(2*n-1,2*j) = K(2*n-1,2*j) + k(7,4);
    K(2*n-1,2*m-1) = K(2*n-1,2*m-1) + k(7,5);
    K(2*n-1,2*m) = K(2*n-1,2*m) + k(7,6);
    K(2*n-1,2*n-1) = K(2*n-1,2*n-1) + k(7,7);
    K(2*n-1,2*n) = K(2*n-1,2*n) + k(7,8);
    K(2*n,2*i-1) = K(2*n,2*i-1) + k(8,1);
    K(2*n,2*i) = K(2*n,2*i) + k(8,2);
    K(2*n,2*j-1) = K(2*n,2*j-1) + k(8,3);
    K(2*n,2*j) = K(2*n,2*j) + k(8,4);
    K(2*n,2*m-1) = K(2*n,2*m-1) + k(8,5);
    K(2*n,2*m) = K(2*n,2*m) + k(8,6);
    K(2*n,2*n-1) = K(2*n,2*n-1) + k(8,7);
    K(2*n,2*n) = K(2*n,2*n) + k(8,8);
    y = K;
end

%% 三角形面积
function y = BilinearQuadElementArea(x1,y1,x2,y2,x3,y3,x4,y4)
%BilinearQuadElementArea   This function returns the area 
%                          of the bilinear quadrilateral 
%                          element whose first node has 
%                          coordinates (x1,y1), second 
%                          node has coordinates (x2,y2), 
%                          third node has coordinates 
%                          (x3,y3), and fourth node has
%                          coordinates (x4,y4) .
yfirst = (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;
ysecond = (x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3))/2;
y = yfirst + ysecond;
end

%% 单元主应力
function y = BilinearQuadElementPStresses(sigma)
%BilinearQuadElementPStresses   This function returns the element 
%                               principal stresses and their 
%                               angle given the element 
%                               stress vector.
R = (sigma(1) + sigma(2))/2;
Q = ((sigma(1) - sigma(2))/2)^2 + sigma(3)*sigma(3);
M = 2*sigma(3)/(sigma(1) - sigma(2));
s1 = R + sqrt(Q);
s2 = R - sqrt(Q);
theta = (atan(M)/2)*180/pi;
y = [s1 ; s2 ; theta];
end

%% 单元应力（单元力矢量？）
function w = BilinearQuadElementStresses(E,NU,x1,y1,x2,y2,x3,y3,x4,y4,p,u)
%BilinearQuadElementStresses   This function returns the element 
%                              stress vector for a bilinear   
%                              quadrilateral element with modulus 
%                              of elasticity E, Poisson's ratio 
%                              NU, coordinates of 
%                              node 1 (x1,y1), coordinates 
%                              of node 2 (x2,y2), coordinates of 
%                              node 3 (x3,y3), and coordinates of 
%                              node 4 (x4,y4). Use p = 1 for cases 
%                              of plane stress, and p = 2 for 
%                              cases of plane strain.
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
if p == 1 
   D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];
elseif p == 2
   D = (E/(1+NU)/(1-2*NU))*[1-NU, NU, 0 ; NU, 1-NU, 0 ; 0, 0, (1-2*NU)/2];
end
w = D*B*u
%
% We also calculate the stresses at the centroid of the element
%
wcent = subs(w, {s,t}, {0,0});
w = double(wcent);
end