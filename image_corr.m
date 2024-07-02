%% Image correction algorithm

%The code starts out from observed y- and z- position and moves back
%to true y- and z- positions by using refraction rules. You can find full
%derivation in:

%Caldag, H. O., Acemoglu, A., Yesilyurt, S. (2017). Experimental
%characterization of helical swimming trajectories in circular channels.
%Microfluidics and Nanofluidics, 21:136.

%Inputs:
%xi: Coordinate in horizontal axis (y- coordinate for our case)
%yi: Coordinate in vertical axis (z- coordinate for our case)
%R: Total circular channel radius (including its thickness)
%n: Refractivity index

%Outputs:
%x0: Corrected coordinate on horizontal axis.
%y0: Corrected coordinate on vertical axis.
%ro: Corrected radial position.

function [xo,yo,ro] = image_corr(xi,yi,R,n) 

eps = 0.001;
xa = xi;
yb = yi;
tha = acos(xi/R);
ya = R*sin(tha);
th2a = abs(pi/2-tha);
th1a = asin(sin(th2a)/n);   
alpha = th2a-th1a;
yc = (sin(th1a)./sin(alpha+eps)*R);
thb = acos(yi/R);
xb = R*sin(thb);
th2b = abs(pi/2-thb);
th1b = asin(sin(th2b)/n);
beta = th2b-th1b;
xc = (sin(th1b)./sin(beta+eps)*R);
sa = (ya+yc)./xa;
sb = yb./(-xb-xc);
bb = -sb.*xc;
ba = -yc;
xo = (bb-ba)./(sa-sb);
yo = sb.*xo + bb;
ro = sqrt(xo.*xo + yo.*yo);
end
