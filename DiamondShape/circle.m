function h = circle(xc,yc,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + xc;
yunit = r * sin(th) + yc;
h = fill(xunit, yunit,'w');