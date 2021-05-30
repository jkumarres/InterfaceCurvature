function [a,b,r] = FitCircle(x,y)

N = length(x);
NN = 10;
eta = 0.001;

n = floor(N/2);

x0 = sum(x)/N;
y0 = sum(y)/N;

x = x - x0;
y = y - y0;

x1 = x(1);
y1 = y(1);
x2 = x(n);
y2 = y(n);
x3 = x(N);
y3 = y(N);
    
x12 = x1-x2;
x23 = x2-x3;
x31 = x3-x1;
    
y12 = y1-y2;
y23 = y2-y3;
y31 = y3-y1;
    
Num = (x12^2 + y12^2) * (x23^2 + y23^2) * (x31^2 + y31^2);
Den = (x1 * y23 + x2 * y31 + x3 * y12);
    
r = abs( 0.5 * sqrt(Num) / Den );

Num = x1^2 * y23 + x2^2 * y31 + x3^2 * y12 - y12 * y23 * y31;
a = 0.5 * Num / Den;

Num = y1^2 * x23 + y2^2 * x31 + y3^2 * x12 - x12 * x23 * x31;
b = 0.5 * Num / (-Den);

s = r;

r = r/s;
x = x/s;
y = y/s;
a = a/s;
b = b/s;

for i=1:NN
    E = 0;
    dEda = 0;
    dEdb = 0;
    dEdr = 0;
    
    for j=1:N
        dE = (x(j)-a)^2 + (y(j)-b)^2 - r^2 ;
        E = E + dE^2;
        
        dEda = dEda - 4 * dE * (x(j)-a);
        dEdb = dEdb - 4 * dE * (y(j)-b);
        dEdr = dEdr - 4 * dE * r;
    end
    
    r = r - eta * dEdr;
    a = a - eta * dEda;
    b = b - eta * dEdb;    
end

a = s * a;
b = s * b;
r = s * r;

a = a + x0;
b = b + y0;

return;