%% 1.2 
j = 10; o = 15; n = 14; a = 1; t = 20; h = 8; v = 22; y = 25; g = 7; i = 9;
jonathanvanhyning = j+o+n+a+t+h+a+n+v+a+n+h+y+n+i+n+g

%% 1.3.1
% z = 25-(100-7exp(5+cos(pi/3)
z = 25-(100-7*exp(5+cos(pi/3)));

%% 1.4.1 
syms t s
a = log10(cos(t)-3*s);
diff(a, t)
diff(a, s)
%% 1.4.2 
asin(1)
%% 1.6
size = 4; 
Fibonacci = zeros(size,size);
Fibonacci(1) = 1;
Fibonacci(2) = 1;
for i = 2:1:(size^2)-1
    Fibonacci(i+1) = Fibonacci(i) + Fibonacci(i-1);
end
Fibonacci
