clc; clear all ;

% Meshing
[T,R] = meshgrid(linspace(0,2*pi,200),linspace(0,2,16)); 
X = R.*cos(T);  % Circular X
Y = R.*sin(T);  % Circular Y
z = ones(size(X));


% Time Vector
stept = 0.1;
limt = 10;
mint = 1;
t = mint:stept:limt; 

% Wave Properties
A = 1;              % Amplitude
w = 5;              % Frequency

%%
Z = cell(size(t));
K = cell(size(X));
H = cell(size(X));

for k = 1:(length(t)-1)
        Theta = atan(abs(Y)./X);
        Theta(1,:) = 0;
        Z{k} = A*(abs(Y).*(sin(Theta).^2)).*sin(w*Theta+t(k)).^2;
        [K{k},H{k}] = surfature(X,Y,Z{k});
end
    

%%
figure(1) % Animated Plot of Motion

for k = 1:(length(t)-1)
    surf(X,Y,Z{k})
    zlim([-1 5])
    drawnow
    pause(0.2)
end
%%
figure(2)
surf(X,Y,Z{2})
figure(3)
surf(X,Y,K{2})