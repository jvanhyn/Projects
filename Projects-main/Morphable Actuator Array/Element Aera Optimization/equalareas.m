
function F = equalareas(M)

h1 = M(1)
h2 = M(2)
h3 = M(3)
h4 = M(4)
h5 = M(5)
A  = M(6)

theta = 3;
w = 4.51;
s = 1.5;

A_0 = 0.5*w^2/tand(theta);

L_0 = 2*A_0/w;
L_1 = s + h1 + L_0;
L_2 = 2*s + h1 + h2 + L_0;
L_3 = 3*s + h1 + h2 + h3 + L_0;
L_4 = 3*s + h1 + h2 + h3 + + h4 + L_0;

L_a = h1 + L_0;
L_b = h2 + L_1;
L_c = h3 + L_2;
L_d = h4 + L_3;
L_e = h5 + L_4;

F(1) = 0.5*L_a^2*tand(theta) - A_0 - A
F(2) = 0.5*L_b^2*tand(theta) - 0.5*L_1^2*tand(theta) - A
F(3) = 0.5*L_c^2*tand(theta) - 0.5*L_2^2*tand(theta) - A
F(4) = 0.5*L_d^2*tand(theta) - 0.5*L_3^2*tand(theta) - A
F(5) = 0.5*L_e^2*tand(theta) - 0.5*L_4^2*tand(theta) - A
F(6) = 4*s + h1 + h2 + h3 + + h4 + h5 - 64;

end

