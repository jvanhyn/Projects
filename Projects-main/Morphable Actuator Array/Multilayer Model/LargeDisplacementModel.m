function [XL,YL,X,Y] = LargeDisplacementModel(Ex,tx,Es,ts,Ep,tp,Ec,tc,Em,tm,V,L,n,w,d31)
%LARGEDISPLACEMENTMODEL Summary of this function goes here
%   Detailed explanation goes here
lseg = L/n;

G=w*(Ex*tx+2*Es*ts+4*Ec*tc+Ep*tp+Em*tm);
H=0.5*w*(Ep*tp^2+2*Ep*(tx/2+ts+tc)*tp-Em*tm^2-2*Em*(tx/2+ts+tc)*tm);
F=(w/3)*(2*Ex*(tx/2)^3+2*Es*(3*(tx/2)*ts*(tx/2+ts)+ts^3)+2*Ec*(6*(tx/2+ts)*tc*(tx/2+ts+tc)+3*(tc+tp)*tc*(tx/2+ts+tc)+3*(tx/2+ts)*tc*(tp+tc)+3*((tc+tp)^2)*tc+2*tc^3)+Ep*(3*(tx/2+ts+tp)*tp*(tx/2+ts+tp+tc)+tp^3)+Em*(3*(tx/2+ts+tc)*tm*(tx/2+ts+tc+tm)+tm^3));
E3=-V/tp;
N1E=w*d31*Ep*E3*tp;
M2E=0.5*w*(d31*E3*Ep*(2*tp*(tx/2+ts+tc)+tp^2));

for i=1:1:n
    X(i)=((F*G-H^2)/(H*N1E-G*M2E))*sin(((H*N1E-G*M2E)/(F*G-H^2))*i*lseg)-i*lseg;                    % point i X coordinate
    XL(i)=(X(i)+i*(L/n))/L;                                                                         % point i normalized X coordinate
    Y(i)=((H^2-F*G)/(H*N1E-G*M2E))*cos(((H*N1E-G*M2E)/(F*G-H^2))*i*lseg)+((F*G-H^2)/(H*N1E-G*M2E)); % point i Y coordinate
    YL(i)=Y(i)/L;                                                                                   % point i normalized Y coordinate
end
end

