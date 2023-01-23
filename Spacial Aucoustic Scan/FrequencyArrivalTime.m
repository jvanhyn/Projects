function [atime,ai, sign] = FrequencyArrivalTime(t,sig,thresh,cutoff,offset)

domain = cutoff:(length(sig)-offset);
fs = 1/mean(diff(t));
sigbp = bandpass(sig,[4500,6000], fs); 
sign = abs(sigbp/max(sigbp(domain)));

atime = 0.0;
ai = 0;

for i = domain
    if (sign(i)) - thresh >= 0
        atime = t(i);
        ai = i;
        break;
    end
end

end