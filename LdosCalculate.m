% To continue with dispersion calculation:

k = input('Enter k point interested : ');
band = input('Enter number(order) of bands interested(1-6): ');

kk = 3.14*k;

HH = Alpha + Beta*exp(1i*kk)+ Beta'*exp(-1i*kk);
[VV,DD] = eig(HH);
eigEE = diag(DD);

Band = NW*NU/2-3+band;

DOS_k = zeros(NU*NW,1);
for i = 1:NU*NW
    DOS_k(i) = VV(i,Band)'*VV(i,Band); 
end

