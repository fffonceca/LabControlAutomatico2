function [prbs, estadoN] = prbs(N, estado0)

prbs = zeros([1,N]);
estadoN = estado0;

for x = 1:N
    resultado = xor(estadoN(1),estadoN(15));
    prbs(x) = resultado;
    estadoN = estadoN([15 1:14]);
    estadoN(1) = resultado;
end
end