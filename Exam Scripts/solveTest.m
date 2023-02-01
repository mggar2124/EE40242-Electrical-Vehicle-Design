clear
clc


syms x y z a b 

y = 1;
%z = 1;
x = 3;
b = 0.5;
a = 3;

%% Equations

x = (2*a) + b;

equation = x + y + z == 10;

%% Solver
symbol = z;

solved = vpa(solve(equation, symbol))