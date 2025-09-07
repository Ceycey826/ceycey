function y = func(x)
    % Camel Back - Six Hump Function
    % f(x) = 4*x1^2 - 2.1*x1^4 + (1/3)*x1^6 + x1*x2 - 4*x2^2 + 4*x2^4

    x1 = x(1);
    x2 = x(2);

    y = 4*x1^2 - 2.1*x1^4 + (1/3)*x1^6 + x1*x2 - 4*x2^2 + 4*x2^4;
end

