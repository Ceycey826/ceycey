function yprpr = hessianfunc(x)
    % Hessian of the Camel Back (Six-Hump) Function

    x1 = x(1);
    x2 = x(2);

    d2f_dx1x1 = 8 - 25.2*x1^2 + 10*x1^4;
    d2f_dx1x2 = 1;
    d2f_dx2x2 = -8 + 48*x2^2;

    yprpr = [d2f_dx1x1, d2f_dx1x2;
             d2f_dx1x2, d2f_dx2x2];
end
