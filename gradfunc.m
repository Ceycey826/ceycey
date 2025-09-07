function ypr = gradfunc(x)
    % Gradient of the Camel Back (Six-Hump) Function

    x1 = x(1);
    x2 = x(2);

    df_dx1 = 8*x1 - 8.4*x1^3 + 2*x1^5 + x2;
    df_dx2 = x1 - 8*x2 + 16*x2^3;

    ypr = [df_dx1; df_dx2];
end

