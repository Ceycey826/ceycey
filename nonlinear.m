clear all
close all
clc

X1 = -2:0.01:2;        
X2 = -2:0.01:2;
[x1, x2] = meshgrid(X1, X2);

% Camel Back (Six Hump) fonksiyonu
F = 4 * x1.^2 - 2.1 * x1.^4 + (1/3) * x1.^6 + x1 .* x2 - 4 * x2.^2 + 4 * x2.^4;

realFMin = min(min(F));
fprintf('Fonksiyonun grid üzerindeki minimum değeri: %.4f\n', realFMin)

figure
mesh(x1, x2, F)
xlabel('X1');
ylabel('X2');
zlabel('f(X1,X2)');
title("Camel (Six-Hump) Function - Mesh Plot")


start_points = -2 + 4 * rand(2,3);
colors = ['r', 'g', 'w'];


%% Newton-Raphson Algorithm on CB6 with 3 Random Starts
epsilon = 1e-4;

figure
contourf(x1, x2, F, 50)
xlabel('X1');
ylabel('X2');
title('Newton-Raphson on CB6 - 3 Random Starts in [-2,2]')
colorbar
hold on

for i = 1:3
    fprintf('\n--- Starting Point %d ---\n', i);

    x0 = start_points(:, i); 
    x = x0;
    color = colors(i);

    tic;

    fprintf('k=1, x1=%.4f, x2=%.4f, f(x)=%.4f\n', x(1), x(2), func(x));
    plot(x(1), x(2), [color '.']); hold on;

    % İlk adımdan önce Hessian kontrolü
    if rcond(hessianfunc(x)) < 1e-10
        warning('Hessian matrix is near-singular at initial point. Skipping this run.');
        continue
    end

    x_next = x - hessianfunc(x) \ gradfunc(x);
    fprintf('k=2, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
        x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), [color '*'])

    k = 3;
    while norm(gradfunc(x_next)) > epsilon && abs(func(x_next) - func(x)) > epsilon
        x = x_next;

        % Her iterasyonda güvenlik kontrolü
        H = hessianfunc(x);
        if rcond(H) < 1e-10
            warning('Hessian matrix became near-singular. Stopping this run.');
            break
        end

        x_next = x - H \ gradfunc(x);

        fprintf('k=%d, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
            k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
        plot(x_next(1), x_next(2), [color '*'])
        k = k + 1;
    end

    elapsed = toc;
    fprintf('Final Result: x1=%.4f, x2=%.4f, f(x)=%.4f (after %d steps, %.4f seconds)\n', ...
        x_next(1), x_next(2), func(x_next), k-1, elapsed);
end

set(gca, 'fontsize', 20)


%% Hestenes-Stiefel on CB6 - 3 Random Starting Points in [-2, 2]
epsilon = 1e-4;

figure
contourf(x1, x2, F, 50)
hold on
title('Hestenes-Stiefel on CB6 - 3 Random Starts (Safe)')
colorbar

for run = 1:3
    fprintf('\n--- Starting Point %d ---\n', run);
    
    x = start_points(:, run); 
    color = colors(run);
    fprintf('k=1, x1=%.4f, x2=%.4f, f(x)=%.4f\n', x(1), x(2), func(x));
    plot(x(1), x(2), [color '.'])

    tic 

    g = gradfunc(x);
    d = -g;

    alpha_range = 0:0.01:1;
    funcalpha = zeros(length(alpha_range),1);
    for i = 1:length(alpha_range)
        funcalpha(i) = func(x + alpha_range(i)*d);
    end
    [~, ind] = min(funcalpha);
    alpha = alpha_range(ind);

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);

    denom = d' * (g_next - g);
    if abs(denom) < 1e-8
        beta = 0;
    else
        beta = (g_next' * (g_next - g)) / denom;
    end
    d_next = -g_next + beta * d;

    fprintf('k=2, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
        x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), [color '*'])

    k = 3;
    while norm(gradfunc(x_next)) > epsilon || abs(func(x_next) - func(x)) > epsilon
        if any(isnan(x_next)) || any(isnan(g_next)) || isnan(func(x_next))
            fprintf('NaN detected — terminating this run.\n');
            break
        end

        x = x_next;
        g = g_next;
        d = d_next;

        funcalpha = zeros(length(alpha_range),1);
        for i = 1:length(alpha_range)
            funcalpha(i) = func(x + alpha_range(i)*d);
        end
        [~, ind] = min(funcalpha);
        alpha = alpha_range(ind);

        x_next = x + alpha * d;
        g_next = gradfunc(x_next);

        denom = d' * (g_next - g);
        if abs(denom) < 1e-8
            beta = 0;
        else
            beta = (g_next' * (g_next - g)) / denom;
        end
        d_next = -g_next + beta * d;

        fprintf('k=%d, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
            k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
        plot(x_next(1), x_next(2), [color '*'])
        k = k + 1;
    end

    elapsed = toc; 
    if ~any(isnan(x_next))
        fprintf('Final Result: x1=%.4f, x2=%.4f, f(x)=%.4f (after %d steps)\n', ...
            x_next(1), x_next(2), func(x_next), k-1);
    end
    fprintf('Elapsed time: %.6f seconds\n', elapsed); 
end

set(gca, 'fontsize', 20)

%% Polak–Ribiere on CB6 - 3 Random Starts
epsilon = 1e-4;
max_iter = 100;

figure
contourf(x1, x2, F, 50)
hold on
title('Polak-Ribiere on CB6 - 3 Random Starts')
colorbar

for run = 1:3
    fprintf('\n--- Starting Point %d ---\n', run);

    x = start_points(:, run);  
    color = colors(run);
    fprintf('k=1, x1=%.4f, x2=%.4f, f(x)=%.4f\n', x(1), x(2), func(x));
    plot(x(1), x(2), [color '.'])

    tic 

    g = gradfunc(x);
    d = -g;

    alpha_range = 0:0.01:1;
    funcalpha = zeros(length(alpha_range),1);
    for i = 1:length(alpha_range)
        funcalpha(i) = func(x + alpha_range(i)*d);
    end
    [~, ind] = min(funcalpha);
    alpha = alpha_range(ind);

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);

    denom = g' * g;
    if abs(denom) < 1e-8
        beta = 0;
    else
        beta = (g_next' * (g_next - g)) / denom;
    end
    d_next = -g_next + beta * d;

    fprintf('k=2, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
        x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), [color '*'])

    k = 3;
    while (norm(gradfunc(x_next)) > epsilon || abs(func(x_next) - func(x)) > epsilon) && k <= max_iter
        if any(isnan(x_next)) || any(isnan(g_next)) || isnan(func(x_next))
            fprintf('NaN detected — terminating.\n');
            break
        end
        if norm(x_next - x) < 1e-6
            fprintf('Step too small — assumed convergence.\n');
            break
        end
        if alpha < 1e-5 || norm(d_next) < 1e-8
            fprintf('Direction/alpha too small — stopping.\n');
            break
        end

        x = x_next;
        g = g_next;
        d = d_next;

        for i = 1:length(alpha_range)
            funcalpha(i) = func(x + alpha_range(i)*d);
        end
        [~, ind] = min(funcalpha);
        alpha = alpha_range(ind);

        x_next = x + alpha * d;
        g_next = gradfunc(x_next);

        denom = g' * g;
        if abs(denom) < 1e-8
            beta = 0;
        else
            beta = (g_next' * (g_next - g)) / denom;
        end
        d_next = -g_next + beta * d;

        fprintf('k=%d, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
            k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
        plot(x_next(1), x_next(2), [color '*'])
        k = k + 1;
    end

    elapsed = toc; 
    if ~any(isnan(x_next))
        fprintf('Final Result: x1=%.4f, x2=%.4f, f(x)=%.4f (after %d steps)\n', ...
            x_next(1), x_next(2), func(x_next), k-1);
    end
    fprintf('Elapsed time: %.6f seconds\n', elapsed); 
end

set(gca, 'fontsize', 20)

%% Fletcher–Reeves on CB6 - 3 Random Starts
epsilon = 1e-4;
max_iter = 100;

figure
contourf(x1, x2, F, 50)
hold on
title('Fletcher–Reeves on CB6 - 3 Random Starts')
colorbar

for run = 1:3
    fprintf('\n--- Starting Point %d ---\n', run);

    x = start_points(:, run); 
    color = colors(run);
    fprintf('k=1, x1=%.4f, x2=%.4f, f(x)=%.4f\n', x(1), x(2), func(x));
    plot(x(1), x(2), [color '.'])

    tic 

    g = gradfunc(x);
    d = -g;

    alpha_range = 0:0.01:1;
    funcalpha = zeros(length(alpha_range),1);
    for i = 1:length(alpha_range)
        funcalpha(i) = func(x + alpha_range(i)*d);
    end
    [~, ind] = min(funcalpha);
    alpha = alpha_range(ind);

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);

    denom = g' * g;
    if abs(denom) < 1e-8
        beta = 0;
    else
        beta = (g_next' * g_next) / denom;
    end
    d_next = -g_next + beta * d;

    fprintf('k=2, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
        x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
    plot(x_next(1), x_next(2), [color '*'])
    k = 3;

    while (norm(g_next) > epsilon || abs(func(x_next) - func(x)) > epsilon) && k <= max_iter
        if any(isnan(x_next)) || any(isnan(g_next)) || isnan(func(x_next))
            fprintf('NaN detected — terminating.\n'); break;
        end
        if norm(x_next - x) < 1e-6 || alpha < 1e-5 || norm(d_next) < 1e-8
            fprintf('Too small step/direction — terminating.\n'); break;
        end

        x = x_next;
        g = g_next;
        d = d_next;

        for i = 1:length(alpha_range)
            funcalpha(i) = func(x + alpha_range(i)*d);
        end
        [~, ind] = min(funcalpha);
        alpha = alpha_range(ind);

        x_next = x + alpha * d;
        g_next = gradfunc(x_next);

        denom = g' * g;
        if abs(denom) < 1e-8
            beta = 0;
        else
            beta = (g_next' * g_next) / denom;
        end
        d_next = -g_next + beta * d;

        fprintf('k=%d, x1=%.4f, x2=%.4f, f(x)=%.4f, abs. error=%.4f\n', ...
            k, x_next(1), x_next(2), func(x_next), abs(func(x_next) - func(x)));
        plot(x_next(1), x_next(2), [color '*'])
        k = k + 1;
    end

    elapsed = toc; 
    if ~any(isnan(x_next))
        fprintf('Final Result: x1=%.4f, x2=%.4f, f(x)=%.4f (after %d steps)\n', ...
            x_next(1), x_next(2), func(x_next), k-1);
    end
    fprintf('Elapsed time: %.6f seconds\n', elapsed);
end

set(gca, 'fontsize', 20)



