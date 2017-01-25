
function [] = dump_window_function(fo, name, f, M)
    w = f(M);
    for i = 1:M
        fprintf(fo, '%-29s %9d %9d %49.30f\n', name, M, i, w(i));
    end
end
