
function [] = dump_window_function(fo, source, name, f, M)
    w = f(M);
    for i = 1:M
        fprintf(fo, '%-6s %-24s %6d %6d %53.40f\n', source, name, M, i, w(i));
    end
end
