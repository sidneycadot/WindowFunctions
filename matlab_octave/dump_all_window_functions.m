
if exist('OCTAVE_VERSION', 'builtin')
    disp('running in OCTAVE');
    pkg('load', 'signal');
    filename = 'octave_windows.txt';
else
    disp('running in MATLAB');
    filename = 'matlab_windows.txt';
end

fo = fopen(filename, 'w');

minM = 1;
maxM = 100;

for M = minM:maxM

    disp(M);

    dump_window_function(fo, 'barthannwin'              , @(M) barthannwin    (M             ), M);

    dump_window_function(fo, 'bartlett'                 , @(M) bartlett       (M             ), M);

    dump_window_function(fo, 'blackman'                 , @(M) blackman       (M             ), M);
    dump_window_function(fo, 'blackman_periodic'        , @(M) blackman       (M, 'periodic' ), M);
    dump_window_function(fo, 'blackman_symmetric'       , @(M) blackman       (M, 'symmetric'), M);

    dump_window_function(fo, 'blackmanharris'           , @(M) blackmanharris (M             ), M);
    dump_window_function(fo, 'blackmanharris_periodic'  , @(M) blackmanharris (M, 'periodic' ), M);
    dump_window_function(fo, 'blackmanharris_symmetric' , @(M) blackmanharris (M, 'symmetric'), M);

    dump_window_function(fo, 'bohmanwin'                , @(M) bohmanwin      (M             ), M);

    % The chebwin() function has an optional 'R' parameter that defaults to 100 if not specified.

    dump_window_function(fo, 'chebwin'                  , @(M) chebwin        (M             ), M);
    dump_window_function(fo, 'chebwin_100p0'            , @(M) chebwin        (M,  100.0     ), M);
    dump_window_function(fo, 'chebwin_120p0'            , @(M) chebwin        (M,  120.0     ), M);

    dump_window_function(fo, 'flattopwin'               , @(M) flattopwin     (M             ), M);
    dump_window_function(fo, 'flattopwin_periodic'      , @(M) flattopwin     (M, 'periodic' ), M);
    dump_window_function(fo, 'flattopwin_symmetric'     , @(M) flattopwin     (M, 'symmetric'), M);

    % The gausswin() function has an optional 'Alpha' parameter that defaults to 2.5 if not specified.

    dump_window_function(fo, 'gausswin'                 , @(M) gausswin       (M             ), M);
    dump_window_function(fo, 'gausswin_2p5'             , @(M) gausswin       (M,  2.5       ), M);
    dump_window_function(fo, 'gausswin_3p2'             , @(M) gausswin       (M,  3.2       ), M);

    dump_window_function(fo, 'hamming'                  , @(M) hamming        (M             ), M);
    dump_window_function(fo, 'hamming_periodic'         , @(M) hamming        (M, 'periodic' ), M);
    dump_window_function(fo, 'hamming_symmetric'        , @(M) hamming        (M, 'symmetric'), M);

    dump_window_function(fo, 'hann'                     , @(M) hann           (M             ), M);
    dump_window_function(fo, 'hann_periodic'            , @(M) hann           (M, 'periodic' ), M);
    dump_window_function(fo, 'hann_symmetric'           , @(M) hann           (M, 'symmetric'), M);

    dump_window_function(fo, 'hanning'                  , @(M) hanning        (M             ), M);
    dump_window_function(fo, 'hanning_periodic'         , @(M) hanning        (M, 'periodic' ), M);
    dump_window_function(fo, 'hanning_symmetric'        , @(M) hanning        (M, 'symmetric'), M);

    % The kaiser() function has an optional 'beta' parameter that defaults to 0.5 if not specified.

    dump_window_function(fo, 'kaiser'                   , @(M) kaiser         (M             ), M);
    dump_window_function(fo, 'kaiser_0p5'               , @(M) kaiser         (M,  0.5       ), M);
    dump_window_function(fo, 'kaiser_0p8'               , @(M) kaiser         (M,  0.8       ), M);

    dump_window_function(fo, 'nuttallwin'               , @(M) nuttallwin     (M             ), M);
    dump_window_function(fo, 'nuttallwin_periodic'      , @(M) nuttallwin     (M, 'periodic' ), M);
    dump_window_function(fo, 'nuttallwin_symmetric'     , @(M) nuttallwin     (M, 'symmetric'), M);

    dump_window_function(fo, 'parzenwin'                , @(M) parzenwin      (M             ), M);

    dump_window_function(fo, 'rectwin'                  , @(M) rectwin        (M             ), M);

    if exist('taylorwin', 'file')

        % The taylorwin() function has optional 'nbar' and 'sll' arguments that default to 4 and -30
        % if not specified.
        %
        % The taylorwin() function is not available in Octave (v6.2.0).

        dump_window_function(fo, 'taylorwin'            , @(M) taylorwin      (M             ), M);
        dump_window_function(fo, 'taylorwin_4'          , @(M) taylorwin      (M,  4         ), M);
        dump_window_function(fo, 'taylorwin_5'          , @(M) taylorwin      (M,  5         ), M);
        dump_window_function(fo, 'taylorwin_6'          , @(M) taylorwin      (M,  6         ), M);
        dump_window_function(fo, 'taylorwin_4_m20'      , @(M) taylorwin      (M,  4, -20    ), M);
        dump_window_function(fo, 'taylorwin_4_m30'      , @(M) taylorwin      (M,  4, -30    ), M);
        dump_window_function(fo, 'taylorwin_4_m40'      , @(M) taylorwin      (M,  4, -40    ), M);
        dump_window_function(fo, 'taylorwin_5_m20'      , @(M) taylorwin      (M,  5, -20    ), M);
        dump_window_function(fo, 'taylorwin_5_m30'      , @(M) taylorwin      (M,  5, -30    ), M);
        dump_window_function(fo, 'taylorwin_5_m40'      , @(M) taylorwin      (M,  5, -40    ), M);
        dump_window_function(fo, 'taylorwin_6_m20'      , @(M) taylorwin      (M,  6, -20    ), M);
        dump_window_function(fo, 'taylorwin_6_m30'      , @(M) taylorwin      (M,  6, -30    ), M);
        dump_window_function(fo, 'taylorwin_6_m40'      , @(M) taylorwin      (M,  6, -40    ), M);

    end

    dump_window_function(fo, 'triang'                   , @(M) triang         (M             ), M);

    dump_window_function(fo, 'tukeywin'                 , @(M) tukeywin       (M             ), M);
    dump_window_function(fo, 'tukeywin_0p0'             , @(M) tukeywin       (M,  0.0       ), M);
    dump_window_function(fo, 'tukeywin_0p2'             , @(M) tukeywin       (M,  0.2       ), M);
    dump_window_function(fo, 'tukeywin_0p5'             , @(M) tukeywin       (M,  0.5       ), M);
    dump_window_function(fo, 'tukeywin_0p8'             , @(M) tukeywin       (M,  0.8       ), M);
    dump_window_function(fo, 'tukeywin_1p0'             , @(M) tukeywin       (M,  1.0       ), M);

end

fclose(fo);
