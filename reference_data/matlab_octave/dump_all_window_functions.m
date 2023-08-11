
if exist('OCTAVE_VERSION', 'builtin')
    disp('Running in OCTAVE ...');
    pkg('load', 'signal'); % In Octave, the 'signal' package contains the window-function functiononality.
    program = 'Octave';
    source = 'octave';
    filename = 'octave_windows.txt';
else
    program = 'Matlab';
    source = 'matlab';
    filename = 'matlab_windows.txt';
end

fo = fopen(filename, 'w');

fprintf(fo, "# %s version %s\n", program, version);
fprintf(fo, "# %s\n", datestr(now));

minM = 1;
maxM = 100;

for M = minM:maxM

    disp(M);

    dump_window_function(fo, source, 'barthannwin'              , @(M) barthannwin    (M             ), M);

    dump_window_function(fo, source, 'bartlett'                 , @(M) bartlett       (M             ), M);

    dump_window_function(fo, source, 'blackman'                 , @(M) blackman       (M             ), M);
    dump_window_function(fo, source, 'blackman_periodic'        , @(M) blackman       (M, 'periodic' ), M);
    dump_window_function(fo, source, 'blackman_symmetric'       , @(M) blackman       (M, 'symmetric'), M);

    dump_window_function(fo, source, 'blackmanharris'           , @(M) blackmanharris (M             ), M);
    dump_window_function(fo, source, 'blackmanharris_periodic'  , @(M) blackmanharris (M, 'periodic' ), M);
    dump_window_function(fo, source, 'blackmanharris_symmetric' , @(M) blackmanharris (M, 'symmetric'), M);

    dump_window_function(fo, source, 'bohmanwin'                , @(M) bohmanwin      (M             ), M);

    dump_window_function(fo, source, 'boxcar'                   , @(M) boxcar         (M             ), M);

    % The chebwin() function has an optional 'R' parameter that defaults to 100 if not specified.

    dump_window_function(fo, source, 'chebwin'                  , @(M) chebwin        (M             ), M);
    dump_window_function(fo, source, 'chebwin_100p0'            , @(M) chebwin        (M,  100.0     ), M);
    dump_window_function(fo, source, 'chebwin_120p0'            , @(M) chebwin        (M,  120.0     ), M);

    dump_window_function(fo, source, 'flattopwin'               , @(M) flattopwin     (M             ), M);
    dump_window_function(fo, source, 'flattopwin_periodic'      , @(M) flattopwin     (M, 'periodic' ), M);
    dump_window_function(fo, source, 'flattopwin_symmetric'     , @(M) flattopwin     (M, 'symmetric'), M);

    % The gausswin() function has an optional 'Alpha' parameter that defaults to 2.5 if not specified.

    dump_window_function(fo, source, 'gausswin'                 , @(M) gausswin       (M             ), M);
    dump_window_function(fo, source, 'gausswin_2p5'             , @(M) gausswin       (M,  2.5       ), M);
    dump_window_function(fo, source, 'gausswin_3p2'             , @(M) gausswin       (M,  3.2       ), M);

    dump_window_function(fo, source, 'hamming'                  , @(M) hamming        (M             ), M);
    dump_window_function(fo, source, 'hamming_periodic'         , @(M) hamming        (M, 'periodic' ), M);
    dump_window_function(fo, source, 'hamming_symmetric'        , @(M) hamming        (M, 'symmetric'), M);

    dump_window_function(fo, source, 'hann'                     , @(M) hann           (M             ), M);
    dump_window_function(fo, source, 'hann_periodic'            , @(M) hann           (M, 'periodic' ), M);
    dump_window_function(fo, source, 'hann_symmetric'           , @(M) hann           (M, 'symmetric'), M);

    dump_window_function(fo, source, 'hanning'                  , @(M) hanning        (M             ), M);
    dump_window_function(fo, source, 'hanning_periodic'         , @(M) hanning        (M, 'periodic' ), M);
    dump_window_function(fo, source, 'hanning_symmetric'        , @(M) hanning        (M, 'symmetric'), M);

    % The kaiser() function has an optional 'beta' parameter that defaults to 0.5 if not specified.

    dump_window_function(fo, source, 'kaiser'                   , @(M) kaiser         (M             ), M);
    dump_window_function(fo, source, 'kaiser_0p5'               , @(M) kaiser         (M,  0.5       ), M);
    dump_window_function(fo, source, 'kaiser_0p8'               , @(M) kaiser         (M,  0.8       ), M);

    dump_window_function(fo, source, 'nuttallwin'               , @(M) nuttallwin     (M             ), M);
    dump_window_function(fo, source, 'nuttallwin_periodic'      , @(M) nuttallwin     (M, 'periodic' ), M);
    dump_window_function(fo, source, 'nuttallwin_symmetric'     , @(M) nuttallwin     (M, 'symmetric'), M);

    dump_window_function(fo, source, 'parzenwin'                , @(M) parzenwin      (M             ), M);

    dump_window_function(fo, source, 'rectwin'                  , @(M) rectwin        (M             ), M);

    if exist('taylorwin', 'file')

        % The taylorwin() function has optional 'nbar' and 'sll' arguments that default to 4 and -30
        % if not specified.
        %
        % Note: the taylorwin() function is not available in Octave (v6.2.0).

        dump_window_function(fo, source, 'taylorwin'            , @(M) taylorwin      (M             ), M);
        dump_window_function(fo, source, 'taylorwin_3'          , @(M) taylorwin      (M,  3         ), M);
        dump_window_function(fo, source, 'taylorwin_4'          , @(M) taylorwin      (M,  4         ), M);
        dump_window_function(fo, source, 'taylorwin_5'          , @(M) taylorwin      (M,  5         ), M);
        dump_window_function(fo, source, 'taylorwin_3_m20'      , @(M) taylorwin      (M,  3, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_3_m30'      , @(M) taylorwin      (M,  3, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_3_m40'      , @(M) taylorwin      (M,  3, -40    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m20'      , @(M) taylorwin      (M,  4, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m30'      , @(M) taylorwin      (M,  4, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m40'      , @(M) taylorwin      (M,  4, -40    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m20'      , @(M) taylorwin      (M,  5, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m30'      , @(M) taylorwin      (M,  5, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m40'      , @(M) taylorwin      (M,  5, -40    ), M);

    end

    dump_window_function(fo, source, 'triang'                   , @(M) triang         (M             ), M);

    dump_window_function(fo, source, 'tukeywin'                 , @(M) tukeywin       (M             ), M);
    dump_window_function(fo, source, 'tukeywin_0p0'             , @(M) tukeywin       (M,  0.0       ), M);
    dump_window_function(fo, source, 'tukeywin_0p2'             , @(M) tukeywin       (M,  0.2       ), M);
    dump_window_function(fo, source, 'tukeywin_0p5'             , @(M) tukeywin       (M,  0.5       ), M);
    dump_window_function(fo, source, 'tukeywin_0p8'             , @(M) tukeywin       (M,  0.8       ), M);
    dump_window_function(fo, source, 'tukeywin_1p0'             , @(M) tukeywin       (M,  1.0       ), M);

end

fclose(fo);
