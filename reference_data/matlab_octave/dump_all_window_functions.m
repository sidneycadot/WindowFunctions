
if exist('OCTAVE_VERSION', 'builtin')
    disp('Running in OCTAVE ...');
    pkg('load', 'signal'); % In Octave, most window functions require that the 'signal' package is loaded.
    program = 'Octave';
    source = 'octave';
    filename = 'octave_windows.txt';
else
    disp('Running in MATLAB ...');
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

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'barthannwin'               , @(M) barthannwin     (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'bartlett'                  , @(M) bartlett        (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'blackman'                  , @(M) blackman        (M             ), M);
        dump_window_function(fo, source, 'blackman_periodic'         , @(M) blackman        (M, 'periodic' ), M);
        dump_window_function(fo, source, 'blackman_symmetric'        , @(M) blackman        (M, 'symmetric'), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'blackmanharris'            , @(M) blackmanharris  (M             ), M);
        dump_window_function(fo, source, 'blackmanharris_periodic'   , @(M) blackmanharris  (M, 'periodic' ), M);
        dump_window_function(fo, source, 'blackmanharris_symmetric'  , @(M) blackmanharris  (M, 'symmetric'), M);
    end

    if ismember(source, {'not-in-matlab', 'octave'})
        dump_window_function(fo, source, 'blackmannuttall'           , @(M) blackmannuttall (M             ), M);
        dump_window_function(fo, source, 'blackmannuttall_periodic'  , @(M) blackmannuttall (M, 'periodic' ), M);
        dump_window_function(fo, source, 'blackmannuttall_symmetric' , @(M) blackmannuttall (M, 'symmetric'), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'bohmanwin'                 , @(M) bohmanwin       (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'boxcar'                    , @(M) boxcar          (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        % The chebwin() function has an optional 'R' parameter that defaults to 100 if not specified.
        dump_window_function(fo, source, 'chebwin'                   , @(M) chebwin         (M             ), M);
        dump_window_function(fo, source, 'chebwin_100p0'             , @(M) chebwin         (M,  100.0     ), M);
        dump_window_function(fo, source, 'chebwin_120p0'             , @(M) chebwin         (M,  120.0     ), M);
    end

    if ismember(source, {'not-in-matlab', 'octave'})
        % The expwin() function has an optional 'ALPHA' -or- 'SLL' parameter. If it is nonnegative, it is interpreted as
        % an 'alpha' value; if it is negative, it is interpreted as an SLL value in decibels.
        % It also has an optional "canonical" parameter.
        % It is only available in Octave, not in Matlab.
        dump_window_function(fo, source, 'expwin_2p5'                , @(M) expwin          (M,   2.5              ), M);
        dump_window_function(fo, source, 'expwin_3p2'                , @(M) expwin          (M,   3.2              ), M);
        dump_window_function(fo, source, 'expwin_minus_10'           , @(M) expwin          (M, -10.0              ), M);
        dump_window_function(fo, source, 'expwin_minus_20'           , @(M) expwin          (M, -20.0              ), M);
        dump_window_function(fo, source, 'expwin_2p5_canonical'      , @(M) expwin          (M,   2.5, 'canonical' ), M);
        dump_window_function(fo, source, 'expwin_3p2_canonical'      , @(M) expwin          (M,   3.2, 'canonical' ), M);
        dump_window_function(fo, source, 'expwin_minus_10_canonical' , @(M) expwin          (M, -10.0, 'canonical' ), M);
        dump_window_function(fo, source, 'expwin_minus_20_canonical' , @(M) expwin          (M, -20.0, 'canonical' ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'flattopwin'                , @(M) flattopwin      (M             ), M);
        dump_window_function(fo, source, 'flattopwin_periodic'       , @(M) flattopwin      (M, 'periodic' ), M);
        dump_window_function(fo, source, 'flattopwin_symmetric'      , @(M) flattopwin      (M, 'symmetric'), M);
    end

    if ismember(source, {'not-in-matlab', 'octave'})
        % The gaussian() function has an optional 'A' parameter that defaults to 1.0 if not specified.
        % It is only available in Octave, not in Matlab.
        dump_window_function(fo, source, 'gaussian'                  , @(M) gaussian        (M             ), M);
        dump_window_function(fo, source, 'gaussian_1p0'              , @(M) gaussian        (M,  1.0       ), M);
        dump_window_function(fo, source, 'gaussian_3p2'              , @(M) gaussian        (M,  2.5       ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        % The gausswin() function has an optional 'alpha' parameter that defaults to 2.5 if not specified.
        dump_window_function(fo, source, 'gausswin'                  , @(M) gausswin        (M             ), M);
        dump_window_function(fo, source, 'gausswin_2p5'              , @(M) gausswin        (M,  2.5       ), M);
        dump_window_function(fo, source, 'gausswin_3p2'              , @(M) gausswin        (M,  3.2       ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'hamming'                   , @(M) hamming         (M             ), M);
        dump_window_function(fo, source, 'hamming_periodic'          , @(M) hamming         (M, 'periodic' ), M);
        dump_window_function(fo, source, 'hamming_symmetric'         , @(M) hamming         (M, 'symmetric'), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'hann'                      , @(M) hann            (M             ), M);
        dump_window_function(fo, source, 'hann_periodic'             , @(M) hann            (M, 'periodic' ), M);
        dump_window_function(fo, source, 'hann_symmetric'            , @(M) hann            (M, 'symmetric'), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'hanning'                   , @(M) hanning         (M             ), M);
        dump_window_function(fo, source, 'hanning_periodic'          , @(M) hanning         (M, 'periodic' ), M);
        dump_window_function(fo, source, 'hanning_symmetric'         , @(M) hanning         (M, 'symmetric'), M);
    end


    if ismember(source, {'matlab', 'octave'})
        % The kaiser() function has an optional 'beta' parameter that defaults to 0.5 if not specified.
        dump_window_function(fo, source, 'kaiser'                    , @(M) kaiser          (M             ), M);
        dump_window_function(fo, source, 'kaiser_0p5'                , @(M) kaiser          (M,  0.5       ), M);
        dump_window_function(fo, source, 'kaiser_0p8'                , @(M) kaiser          (M,  0.8       ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'nuttallwin'                , @(M) nuttallwin      (M             ), M);
        dump_window_function(fo, source, 'nuttallwin_periodic'       , @(M) nuttallwin      (M, 'periodic' ), M);
        dump_window_function(fo, source, 'nuttallwin_symmetric'      , @(M) nuttallwin      (M, 'symmetric'), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'parzenwin'                 , @(M) parzenwin       (M             ), M);
    end

    if ismember(source, {'not-in-matlab', 'octave'})
        dump_window_function(fo, source, 'poisswin_2p5'              , @(M) poisswin        (M,  2.5       ), M);
        dump_window_function(fo, source, 'poisswin_3p2'              , @(M) poisswin        (M,  3.2       ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'rectwin'                   , @(M) rectwin         (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        % The taylorwin() function has optional 'nbar' and 'sll' arguments that default to 4 and -30
        % if not specified.
        dump_window_function(fo, source, 'taylorwin'                 , @(M) taylorwin       (M             ), M);
        dump_window_function(fo, source, 'taylorwin_3'               , @(M) taylorwin       (M,  3         ), M);
        dump_window_function(fo, source, 'taylorwin_4'               , @(M) taylorwin       (M,  4         ), M);
        dump_window_function(fo, source, 'taylorwin_5'               , @(M) taylorwin       (M,  5         ), M);
        dump_window_function(fo, source, 'taylorwin_3_m20'           , @(M) taylorwin       (M,  3, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_3_m30'           , @(M) taylorwin       (M,  3, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_3_m40'           , @(M) taylorwin       (M,  3, -40    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m20'           , @(M) taylorwin       (M,  4, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m30'           , @(M) taylorwin       (M,  4, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_4_m40'           , @(M) taylorwin       (M,  4, -40    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m20'           , @(M) taylorwin       (M,  5, -20    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m30'           , @(M) taylorwin       (M,  5, -30    ), M);
        dump_window_function(fo, source, 'taylorwin_5_m40'           , @(M) taylorwin       (M,  5, -40    ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'triang'                    , @(M) triang          (M             ), M);
    end

    if ismember(source, {'matlab', 'octave'})
        dump_window_function(fo, source, 'tukeywin'                  , @(M) tukeywin        (M             ), M);
        dump_window_function(fo, source, 'tukeywin_0p0'              , @(M) tukeywin        (M,  0.0       ), M);
        dump_window_function(fo, source, 'tukeywin_0p2'              , @(M) tukeywin        (M,  0.2       ), M);
        dump_window_function(fo, source, 'tukeywin_0p5'              , @(M) tukeywin        (M,  0.5       ), M);
        dump_window_function(fo, source, 'tukeywin_0p8'              , @(M) tukeywin        (M,  0.8       ), M);
        dump_window_function(fo, source, 'tukeywin_1p0'              , @(M) tukeywin        (M,  1.0       ), M);
    end

end

fclose(fo);
