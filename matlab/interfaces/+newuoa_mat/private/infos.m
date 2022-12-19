% This is a function defining exit flags.
%
% Coded by Zaikun Zhang in July 2020.
%
% Last Modified: Saturday, July 10, 2021 PM10:13:33

function info = infos(name)

    switch name
        case 'invalid_input'
            info = -100;
        case 'small_tr_radius'
            info = 0;
        case 'ftarget_achieved'
            info = 1;
        case 'trsubp_failed'
            info = 2;
        case 'maxfun_reached'
            info = 3;
        case 'maxtr_reached'
            info = 20;
        case 'nan_x'
            info = -1;
        case 'nan_inf_f'
            info = -2;
        case 'nan_model'
            info = -3;
        case 'damaging_rounding'
            info = 7;
        otherwise
            error('Error (info.m): unknown info %s', name);
    end
end
