

def autobalance_at_initial_time(pos, mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict):

    time_index =0
    initial_time = -10
    slope=1
    tolerance=1
    balance_var = 'Demod 2 x'
    bal_axis = 2

    if pos[time_index]==initial_time:
        meas.autobalance(slope, tolerance, offset=0, bal_axis=bal_axis, var=balance_var, ikwargs_dict=ikwargs_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict)

    return mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict

def autobalance(pos, mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict):
    balance_var = 'Demod 2 x'
    bal_axis = 2
    slope_bal = -0.006080756498822003
    tol_bal = 500e-6
    mobj_dict_temp = {'axis_2':mobj_dict['corotate_axes12'][1]}
    meas.autobalance(slope_bal, tol_bal, bal_axis=bal_axis, var=balance_var, ikwargs_dict=ikwargs_dict, mobj_dict=mobj_dict_temp, iobj_dict=iobj_dict)
    return mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict
