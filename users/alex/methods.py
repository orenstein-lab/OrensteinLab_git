
autobalance(slope, tolerance, offset=0, bal_axis=2, var='Demod 1 x', ikwargs_dict={}, mobj_dict={}, iobj_dict={}, print_flag=True, close_devices=True)

def autobalance_at_initial_time(pos, mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict):

    time_index =0
    initial_time = -10
    slope=1
    tolerance=1
    balance_var = 'Demod 2 x'
    bal_axis = 2

    if pos[time_index]==initial_time:
        autobalance(slope, tolerance, offset=0, bal_axis=bal_axis, var=balance_var, ikwargs_dict=ikwargs_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict)

    return 0
