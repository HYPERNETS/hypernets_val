
from OPTIONS import OptionsManager as om
from datetime import datetime as dt
from datetime import timedelta
import json

def get_start_end_date_from_args(args):
    start_date = None
    end_date = None
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except:
            try:
                srel = int(args.start_date)
                start_date = dt.now().replace(hour=0,minute=0,second=0,microsecond=0)+timedelta(days=srel)
            except:
                print(f'[WARNING] --start_date (-sd)  could not be parsed from {args.start_date}. Format should be YYYY-mm-dd or an integer value for relative selection')
                start_date = None


    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except:
            try:
                srel = int(args.end_date)
                end_date = dt.now().replace(hour=0,minute=0,second=0,microsecond=0)+timedelta(days=srel)
            except:
                print(f'[WARNING] --end_date (-sd)  could not be parsed from {args.end_date}. Format should be YYYY-mm-dd or an integer value for relative selection')
                end_date = None
    if start_date is not None and end_date is None:
        end_date = start_date
        print(f'[WARNING] As end date was not specified, it was set equal to start date.')

    if args.verbose and start_date is not None and end_date is not None:
        print(f'[INFO] Start date: {start_date.strftime("%Y-%m-%d")} End date: {end_date.strftime("%Y-%m-%d")}')

    return start_date,end_date



def check_args_required_from_list(args,list_required):
    args_dict = vars(args)
    for arg in args_dict:
        if arg not in list_required:
            return False
    return True

def check_arg_type(args,arg,type):
    args_dict = vars(args)
    if not arg in args_dict:
        return False
    arg_value = args_dict[arg]
    return check_arg_type_impl(arg_value,type)

def check_arg_type_impl(arg_value,type):
    res = om.get_value_param_impl(arg_value, type, None)
    if res is None:
        return False
    return True

def get_args_as_dict(args,required_dict,include_dates):
    args_dict = vars(args)
    res = {}
    for arg in required_dict:
        if not arg in args_dict:
            print(f'[ERROR] Argument {arg} is required. Please add it in your script call.')
            return None
        type = required_dict[arg]['type']
        pvalues = required_dict[arg]['type']['potential_values'] if 'potential_values' in required_dict[arg] else None
        arg_value = args_dict[arg]
        res[arg] = om.get_value_param_impl(arg_value,type,None)
        if res[arg] is None:
            print(f'[ERROR] Argument {arg} is not valid. Please check it.')
            return None
        if pvalues is not None and res[art] not in pvalues:
            print(f'[ERROR] Argument {arg} should take a value among {pvalues}')
            return None
    if include_dates:
        start_date, end_date = get_start_end_date_from_args(args)
        if start_date is None or end_date is None:
            return None
        res['start_date'] = start_date
        res['end_date'] = end_date

    return res
