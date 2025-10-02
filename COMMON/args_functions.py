
from datetime import datetime as dt
from datetime import timedelta
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