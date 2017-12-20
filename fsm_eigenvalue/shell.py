import argparse
import logging
import os

from . import __version__, DEFAULT_PAGINATE_BY
from .main import do_everything


def main():
    # Setup command line option parser
    parser = argparse.ArgumentParser(
        description='Parametric modeling of buckling and free vibration in '\
                    'prismatic shell structures, performed by solving the '\
                    'eigenvalue problem in HCFSM.'
    )
    parser.add_argument(
        'data_file',
        help="Data file describing the parametric model, please see "\
             "'examples/data-files/barbero-viscoelastic.yaml' for an example"
    )
    parser.add_argument(
        '-r',
        '--results-file',
        metavar='FILENAME',
        help="Store results to the selected FILENAME, uses '<data_file>.hdf5' by default"
    )
    parser.add_argument(
        '-d',
        '--purge-integral-db-cache',
        action='store_true',
        help='Purge the integral db cache, forcing it to redownload'
    )
    parser.add_argument(
        '-p',
        '--paginate-by',
        metavar='NUM',
        type=int,
        default=DEFAULT_PAGINATE_BY,
        help="Show progress every NUM iterations, %d by default" % DEFAULT_PAGINATE_BY
    )
    parser.add_argument(
        '-q',
        '--quiet',
        action='store_const',
        const=logging.WARN,
        dest='verbosity',
        help='Be quiet, show only warnings and errors'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_const',
        const=logging.DEBUG,
        dest='verbosity',
        help='Be very verbose, show debug information'
    )
    parser.add_argument(
        '--version',
        action='version',
        version="%(prog)s " + __version__
    )
    args = parser.parse_args()

    # Configure logging
    log_level = args.verbosity or logging.INFO
    logging.basicConfig(level=log_level, format="%(asctime)s [%(levelname)s] %(message)s")

    if not args.results_file:
        args.results_file = os.path.splitext(args.data_file)[0] + '.hdf5'

    do_everything(
        data_file=args.data_file,
        results_file=args.results_file,
        purge_integral_db_cache=args.purge_integral_db_cache,
        paginate_by=args.paginate_by,
    )

if __name__ == '__main__':
    main()
