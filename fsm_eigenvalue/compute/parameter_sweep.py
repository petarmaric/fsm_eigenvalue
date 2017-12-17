from contextlib import contextmanager
import itertools
import multiprocessing

from beam_integrals.beam_types import BaseBeamType
from simple_plugins import AttrDict

from .core import get_modal_composite, perform_iteration
from .integral_db import check_for_integral_db, open_integral_db


def _init_pool(*data):
    global _pool_data

    data_keys = 'beam_type_id, search_space, strip_data, materials, astiff_shape'.split(', ')
    _pool_data = AttrDict(zip(data_keys, data))

    _pool_data.beam_type = BaseBeamType.coerce(_pool_data.beam_type_id)
    _pool_data.integral_db = open_integral_db(_pool_data.beam_type_id)

def _worker(args):
    a, t_b = args
    c = _pool_data

    raw_results = [
        perform_iteration(c.integral_db, c.beam_type, c.strip_data, c.materials, c.astiff_shape, a, t_b, m)
        for m in c.search_space['m']
    ]
    modal_composite = get_modal_composite(raw_results)

    return a, t_b, raw_results, modal_composite

@contextmanager
def parameter_sweep(beam_type_id, search_space, strip_data, materials, astiff_shape, purge_integral_db_cache=False):
    check_for_integral_db(beam_type_id, purge_cache=purge_integral_db_cache)

    try:
        pool = multiprocessing.Pool(
            initializer=_init_pool,
            initargs=(beam_type_id, search_space, strip_data, materials, astiff_shape),
        )

        yield pool.imap(
            func=_worker,
            iterable=itertools.product(search_space['a'], search_space['t_b']),
            chunksize=len(search_space['t_b'])
        )
    finally:
        pool.terminate()
