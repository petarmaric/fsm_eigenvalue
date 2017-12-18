from . import DEFAULT_PAGINATE_BY
from .compute import parameter_sweep
from .load import load_data_from
from .store import store_results_to


def do_everything(data_file, results_file, purge_integral_db_cache=False, paginate_by=DEFAULT_PAGINATE_BY):
    beam_type_id, search_space, _, strip_data, materials, astiff_shape = load_data_from(data_file)

    with parameter_sweep(beam_type_id, search_space, strip_data, materials, astiff_shape, purge_integral_db_cache) as results_iterator:
        store_results_to(results_file, data_file, search_space, astiff_shape, results_iterator, paginate_by)
