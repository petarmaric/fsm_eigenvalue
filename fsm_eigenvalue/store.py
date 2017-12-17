from contextlib import contextmanager
from datetime import datetime
import logging
from timeit import default_timer as timer

import numpy as np
import tables as tb
from tzlocal import get_localzone
import yaml

from . import __version__, DEFAULT_PAGINATE_BY


logger = logging.getLogger(__name__)


def vector_f64(shape):
    return np.float64, shape

RAW_RESULTS_TABLE_SPEC = [
    # column_name,       dtype,      unit,    description
    ('a',                np.float64, 'mm',    'strip length'),
    ('t_b',              np.float64, 'mm',    'base strip thickness'),
    ('m',                np.int32,   '',      'mode'),
    ('omega',            np.float64, 'rad/s', 'natural frequency'),
    ('omega_approx',     np.float64, 'rad/s', 'natural frequency approximated from critical buckling stress'),
    ('omega_rel_err',    np.float64, '',      'natural frequency relative approximation error'),
    ('sigma_cr',         np.float64, 'MPa',   'critical buckling stress'),
    ('sigma_cr_approx',  np.float64, 'MPa',   'critical buckling stress approximated from natural frequency'),
    ('sigma_cr_rel_err', np.float64, '',      'critical buckling stress relative approximation error'),
    ('Phi_omega',        vector_f64, '',      'natural frequency mode shape'),
    ('Phi_sigma_cr',     vector_f64, '',      'critical buckling stress mode shape'),
    ('Phi_rel_err',      vector_f64, '',      'mode shape relative error'),
]

MODAL_COMPOSITES_TABLE_SPEC = [
    # column_name,       dtype,      unit,    description
    ('a',                np.float64, 'mm',    'strip length'),
    ('t_b',              np.float64, 'mm',    'base strip thickness'),
    ('m_dominant',       np.int32,   '',      'dominant mode, modal composite via sigma_cr'),
    ('omega',            np.float64, 'rad/s', 'natural frequency'),
    ('omega_approx',     np.float64, 'rad/s', 'natural frequency approximated from critical buckling stress'),
    ('omega_rel_err',    np.float64, '',      'natural frequency relative approximation error'),
    ('sigma_cr',         np.float64, 'MPa',   'critical buckling stress'),
    ('sigma_cr_approx',  np.float64, 'MPa',   'critical buckling stress approximated from natural frequency'),
    ('sigma_cr_rel_err', np.float64, '',      'critical buckling stress relative approximation error'),
]


def get_hdf5_table_description(table_spec, vector_shape):
    return np.dtype([
        (column_name, dtype(vector_shape) if callable(dtype) else dtype) # vectorize the ``dtype``, if callable
        for column_name, dtype, _, _ in table_spec
    ])

def get_column_units(table_spec):
    return {
        column_name: unit
        for column_name, _, unit, _ in table_spec
    }

def get_column_descriptions(table_spec):
    return {
        column_name: description
        for column_name, _, _, description in table_spec
    }

@contextmanager
def create_table(file, group, name, table_spec, vector_shape, expectedrows, indexes=None):
    # Use `expectedrows` to help PyTables determine the optimal chunk size
    table_description = get_hdf5_table_description(table_spec, vector_shape)
    table = file.create_table(group, name, table_description, expectedrows=expectedrows)

    # Add table metadata
    table.attrs.column_units_as_yaml = yaml.dump(
        get_column_units(table_spec),
        default_flow_style=False
    )
    table.attrs.column_descriptions_as_yaml = yaml.dump(
        get_column_descriptions(table_spec),
        default_flow_style=False
    )

    yield table

    # Flush the table manually, or the last data chunk won't be filled with correct values
    table.flush()

    if indexes:
        logger.info("Creating a completely sorted index (CSI) on %s columns to speed up '%s' table lookups... ", indexes, name)
        start = timer()
        for col in indexes:
            table.cols._f_col(col).create_csindex()
        logger.info("Index creation completed in %f second(s)", timer() - start)

    table.close()

def store_results_to(results_file, data_file, search_space, astiff_shape, results_iterator, paginate_by=DEFAULT_PAGINATE_BY):
    # These filters will be applied to all the datasets created immediately under the root group:
    #   * 'complib': Specifies the compression library to be used. Although PyTables
    #     supports many interesting compression libraries, HDF5 itself provides
    #     only 2 pre-defined filters for compression: ZLIB and SZIP. We can't
    #     use SZIP due to licensing issues, therefore ZLIB has been chosen by default as
    #     it's supported by all major HDF5 viewers (HDFView, HDF Compass, ViTables,
    #     HDF Explorer).
    #   * 'complevel': Specifies a compression level for data. Using the lowest
    #     level (1) by default, per PyTables optimization recommendations (see references).
    #   * 'shuffle': Enable the Shuffle filter to improve the compression ratio.
    #   * 'fletcher32': Enable the Fletcher32 filter to add a checksum on each
    #     data chunk.
    #
    # References:
    #   * https://www.hdfgroup.org/services/filters.html
    #   * https://www.hdfgroup.org/hdf5-quest.html#gcomp
    #   * https://www.hdfgroup.org/HDF5/faq/compression.html
    #   * http://www.pytables.org/usersguide/libref/helper_classes.html#the-filters-class
    #   * http://www.pytables.org/usersguide/optimization.html#compression-issues
    #   * http://www.pytables.org/usersguide/optimization.html#shuffling-or-how-to-make-the-compression-process-more-effective
    filters = tb.Filters(complib='zlib', complevel=1, shuffle=True, fletcher32=True)

    start = timer()
    with tb.open_file(results_file, 'w', filters=filters) as out:
        # Add the root group metadata
        out.root._v_attrs.generator_name = 'fsm_eigenvalue'
        out.root._v_attrs.generator_version = __version__
        out.root._v_attrs.created_at = datetime.now(get_localzone()).replace(
            microsecond=0
        ).isoformat()

        # Add data file contents
        with open(data_file, 'r') as fp:
            out.root._v_attrs.data_file = fp.read()

        astiff_size = astiff_shape[0]
        num_iterations = len(search_space['a']) * len(search_space['t_b'])

        logger.info('Performing a multi-dimensional parameter sweep and storing its results...')
        parameter_sweep_group = out.create_group(out.root, 'parameter_sweep')
        with create_table(
            out, parameter_sweep_group, 'raw_results',
            table_spec=RAW_RESULTS_TABLE_SPEC,
            vector_shape=astiff_size,
            expectedrows=num_iterations * len(search_space['m']),
            indexes=['a', 't_b', 'm']
        ) as raw_results_table, \
        create_table(
            out, parameter_sweep_group, 'modal_composites',
            table_spec=MODAL_COMPOSITES_TABLE_SPEC,
            vector_shape=astiff_size,
            expectedrows=num_iterations,
            indexes=['a', 't_b']
        ) as modal_composites_table:
            num_iterations_digits = np.ceil(np.log10(num_iterations))
            progress_fmt = "%6.2f%% (%{0}d/%{0}d iterations)".format(num_iterations_digits)
            for index, (_, _, raw_results, modal_composite) in enumerate(results_iterator, start=1):
                raw_results_table.append(raw_results) # Bulk insert
                modal_composites_table.append([modal_composite])

                if index % paginate_by == 0:
                    logger.info(progress_fmt, 100.0 * index / num_iterations, index, num_iterations)

            elapsed = timer() - start
            logger.info("Completed in %.2f second(s), %.3f millisecond(s) per iteration", elapsed, 1000.0 * elapsed/num_iterations)
