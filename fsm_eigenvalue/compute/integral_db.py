import logging
import os
import shutil
from timeit import default_timer as timer

from beam_integrals.beam_types import BaseBeamType
from beam_integrals.integrals import BaseIntegral
import requests
import tables as tb

from .. import BASE_CACHE_DIR


logger = logging.getLogger(__name__)


INTEGRAL_DB_CACHE_DIR = os.path.join(BASE_CACHE_DIR, 'integrals')
INTEGRAL_DB_URL_FMT = "https://bitbucket.org/petar/export_beam_integrals/downloads/%s.hdf5"
INTEGRAL_DB_DOWNLOAD_BUFFER_SIZE = 100 * 1024


def get_integral_db_url(beam_type_id):
    beam_type = BaseBeamType.coerce(beam_type_id)
    return INTEGRAL_DB_URL_FMT % beam_type.filename

def get_integral_db_filename(beam_type_id):
    beam_type = BaseBeamType.coerce(beam_type_id)
    return os.path.join(INTEGRAL_DB_CACHE_DIR, "%s.hdf5" % beam_type.filename)

def download_integral_db(beam_type_id):
    db_url = get_integral_db_url(beam_type_id)
    db_filename = get_integral_db_filename(beam_type_id)
    db_dirname = os.path.dirname(db_filename)

    if not os.path.exists(db_dirname):
        logger.info("Creating the integral db cache directory '%s'...", INTEGRAL_DB_CACHE_DIR)
        os.makedirs(db_dirname)

    logger.info("Downloading %s...", db_url)
    start = timer()
    with requests.get(db_url, stream=True) as r:
        with open(db_filename, 'wb') as fp:
            fp.writelines(r.iter_content(INTEGRAL_DB_DOWNLOAD_BUFFER_SIZE))
    logger.info("Download completed in %f second(s)", timer() - start)

def purge_integral_db_cache():
    logger.warn("Purging the integral db cache directory '%s'...", INTEGRAL_DB_CACHE_DIR)
    shutil.rmtree(INTEGRAL_DB_CACHE_DIR)

def check_for_integral_db(beam_type_id, purge_cache=False):
    if purge_cache:
        purge_integral_db_cache()

    db_filename = get_integral_db_filename(beam_type_id)
    try:
        if tb.is_hdf5_file(db_filename):
            logger.info("Valid integral db file '%s' found in cache", db_filename)
            return
    except (IOError, tb.HDF5ExtError):
        pass

    logger.warn("'%s' is not a valid integral db file!", db_filename)
    download_integral_db(beam_type_id)

def open_integral_db(beam_type_id):
    return tb.open_file(get_integral_db_filename(beam_type_id))

def lookup_normalized_integral(integral_db, integral_id, m=None, t=None, v=None, n=None):
    integral_name = BaseIntegral.coerce(integral_id).name
    integral_table = integral_db.root._f_get_child(integral_name)

    d = locals()
    condition_vars = {
        var: d[var]
        for var in integral_table._v_attrs.used_variables_list
    }
    assert all(condition_vars.values())
    condition_str = ' & '.join(
        "(%s == %d)" % (k, v)
        for k, v in sorted(condition_vars.items())
    )

    matches = integral_table.read_where(condition_str)
    assert len(matches) == 1 # Ensure that our 'primary key' is indeed unique
    row = matches[0]

    normalized_integral = row['integral_float64']
    scale_factor = row['scale_factor']
    return normalized_integral, scale_factor

def get_scaled_integral(integral_db, integral_id, a, m=None, t=None, v=None, n=None):
    normalized_integral, scale_factor = lookup_normalized_integral(
        integral_db, integral_id, m, t, v, n
    )

    scale_by_value = a**scale_factor
    scaled_integral = normalized_integral * scale_by_value

    return scaled_integral
