About
=====

Console app and Python API implementing a generalization of eigenvalue problem
within the harmonic coupled finite strip method, used for parametric modeling
of static and dynamic inelastic buckling, free vibration, damage and failure in
prismatic shell structures.

This work is a part of the investigation within the research project
[ON174027]_, supported by the Ministry for Science and Technology, Republic of
Serbia. This support is gratefully acknowledged.

References
----------

.. [Milasinovic1997]
   Milašinović, D.D. "The Finite Strip Method in Computational Mechanics".
   Faculties of Civil Engineering: University of Novi Sad, Technical University
   of Budapest and University of Belgrade: Subotica, Budapest, Belgrade. (1997)
.. [ON174027]
   "Computational Mechanics in Structural Engineering"

Installation
============

To install fsm_eigenvalue run::

    $ pip install fsm_eigenvalue

Console app usage
=================

Quick start::

    $ fsm_eigenvalue <filename>

Show help::

    $ fsm_eigenvalue --help

Python API usage
================

Quick start::

    >>> import logging
    >>> logging.basicConfig(level=logging.DEBUG)

    >>> from fsm_eigenvalue.compute import parameter_sweep
    >>> from fsm_eigenvalue.load import load_data_from
    >>> from fsm_eigenvalue.store import store_results_to

    >>> data_file = 'examples/data-files/barbero-viscoelastic.yaml'
    >>> results_file = data_file.replace('.yaml', '.hdf5')

    >>> beam_type_id, search_space, nodal_graph, strip_data, materials, astiff_shape = load_data_from(data_file)
    >>> with parameter_sweep(beam_type_id, search_space, strip_data, materials, astiff_shape) as results_iterator:
    ...     store_results_to(results_file, data_file, search_space, astiff_shape, results_iterator)

Contribute
==========

If you find any bugs, or wish to propose new features `please let us know`_.

If you'd like to contribute, simply fork `the repository`_, commit your changes
and send a pull request. Make sure you add yourself to `AUTHORS`_.

.. _`please let us know`: https://github.com/petarmaric/fsm_eigenvalue/issues/new
.. _`the repository`: https://github.com/petarmaric/fsm_eigenvalue
.. _`AUTHORS`: https://github.com/petarmaric/fsm_eigenvalue/blob/master/AUTHORS
