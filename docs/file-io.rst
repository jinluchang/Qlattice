File IO
-------

.. module:: qlat_utils
.. module:: qlat

QAR
^^^

.. toctree::
   :maxdepth: 1

   contents/qar-format.md

.. autosummary::
   :toctree: generated

   list_qar
   qar_create_info
   qar_extract_info
   get_qar_multi_vol_max_size
   qcat
   qcat_bytes
   qcopy_file_info
   does_file_exist
   does_file_exist_qar
   does_regular_file_exist_qar

LatData
^^^^^^^

Multi-dimension array data structure for IO.

.. toctree::
   :maxdepth: 1

   contents/latio-format.md

.. autosummary::
   :recursive:
   :toctree: generated

   LatData
   LatData.info
   LatData.set_info
   LatData.to_numpy
   LatData.from_numpy
   LatData.load
   LatData.save
   mk_lat_data
   load_lat_data

Example: ``examples-py/lat-io.py``

.. literalinclude:: ../examples-py/lat-io.py
