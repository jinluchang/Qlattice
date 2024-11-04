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

Pickle
^^^^^^

.. autosummary::
   :toctree: generated

   save_pickle_obj
   load_pickle_obj
   pickle_cache_call

LatData
^^^^^^^

Multi-dimension array data structure for IO.

File format description
.......................

.. code::

    FILE-HEADER
    BINARY-DATA

``FILE-HEADER`` example::

    #!/usr/bin/env lat-io-glimpse
    data_size
    128
    ndim: 3
    "tsep"[4]: "0" "1" "2" "3"
    "op"[2]: "0" "1"
    "re-im"[2]: "re" "im"
    crc32: 77A655DB
    END_HEADER

``BINARY-DATA`` description:

Consist of the data stored in double precision (little endian) in sequential order as a standard C multidimensional array.

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

C++ example: ``examples-cpp/latio/main.cpp``

.. literalinclude:: ../examples-cpp/latio/main.cpp

C++ code: ``qlat-utils/include/qlat-utils/lat-io.h``

.. literalinclude:: ../qlat-utils/qlat_utils/include/qlat-utils/lat-io.h

Gauge Field
^^^^^^^^^^^

.. autosummary::
   :recursive:
   :toctree: generated

   GaugeField
   GaugeField.save
   GaugeField.load

Gauge Transform
^^^^^^^^^^^^^^^

.. autosummary::
   :recursive:
   :toctree: generated

   GaugeTransform
   GaugeTransform.save
   GaugeTransform.load
   GaugeTransform.save_cps
   GaugeTransform.load_cps

FieldBase
^^^^^^^^^

Support ``np.asarray(f)``.

.. autosummary::
   :recursive:
   :toctree: generated

   Field
   FieldBase
   FieldBase.save_direct
   FieldBase.load_direct
   FieldBase.save_64
   FieldBase.load_64
   FieldBase.save_double
   FieldBase.load_double
   FieldBase.save_float_from_double
   FieldBase.load_double_from_float
   FieldBase.float_from_double
   FieldBase.double_from_float
   FieldBase.to_from_endianness
   FieldBase.as_field
   FieldBase.from_field

FieldSelection
^^^^^^^^^^^^^^

.. autosummary::
   :recursive:
   :toctree: generated

   FieldSelection
   FieldSelection.save
   FieldSelection.load
   FieldSelection.to_psel
   FieldSelection.to_psel_local

SelectedFieldBase
^^^^^^^^^^^^^^^^^

Support ``np.asarray(sf)``.

.. autosummary::
   :recursive:
   :toctree: generated

   SelectedField
   SelectedFieldBase
   SelectedFieldBase.save_direct
   SelectedFieldBase.load_direct
   SelectedFieldBase.save_64
   SelectedFieldBase.load_64
   SelectedFieldBase.save_double
   SelectedFieldBase.load_double
   SelectedFieldBase.save_float_from_double
   SelectedFieldBase.load_double_from_float
   SelectedFieldBase.float_from_double
   SelectedFieldBase.double_from_float
   SelectedFieldBase.to_from_endianness

PointsSelection
^^^^^^^^^^^^^^^

Support ``np.asarray(psel)``.

.. autosummary::
   :recursive:
   :toctree: generated

   PointsSelection
   PointsSelection.save
   PointsSelection.load
   PointsSelection.xg_arr

SelectedPointsBase
^^^^^^^^^^^^^^^^^^

Support ``np.asarray(sp)``.

.. autosummary::
   :recursive:
   :toctree: generated

   SelectedPoints
   SelectedPointsBase
   SelectedPointsBase.save_str
   SelectedPointsBase.load_str
   SelectedPointsBase.to_numpy
   SelectedPointsBase.from_numpy
   SelectedPointsRealD.save
   SelectedPointsRealD.load
   SelectedPointsRealD.to_lat_data
   SelectedPointsRealD.from_lat_data

