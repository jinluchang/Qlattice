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


SelectedFieldBase
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^

.. autosummary::
   :recursive:
   :toctree: generated

   PointsSelection

SelectedPointsBase
^^^^^^^^^^^^^^^^^

.. autosummary::
   :recursive:
   :toctree: generated

   SelectedPoints
   SelectedPointsBase
   SelectedPointsBase.save
   SelectedPointsBase.load
   SelectedPointsBase.save_complex
   SelectedPointsBase.load_complex
   SelectedPointsBase.to_numpy
   SelectedPointsBase.from_numpy

