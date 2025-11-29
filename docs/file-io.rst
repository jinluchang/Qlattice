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

   QFile
   QFile.copy
   QFile.path
   QFile.mode
   QFile.close
   QFile.null
   QFile.eof
   QFile.tell
   QFile.flush
   QFile.seek_set
   QFile.seek_end
   QFile.seek_cur
   QFile.content
   QFile.content_bytes
   QFile.size
   QFile.remaining_size
   QFile.getline
   QFile.getlines
   QFile.qcat
   QFile.qcat_bytes
   QFile.write
   QFile.compute_crc32
   open_qfile
   open_qfile_str
   QarFile
   QarFile.path
   QarFile.mode
   QarFile.close
   QarFile.null
   QarFile.flush
   QarFile.list
   QarFile.has_regular_file
   QarFile.has
   QarFile.__contains__
   QarFile.read
   QarFile.read_data
   QarFile.read_data_bytes
   QarFile.read_info
   QarFile.verify_index
   QarFile.write
   QarFile.show_index
   QarFile.read_index
   QarFile.index_size
   QarFile.index_size_saved
   QarFile.save_index
   open_qar
   open_qar_info
   get_qar_multi_vol_max_size
   set_qar_multi_vol_max_size
   clean_up_qfile_map
   show_all_qfile
   properly_truncate_qar_file
   does_regular_file_exist_qar
   does_file_exist_qar
   qar_build_index
   qar_create
   qar_extract
   qcopy_file
   list_qar
   qcat
   qcat_bytes
   qtouch
   qappend
   qload_datatable
   compute_crc32
   qar_build_index_info
   qar_create_info
   qar_extract_info
   qcopy_file_info
   qtouch_info
   qappend_info
   check_all_files_crc32_info
   does_regular_file_exist_qar_sync_node
   does_file_exist_qar_sync_node
   qar_create_sync_node
   qar_extract_sync_node
   qcopy_file_sync_node
   qcat_sync_node
   qcat_bytes_sync_node
   qload_datatable_sync_node

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
   LatData.to_list
   LatData.from_list
   LatData.copy
   LatData.load
   LatData.save
   LatData.load_str
   LatData.save_str
   LatData.glb_sum
   LatData.glb_sum_in_place
   LatData.bcast
   LatData.set_zero
   LatData.is_complex
   LatDataRealF
   LatDataInt
   LatDataLong
   mk_lat_data
   mk_lat_data_real_f
   mk_lat_data_int
   mk_lat_data_long
   load_lat_data
   load_lat_data_real_f
   load_lat_data_int
   load_lat_data_long

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

