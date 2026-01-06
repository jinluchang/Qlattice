Jackknife
---------

.. module:: qlat_utils

Jackknife method
^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   contents/jackknife-random.md

Jackknife implementation
^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated

   g_mk_jk
   g_mk_jk_val
   g_jk_avg
   g_jk_err
   g_jk_avg_err
   g_jk_avg_err_arr
   g_jk_size
   g_jk_blocking_func

   get_jk_state
   set_jk_state

Utilities
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated

   average
   avg_err
   err_sum
   block_data
   fsqr
   fsqrt

   jackknife
   jk_avg
   jk_err
   jk_avg_err

   sjackknife
   sjk_mk_jk_val
   sjk_avg
   sjk_err
   sjk_avg_err

   rjackknife
   rjk_mk_jk_val
   rjk_avg
   rjk_err
   rjk_avg_err


Example for the Jackknife-bootstrap hybrid method (described in the Jackknife method section): ``examples-py/jackknife-random.py``

.. literalinclude:: ../examples-py/jackknife-random.py

Example for the conventional Super-Jackknife method: ``examples-py/jackknife-super.py``

.. literalinclude:: ../examples-py/jackknife-super.py

Example for a variant of the conventional Super-Jackknife method: ``examples-py/jackknife-super-hash.py``

.. literalinclude:: ../examples-py/jackknife-super-hash.py
