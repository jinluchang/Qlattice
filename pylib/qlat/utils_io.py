import cqlat as c

import os

from cqlat import qremove, qremove_info
from cqlat import qremove_all, qremove_all_info
from cqlat import qmkdir, qmkdir_info, qmkdir_sync_node
from cqlat import obtain_lock, release_lock
from cqlat import does_file_exist, does_file_exist_sync_node
from cqlat import is_directory, is_directory_sync_node
from cqlat import qtouch, qtouch_info
from cqlat import qappend, qappend_info
from cqlat import qrename, qrename_info
from cqlat import qcat, qcat_sync_node
from cqlat import qls, qls_sync_node
from cqlat import qload_datatable, qload_datatable_sync_node
from cqlat import check_time_limit, check_stop

from qlat.mpi import *

@timer
def qmkdirs(path):
    os.makedirs(path, exist_ok=True)

@timer
def qmkdirs_info(path):
    if get_id_node() == 0:
        displayln(f"qmkdirs_info: '{path}'.")
        qmkdirs(path)

@timer
def mk_dirs(path):
    os.makedirs(path, exist_ok=True)

@timer
def mk_dirs_info(path):
    if get_id_node() == 0:
        displayln(f"mk_dirs_info: '{path}'.")
        mk_dirs(path)

@timer
def mk_file_dirs(fn):
    path = os.path.dirname(fn)
    os.makedirs(path, exist_ok=True)

@timer
def mk_file_dirs_info(path):
    if get_id_node() == 0:
        displayln(f"mk_file_dirs_info: '{path}'.")
        mk_file_dirs(path)
