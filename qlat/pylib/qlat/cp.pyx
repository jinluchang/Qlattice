# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

def qremove_info(path):
    return cc.qremove_info(path)

def qremove_all_info(path):
    return cc.qremove_all_info(path)
