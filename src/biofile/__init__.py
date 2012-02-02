# coding: utf-8

# **********************************************************************
# file: __init__.py
#
# This Module is a dynamic import class to parse bio format file
# default module parse detected by file ext or file format
# **********************************************************************

import os


def _get_ftype(fname):
    ext = os.path.splitext(fname)[-1].replace('.', '')
    if ext == 'pls':
        return 'blat'
    elif ext in ('fa', 'fna', 'fsa', 'fasta'):
        return 'fasta'
    elif ext in ('fq', 'fastq'):
        return 'fastq'
    else:
        return ext


def _get_mod_from_ftype(ftype):
    if not ftype:
        raise ValueError('Unkown file type: {0}'.format(ftype))

    try:
        mod = __import__(__package__, fromlist=[ftype]) # loading ftype package
        return getattr(mod, ftype)                      # get ftype module
    except:
        raise


def parse(fname, ftype=None, **kwargs):
    """Parse fname, ftype is the module name should be used"""
    if not ftype:
        ftype = _get_ftype(fname)
    mod = _get_mod_from_ftype(ftype)
    if hasattr(mod, 'parse'):
        return mod.parse(fname, **kwargs) # return ftype's parse result
    else:
        raise ValueError('Module: {0} has not parse method'.format(ftype))


def read(fname, ftype=None, **kwargs):
    """read file """
    if not ftype:
        ftype = _get_ftype(fname)
    mod = _get_mod_from_ftype(ftype)
    if hasattr(mod, 'read'):
        return mod.read(fname, **kwargs) # return ftype's parse result
    else:
        raise ValueError('Module: {0} has not read method'.format(ftype))


def open(fname, ftype=None, **kwargs):
    """Open File with file name
    fname means file name(include file path, if no file path is provided,
    then current path is set as the path).
    ftype mean Format, the file type of file, Def. is using ext as ftype
    Def. is False, using parse
    """
    try:
        return parse(fname, ftype=ftype, **kwargs)
    except:
        raise

