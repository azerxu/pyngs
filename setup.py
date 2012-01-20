#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: setup.py
# **********************************************************************

from distutils.core import setup

setup(name='pyngs',
      version='0.1',
      description='Python Next Generation Sequencing Package',
      long_description='Unkown',
      author='azer xu',
      author_email='azer.xu@gmail.com',
      url='http://www.python.org/sigs/distutils-sig/',
      #packages=['distutils', 'distutils.command'],
      package_dir = {'pyngs': 'src'},
      packages=['pyngs', 'pyngs.biofile'],
      license='MIT',
      platforms='any',
     )




