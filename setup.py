#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: setup.py
#
# pyngs package install script
# **********************************************************************

from distutils.core import setup, Extension
import src as pyngs
import os


# name becomes the name of the installation tarball.
name = 'pyngs'

author = pyngs.__author__

author_email = pyngs.__author_email__

# version gets appended to the end of the tarball name.
__version__ = pyngs.__version__

# package_dir option to tell the Distutils about "root package"
package_dir = {'pyngs': 'src'}

# packages are A list naming all the packages to include.
packages = ['pyngs',                    # the root package
            'pyngs.biofile',            # parse all kinds of bio format file
            'pyngs.lib',                # parse all kinds of base modules
            ]

# loading scripts
script_dir = 'src/scripts'
scripts = [os.path.join(script_dir, script).format(script)
           for script in os.listdir(script_dir)]

# "description" may be a short, summary description of the package
description = 'Next Generation Sequencing Toolkits in Python'

# "long_description" may be valid ReStructuredText which will be turned into
# HTML when displayed
long_description = open('README.txt', 'r').read()

# "platforms" a list of platforms
# platforms = ['Linux', 'Windwos', 'Mac']
platforms = 'any'

# license for the package
license_ = 'MIT'

calign_ext = Extension('pyngs.lib.calign', sources=['src/lib/calignmodule.c'])


if __name__ == '__main__':
    setup(
        name=name,                      # Next Generation Sequencing Toolkits
        version=__version__,            # version alpha
        author=author,
        author_email=author_email,
        url='https://github.com/azerxu/pyngs', # github project url
        package_dir=package_dir,
        packages=packages,
        ext_modules=[calign_ext],
        description=description,
        long_description=long_description,
        platforms=platforms,
        scripts=scripts,
        license = license_
        )

