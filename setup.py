#! /usr/bin/env python
# coding: utf-8

# **********************************************************************
# file: setup.py
#
# ngs package install script
# **********************************************************************

from distutils.core import setup
import src as pyngs

# name becomes the name of the installation tarball.
name = 'pyngs'

author = pyngs.author

author_email = pyngs.author_email

# version gets appended to the end of the tarball name.
__version__ = pyngs.__version__

# package_dir option to tell the Distutils about "root package"
package_dir = {'pyngs': 'src'}

# packages are A list naming all the packages to include.
packages = ['pyngs',                    # the root package
            'pyngs.biofile',            # parse all kinds of bio format file
            ]

# "description" may be a short, summary description of the package
description = 'Next Generation Sequencing Toolkits in Python'


# "long_description" may be valid ReStructuredText which will be turned into
# HTML when displayed
long_description = open('readme.rst', 'r').read()


# "platforms" a list of platforms
# platforms = ['Linux', 'Windwos', 'Mac']
platforms = 'any'

# license for the package
license_ = 'MIT'


if __name__ == '__main__':
    setup(
        name=name,                      # Next Generation Sequencing Toolkits
        version=__version__,            # version alpha
        author=author,
        author_email=author_email,
        url='https://github.com/azerxu/pyngs', # github project url
        description=description,
        long_description=long_description,
        platforms=platforms,
        license = license_
        )

