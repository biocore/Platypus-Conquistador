#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2015--, platypus development team.
#
# Distributed under the terms of the BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from distutils.core import setup


classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries :: Application Frameworks
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ("Platypus Conquistador: Confirming specific taxonomic "
               "groups within your metagenomic samples.")

with open('README.rst') as f:
    long_description = f.read()

base = {"click", "scikit-bio >= 0.2.1, < 0.3.0"}
test = {"nose >= 0.10.1", "pep8", "flake8"}
all_deps = base | test

setup(name='platypus-conquistador',
      version='0.9.0',
      description=description,
      author='Antonio Gonzalez Pena',
      author_email='antgonza@gmail.com',
      maintainer='Antonio Gonzalez Pena',
      maintainer_email='antgonza@gmail.com',
      url='http://github.com/biocore/platypus',
      license='BSD',
      packages=['platypus'],
      scripts=['scripts/platypus'],
      install_requires=base,
      extras_require={'test': test, 'all': all_deps},
      package_data={},
      data_files={},
      long_description=long_description,
      classifiers=classifiers)
