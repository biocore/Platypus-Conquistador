#!/usr/bin/env python

from distutils.core import setup
from glob import glob

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011-2013, The Platypus Project"
__credits__ = ["Antonio Gonzalez Pena",]
__license__ = "GPL"
__version__ = "0.0.8-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

# without any of these platypus will not function correctly
required_python_modules = ['qiime'] # qcli will be added later
unavailable_dependencies = []

for module in required_python_modules:
    try:
        exec 'import %s' % module
    except ImportError:
        unavailable_dependencies.append(module)

if unavailable_dependencies:
    print ('Cannot find the following python package(s): %s. Check your '
        'PYTHONPATH environment variable and/or site-packages folder.' %
        ', '.join(unavailable_dependencies))
    exit(1)

# slightly modified from the biom-format setup.py script
qiime_version = tuple(map(int, qiime.__version__.replace('-dev','').split('.')))
if qiime_version < (1, 7, 0):
    print ('The minimum required version of the QIIME libraries is 1.7.0 '
        'please update your version accordingly (your current version %s).' %
        qiime.__version__)
    exit(2)

long_description = """Platypus Conquistador: Confirming specific taxonomic groups within your samples.
"""

setup(name='platytpus',
        version=__version__,
        description='Platypus Conquistador',
        author="Antonio Gonzalez Pena",
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url='http://github.com/qiime/platypus',
        packages=['platypus'],
        scripts=glob('scripts/*py'),
        package_data={},
        data_files={},
        long_description=long_description)

