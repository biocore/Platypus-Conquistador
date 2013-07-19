Platypus Conquistador installation notes
===========================

Platypus Conquistador is a python package that relies on [QIIME](http://www.qiime.org) and [qcli](https://github.com/bipy/qcli). These packages must be installed prior executing the `setup.py` script.

Installation
============

To perform a global installation of Emperor, execute the following command from a terminal session:

    python setup.py install

If you do not want to do a global installation, you will have to add the Emperor scripts and libraries to the `PATH` and `PYTHONPATH` environment variables. To add these variables to your `.bash_profile` issue the following terminal commands:

``` bash
echo "export PATH=$HOME/platypus/scripts/:$PATH" >> ~/.bash_profile
echo "export PYTHONPATH=$HOME/platypus/:PYTHONPATH" >> ~/.bash_profile
python setup.py install --install-scripts=~/platypus/scripts/ --install-purelib=~/platypus/ --install-lib=~/platypus/
```

To test for a correct installation, open a new terminal session and issue the following command to see the help of `platypus_split_db.py`:

    platypus_split_db.py -h