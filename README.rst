Platypus Conquistador
=====================

|Build Status| |Coverage Status|

Platypus Conquistador is a bioinformatic command line tool for the confirmation
of the presence/absence of a specific taxon (or set of taxa) in environmental
shotgun sequence reads. Platypus uses two reference databases: an interest
database (containing sequences of interest) and an other database (containing
all reference library sequences). Platypus uses these databases to generate
data showing the number of sequences to those that are closest to the sequences
of interest, based on similarity, and their sequence match; this information is
also provided for the rest (other) reference database.

To start using Platypus, please refer to the `installation notes <https://github.com/biocore/Platypus-Conquistador/blob/master/INSTALL.md>`__;
also check the `tutorial <https://github.com/biocore/Platypus-Conquistador/blob/master/TUTORIAL.md>`__.

.. |Build Status| image:: https://travis-ci.org/biocore/Platypus-Conquistador.svg
   :target: https://travis-ci.org/biocore/Platypus-Conquistador
.. |Coverage Status| image:: https://coveralls.io/repos/biocore/Platypus-Conquistador/badge.svg
   :target: https://coveralls.io/r/biocore/Platypus-Conquistador
