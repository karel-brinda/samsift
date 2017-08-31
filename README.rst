SAMsift
=======

Getting started
---------------

.. code-block:: bash

       git clone http://github.com/karel-brinda/samsift
       cd samsift
       # filter alignments with with alignment score >94
       ./samsift.sh -i tests/test.bam -o filtered.sam -f 'AS>94'
       # add tag LN with sequence length
       ./samsift.sh -i tests/test.bam -o with_ln.sam -c 'LN=len(SEQ)'


Introduction
------------

SAMsift is a program for filtering and enriching alignments in SAM/BAM format.

Installation
------------

**Using PIP from PyPI:**

.. code-block:: bash

   pip install --upgrade samsift

**Using PIP from Github:**

.. code-block:: bash

   pip install --upgrade git+https://github.com/karel-brinda/samsift


Command-line parameters
-----------------------

.. code-block::

        usage: samsift.py [-h] [-v] [-i file] [-o file] [-f py_expr] [-c py_code]
                          [-d py_expr] [-t py_expr]

        Program: samsift (sift and enrich SAM/BAM alignments using Python expressions)
        Version: 0.1.0
        Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

        optional arguments:
          -h, --help     show this help message and exit
          -v, --version  show program's version number and exit
          -i file        input SAM/BAM file [-]
          -o file        output SAM/BAM file [-]
          -f py_expr     filter [True]
          -c py_code     code to be executed (e.g., assigning new tags) [None]
          -d py_expr     debugging expression to print [None]
          -t py_expr     debugging trigger [True]


Algorithm
---------

.. code-block:: python

        for ALIGNMENT in ALIGNMENTS:
                if eval(DEBUG_TRIGER):
                        print(eval(DEBUG_EXPR))
                if eval(FILTER):
                        exec(CODE)
                        print(ALIGNMENT)



Examples
--------

Author
------

`Karel Brinda <http://brinda.cz>`_ <kbrinda@hsph.harvard.edu>
