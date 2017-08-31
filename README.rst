SAMsift
=======

SAMsift is a program for advanced filtering and tagging of SAM/BAM alignments
using Python expressions.


Getting started
---------------

.. code-block:: bash

       git clone http://github.com/karel-brinda/samsift
       cd samsift
       # keep only alignments with alignment score >94
       samsift/samsift -i tests/test.bam -o filtered.sam -f 'AS>94'
       # add tags 'ln' with sequence length and 'ab' with average base quality
       samsift/samsift -i tests/test.bam -o with_ln_ab.sam -c 'ln=len(SEQ);ab=1.0*sum(QUAL)/ln'


Installation
------------

**Using Bioconda:**

.. code-block:: bash

        # add all necessary Bioconda channels
        conda config --add channels defaults
        conda config --add channels conda-forge
        conda config --add channels bioconda

        # install samsift
        conda install samsift


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

        Program: samsift (advanced filtering and tagging of SAM/BAM alignments using Python expressions)
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


All Python expressions can access variables mirroring the fields from the
alignment section of the `SAM specification
<https://samtools.github.io/hts-specs/SAMv1.pdf>`_, i.e., `QNAME`, `FLAG`,
`RNAME`, `POS` (1-based), `MAPQ`, `CIGAR` , `RNEXT`, `PNEXT`, `TLEN`, `SEQ`,
and `QUAL`.  For instance, keeping only reads with `POS` smaller than  10000
can be done by

.. code-block:: bash

        samsift -i tests/test.bam -f 'POS<10000'


The PySAM representation of current alignment (class `pysam.AlignedSegment
<http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment>`_) is
available through variable `a`. Therefore, the previous example is equivalent
to

.. code-block:: bash

        samsift -i tests/test.bam -f 'a.reference_starts+1<10000'


All SAM tags are translated to variables with equal name. For instance, if
alignment score is provided through the `AS` tag (as it is defined in the
`Sequence Alignment/Map Optional Fields Specification
<https://samtools.github.io/hts-specs/SAMtags.pdf>`_), then alignments with
score smaller or equal to the sequence length can be removed using

.. code-block:: bash

        samsift -i tests/test.bam -f 'AS>len(SEQ)'

If `CODE` is provided, all two-letter variables are back-translated to tags.
For instance, a tag `ab` carrying the average base quality can be added by

.. code-block:: bash

        samsift -i tests/test.bam -c 'ab=1.0*sum(QUAL)/ln'


Similar programs
----------------

* `samtools view <http://www.htslib.org/doc/samtools.html>`_ can filter alignments based on FLAGS, read group tags, and CIGAR strings.
* `sambamba view <http://lomereiter.github.io/sambamba/docs/sambamba-view.html>`_ supports, in addition to SAMtools, filtration using `simple perl expression <https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax>`_. However, it's not possible to compare different tags.

Author
------

`Karel Brinda <http://brinda.cz>`_ <kbrinda@hsph.harvard.edu>
