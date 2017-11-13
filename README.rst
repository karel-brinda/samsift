SAMsift
=======

.. image:: https://travis-ci.org/karel-brinda/samsift.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/samsift

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
	:target: https://anaconda.org/bioconda/samsift

.. image:: https://badge.fury.io/py/samsift.svg
        :target: https://badge.fury.io/py/samsift

SAMsift is a program for advanced filtering and tagging of SAM/BAM alignments
using Python expressions.


Getting started
---------------

.. code-block:: bash

       # clone this repo and add it to PATH
       git clone http://github.com/karel-brinda/samsift
       cd samsift
       export PATH=$(pwd)/samsift:$PATH

       # filtering: keep only alignments with score >94, save them as filtered.bam
       samsift -i tests/test.bam -o filtered.bam -f 'AS>94'
       # filtering: keep only unaligned reads
       samsift -i tests/test.bam -f 'FLAG & 0x04'
       # filtering: keep only aligned reads
       samsift -i tests/test.bam -f 'not(FLAG & 0x04)'
       # filtering: keep only sequences containing ACCAGAGGAT
       samsift -i tests/test.bam -f 'SEQ.find("ACCAGAGGAT")!=-1'
       # filtering: sample alignments with 25% rate
       samsift -i tests/test.bam -f 'random.random()<0.25'
       # filtering: sample alignments with 25% rate with a fixed RNG seed
       samsift -i tests/test.bam -f 'random.random()<0.25' -0 'random.seed(42)'
       # tagging: add tags 'ln' with sequence length and 'ab' with average base quality
       samsift -i tests/test.bam -c 'ln=len(SEQ);ab=1.0*sum(QUALa)/ln'
       # tagging: add a tag 'ii' with the number of the current alignment
       samsift -i tests/test.bam -0 'i=0' -c 'i+=1;ii=i'
       # updating: removing sequences and base qualities
       samsift -i ./tests/test.bam -c 'a.query_sequence=""'
       # updating: switching all reads to unaligned
       samsift -i ./tests/test.bam -c 'a.flag|=0x4;a.reference_start=-1;a.cigarstring="";a.reference_id=-1;a.mapping_quality=0'


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

.. USAGE-BEGIN

.. code-block::

	Program: samsift (advanced filtering and tagging of SAM/BAM alignments using Python expressions)
	Version: 0.2.3
	Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

	Usage:   samsift.py [-i FILE] [-o FILE] [-f PY_EXPR] [-c PY_CODE] [-m STR]
	                    [-0 PY_CODE] [-d PY_EXPR] [-t PY_EXPR]

	Basic options:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -i FILE               input SAM/BAM file [-]
	  -o FILE               output SAM/BAM file [-]
	  -f PY_EXPR            filter [True]
	  -c PY_CODE            code to be executed (e.g., assigning new tags) [None]
	  -m STR                mode: strict (stop on first error)
	                              nonstop-keep (keep alignments causing errors)
	                              nonstop-remove (remove alignments causing errors) [strict]

	Advanced options:
	  -0 PY_CODE            initialization [None]
	  -d PY_EXPR            debugging expression to print [None]
	  -t PY_EXPR            debugging trigger [True]


.. USAGE-END

Algorithm
---------

.. code-block:: python

        exec(INITIALIZATION)
        for ALIGNMENT in ALIGNMENTS:
                if eval(DEBUG_TRIGER):
                        print(eval(DEBUG_EXPR))
                if eval(FILTER):
                        exec(CODE)
                        print(ALIGNMENT)


**Python expressions and code.** All expressions and code should be valid with
respect to `Python 3 <https://docs.python.org/3/>`_. Expressions are evaluated
using the `eval <https://docs.python.org/3/library/functions.html#eval>`_
function and code is executed using the `exec
<https://docs.python.org/3/library/functions.html#exec>`_ function.
Initialization can be used for importing Python modules, setting global
variables (e.g., counters) or loading data from disk. Some modules (e.g.,
``random``) are loaded without an explicit request.

*Example* (printing all alignments):

.. code-block:: bash

        samsift -i tests/test.bam -f 'True'

**SAM fields.** Expressions and code can access variables mirroring the fields
from the alignment section of the `SAM specification
<https://samtools.github.io/hts-specs/SAMv1.pdf>`_, i.e., ``QNAME``, ``FLAG``,
``RNAME``, ``POS`` (1-based), ``MAPQ``, ``CIGAR``, ``RNEXT``, ``PNEXT``,
``TLEN``, ``SEQ``, and ``QUAL``. Several additional variables are defined to
simply accessing some useful information: ``QUALa`` stores the base qualities
as an integer array;  ``SEQs``, ``QUALs``, ``QUALsa`` skip soft-clipped bases;
and ``RNAMEi`` and ``RNEXTi`` store the reference ids as integers.

*Example* (keeping only the alignments with leftmost position <= 10000):

.. code-block:: bash

        samsift -i tests/test.bam -f 'POS<=10000'


SAMsift internally uses the `PySam <http://pysam.readthedocs.io/>`_ library and
the representation of the current alignment (an instance of the class
`pysam.AlignedSegment
<http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment>`_) is
available as a variable ``a``. Therefore, the previous example is equivalent to

.. code-block:: bash

        samsift -i tests/test.bam -f 'a.reference_start+1<=10000'


The ``a`` variable can also be used for modifying the current alignment record.

*Example* (removing the sequence and the bases from every record):

.. code-block:: bash

        samsift -i ./tests/test.bam -c 'a.query_sequence=""'


**SAM tags.** Every SAM tag is translated to a variable with the same name.

*Example* (removing alignments with a score smaller or equal to the sequence length):

.. code-block:: bash

        samsift -i tests/test.bam -f 'AS>len(SEQ)'

If ``CODE`` is provided, all two-letter variables are back-translated after its execution to tags.

*Example* (adding a tag ``ab`` carrying the average base quality):

.. code-block:: bash

        samsift -i tests/test.bam -c 'ab=1.0*sum(QUALa)/len(QUALa)'

**Errors.** If an error occurs during an evalution of an expression or an
execution of a code (e.g., due to accessing an undefined tag), then SAMsift
behavior depends on the specified mode (``-m``).  With the strict mode (``-m
strict``, default), SAMsift will immediately interrupt the computation and
report an error.  With the ``-m nonstop-keep`` option, SAMsift will continue
processing the alignments while keeping the error-causing alignments in the
output.  With the ``-m nonstop-remove`` option, all error-causing alignments
are skipped and ommited from the output.


Similar programs
----------------

* `samtools view <http://www.htslib.org/doc/samtools.html>`_ can filter alignments based on FLAGS, read group tags, and CIGAR strings.
* `sambamba view <http://lomereiter.github.io/sambamba/docs/sambamba-view.html>`_ supports, in addition to SAMtools, filtration using `simple Perl expressions <https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax>`_. However, it's not possible to use floats or compare different tags.
* `bamPals <https://github.com/zeeev/bamPals>`_ adds tags XB, XE, XP and XL.
* `SamJavascript <http://lindenb.github.io/jvarkit/SamJavascript.html>`_ can filter alignments using JavaScript expressions.


Issues
------

Please use `Github issues <https://github.com/karel-brinda/samsift/issues>`_.


Changelog
---------

See `Releases <https://github.com/karel-brinda/samsift/releases>`_.


Licence
-------

`MIT <https://github.com/karel-brinda/samsift/blob/master/LICENSE>`_


Author
------

`Karel Brinda <http://brinda.cz>`_ <kbrinda@hsph.harvard.edu>
