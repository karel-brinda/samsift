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

       git clone http://github.com/karel-brinda/samsift
       cd samsift
       export PATH=$(pwd)/samsift:$PATH
       # filter: score >94, save as filtered.bam
       samsift -i tests/test.bam -o filtered.bam -f 'AS>94'
       # filter: unaligned reads
       samsift -i tests/test.bam -f 'FLAG & 0x04'
       # filter: aligned reads
       samsift -i tests/test.bam -f 'not(FLAG & 0x04)'
       # filter: sequences containing ACCAGAGGAT
       samsift -i tests/test.bam -f 'SEQ.find("ACCAGAGGAT")!=-1'
       # annotation: add tags 'ln' with sequence length and 'ab' with average base quality
       samsift -i tests/test.bam -c 'ln=len(SEQ);ab=1.0*sum(QUAL)/ln'


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
	Version: 0.1.0
	Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

	Usage:   samsift.py [-h] [-v] [-i FILE] [-o FILE] [-f PY_EXPR] [-c PY_CODE] [-d PY_EXPR] [-t PY_EXPR] [-m STR]

	Options:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -i FILE               input SAM/BAM file [-]
	  -o FILE               output SAM/BAM file [-]
	  -f PY_EXPR            filter [True]
	  -c PY_CODE            code to be executed (e.g., assigning new tags) [None]
	  -d PY_EXPR            debugging expression to print [None]
	  -t PY_EXPR            debugging trigger [True]
	  -m STR                mode: strict (stop upon first error)
	                              nonstop-keep (keep alignments causing errors)
	                              nonstop-remove (remove alignments causing errors) [strict]


.. USAGE-END

Algorithm
---------

.. code-block:: python

        for ALIGNMENT in ALIGNMENTS:
                if eval(DEBUG_TRIGER):
                        print(eval(DEBUG_EXPR))
                if eval(FILTER):
                        exec(CODE)
                        print(ALIGNMENT)


**Python expression.** All expressions should be valid `Python 3 expressions
<https://docs.python.org/3/reference/expressions.html>`_. They are evaluated
using the `eval <https://docs.python.org/3/library/functions.html#eval>`_
function.

**Python code.** Code is executed using the `exec
<https://docs.python.org/3/library/functions.html#exec>`_ function.

**SAM fields.** All Python expressions and code can access variables mirroring
all the fields from the alignment section of the `SAM specification
<https://samtools.github.io/hts-specs/SAMv1.pdf>`_, i.e., `QNAME`, `FLAG`,
`RNAME`, `POS` (1-based), `MAPQ`, `CIGAR`, `RNEXT`, `PNEXT`, `TLEN`, `SEQ`,
and `QUAL`.  For instance, we can filter reads, keeping only those with `POS`
smaller than 10000, by

.. code-block:: bash

        samsift -i tests/test.bam -f 'POS<10000'


The PySAM representation of the current alignment (class `pysam.AlignedSegment
<http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment>`_) is
available through the variable `a`. Therefore, the previous example is equivalent
to

.. code-block:: bash

        samsift -i tests/test.bam -f 'a.reference_start+1<10000'


**SAM tags.** All SAM tags are translated to variables with the same name. For
instance, if alignment scores are provided through the `AS` tag (as defined in
the `Sequence Alignment/Map Optional Fields Specification
<https://samtools.github.io/hts-specs/SAMtags.pdf>`_), then alignments with
score smaller or equal to the sequence length can be removed using

.. code-block:: bash

        samsift -i tests/test.bam -f 'AS>len(SEQ)'

If `CODE` is provided, all two-letter variables are back-translated to tags.
For instance, a tag `ab` carrying the average base quality can be added by

.. code-block:: bash

        samsift -i tests/test.bam -c 'ab=1.0*sum(QUAL)/len(QUAL)'


Similar programs
----------------

* `samtools view <http://www.htslib.org/doc/samtools.html>`_ can filter alignments based on FLAGS, read group tags, and CIGAR strings.
* `sambamba view <http://lomereiter.github.io/sambamba/docs/sambamba-view.html>`_ supports, in addition to SAMtools, filtration using `simple perl expression <https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax>`_. However, it's not possible to compare different tags.
* `bamPals <https://github.com/zeeev/bamPals>`_ adds tags XB, XE, XP and XL.
* `SamJavascript <http://lindenb.github.io/jvarkit/SamJavascript.html>`_ can filter alignments using JavaScript expressions.


Author
------

`Karel Brinda <http://brinda.cz>`_ <kbrinda@hsph.harvard.edu>
