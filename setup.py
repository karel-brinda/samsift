# see https://github.com/pypa/sampleproject

import setuptools

import os
import sys

if sys.version_info < (3,2):
    sys.exit('Minimum supported Python version is 3.2')

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Get the current version
exec(open("samsift/version.py").read())

setuptools.setup(
    name='samsift',

    version=VERSION,

    description='SAMsift - sift your alignments',

    long_description=long_description,

    url='https://github.com/karel-brinda/samsift',

    author='Karel Brinda',
    author_email='kbrinda@hsph.harvard.edu',

    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Unix',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    keywords='NGS SAM alignment',

    packages = ["samsift"],

    install_requires=[
            'pysam',
            'wheel',
        ],

    package_data={
        'samsift': [
            '*.py',
        ],
    },

   entry_points={
        'console_scripts': [
            'samsift = samsift.samsift:main',
            'samsift-norm-sam = samsift.samsift_norm_sam:main',
        ],
    },
)
