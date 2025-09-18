# see https://github.com/pypa/sampleproject

import setuptools

import os
import sys

if sys.version_info < (3, 2):
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
    author_email='karel.brinda@inria.fr',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Unix',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='NGS SAM alignment',
    packages=["samsift"],
    install_requires=['pysam'],
    python_requires='>=3.8',
    package_data={
        'samsift': [
            '*.py',
        ],
    },
    long_description_content_type='text/x-rst',
    entry_points={
        'console_scripts': [
            'samsift = samsift.samsift:main',
            'samsift-norm-sam = samsift.samsift_norm_sam:main',
        ],
    },
)
