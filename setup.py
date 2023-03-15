import sys
from pathlib import Path

if sys.version_info < (3, 6):
    sys.exit('scHiCPTR requires Python >= 3.6')

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize



setup(
    name="scKTLD", ##project name, for install and uninstall
    version="1.2",
    author="double-anonymous",
    author_email="double-anonymous",

    packages=find_packages(where="src"),
    package_dir={"":"src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    ext_modules=cythonize([Extension("scKTLD.cpd_utils._detection.ekcpd",sources=["src/scKTLD/cpd_utils/_detection/ekcpd.pyx", "src/scKTLD/cpd_utils/_detection/ekcpd_computation.c","src/scKTLD/cpd_utils/_detection/ekcpd_pelt_computation.c","src/scKTLD/cpd_utils/_detection/kernels.c", ],),Extension(
        "scKTLD.cpd_utils.utils._utils.convert_path_matrix",sources=[ "src/scKTLD/cpd_utils/utils/_utils/convert_path_matrix.pyx","src/scKTLD/cpd_utils/utils/_utils/convert_path_matrix_c.c",],),], language_level="3"),
    description="Identification of TAD-like domains on single-cell Hi-C data by graph embedding and changepoint detection",
    long_description="Identification of TAD-like domains on single-cell Hi-C data by graph embedding and changepoint detection"
)
