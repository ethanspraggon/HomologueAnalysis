from glob import glob
import os
from setuptools import setup, find_packages
import shutil
import subprocess
import sys
# first party
###################################################################################################
try:
    from Utilitor.Buildor.clean_wheels import clean_wheels
    UTILITOR = True
except ImportError:
    print("No Utilitor all ok will just ignore this package")
    UTILITOR = False

# setup date based version number
###################################################################################################
from datetime import datetime

name = "GenomicHomologues"
now = datetime.now()
year = now.strftime("%Y")
month = now.strftime("%m")
month = month[1:] if month[0:1] == "0" else month
version = ".".join([year, month])

setup(
    name=name,
    version=version,
    description="Sequence Queries and manipulation",
    author="Ethan Spraggon",
    author_email="easpraggon@gmail.com",
    packages=find_packages(),
    package_dir={name: name},
    package_data={name: ["Metadata/Uniprot/*"]},
    entry_points={
        "console_scripts": [

            "find_homologues = HomologueAnalysis.ensembl:homologues"
        ]
    },
    # url='http://www.python.org',
    #all the packages and dependecies
    install_requires=["biopython>=1.6", "scipy>=0.18"]
    # extras_require=["scipy>=0.18]  doesnt work yet
)
# stuff to copy to wheels
###################################################################################################
if "bdist_wheel" in sys.argv[1:] and UTILITOR:
    to_path = clean_wheels(name, version, remove=True)
    # install
    ######################################################################
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", str(to_path)])
    # clean
    ######################################################################
    subprocess.check_call([sys.executable, "setup.py", "clean", "--all"])
    files = glob("*egg-info")
    for file in files:
        shutil.rmtree(file)


