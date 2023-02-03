"""
pore_c_py setup script.
"""
import os
import pkg_resources
import re
from setuptools import setup, find_packages
import sys

__pkg_name__ = 'pore-c-py'
__author__ = 'epi2melabs'
__description__ = 'Analyse Pore-C data.'

__pkg_import__ = __pkg_name__.replace("-", "_")

# Use readme as long description and say its github-flavour markdown
from os import path
this_directory = path.abspath(path.dirname(__file__))
kwargs = {'encoding':'utf-8'} if sys.version_info.major == 3 else {}
with open(path.join(this_directory, 'README.md'), **kwargs) as f:
    __long_description__ = f.read()
__long_description_content_type__ = 'text/markdown'

__path__ = os.path.dirname(__file__)
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_import__))

# Get the version number from __init__.py
verstrline = open(os.path.join(__pkg_import__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_import__))

dir_path = os.path.dirname(__file__)
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    install_requires = [
        str(requirement) for requirement in
        pkg_resources.parse_requirements(fh)]

data_files = []
extra_requires = {}
extensions = []


setup(
    name=__pkg_name__,
    version=__version__,
    url='https://github.com/epi2me-labs/{}'.format(__pkg_name__),
    author=__author__,
    author_email='{}@nanoporetech.com'.format(__author__),
    description=__description__,
    long_description=__long_description__,
    long_description_content_type=__long_description_content_type__,
    zip_safe=False,
    install_requires=install_requires,
    packages=find_packages(exclude=("tests",)),
    entry_points={
        "console_scripts": [
            f"{__pkg_name__} = {__pkg_import__}.main:run_main"]},
)
