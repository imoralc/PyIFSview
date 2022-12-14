import pathlib
from setuptools import find_packages, setup

HERE = pathlib.Path(__file__).parent

VERSION = '1.0' 
PACKAGE_NAME = 'PyIFSview' 
AUTHOR = 'Ignacio del Moral Castro' 
AUTHOR_EMAIL = 'ignaciodelmoralcastro@gmail.com' 
URL = 'https://github.com/imoralc/PyIFSview' 

LICENSE = 'MIT' 
DESCRIPTION = 'Library to interactively visualize integral field spectroscopy (IFS) data' 
LONG_DESCRIPTION = (HERE / "README.md").read_text(encoding='utf-8')
LONG_DESC_TYPE = "text/markdown"
INCLUDE_PACKAGE_DATA=True,


INSTALL_REQUIRES = ['numpy', 'astropy', 'matplotlib']

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESC_TYPE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    install_requires=INSTALL_REQUIRES,
    license=LICENSE,
    packages=find_packages(),
    include_package_data=True
)