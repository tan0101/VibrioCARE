import sys

from setuptools import setup

PACKAGE_NAME = 'VibrioCARE'
DESCRIPTION = 'A set of tools in Python developed for the article "Signature genomic traits of the core and accessory genome of Vibrio Cholerae O1 drive lineage transmission and disease severity"'
with open('README.md', encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()
AUTHOR = 'Alexandre Maciel-Guerra, Kubra Babaarslan, Michelle Baker, Aura Rahman, Maqsud Muhammad Hossain, Abdus Sadique, Jahidul Alam, Salim Uzzaman, Mohammad Ferdous Rahman Sarker, Nasrin Sultana, Ashraful Islam Khan, Yasmin Ara Begum, Mokibul Hassan Afrad, Nicola Senin, Zakir Hossain Habib, Tahmina Shirin, Firdausi Qadri, and Tania Dottorini'
AUTHOR_EMAIL = 'tania.dottorini@nottingham.ac.uk'
URL = 'https://github.com/tan0101/VibrioCARE'
MINIMUM_PYTHON_VERSION = 3, 9  # Minimum of Python 3.9
VERSION = '0.0.1'

def check_python_version():
    """Exit when the Python version is too low."""
    if sys.version_info < MINIMUM_PYTHON_VERSION:
        sys.exit("Python {}.{}+ is required.".format(*MINIMUM_PYTHON_VERSION))


check_python_version()

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    install_requires=[
        'python==3.9.15',
        'numpy==1.21.5',
        'pandas==1.4.4',
        'sklearn==1.0.2',
        'scipy==1.9.3',
        'networkx==2.8.4',
        'matplotlib==3.6.2',
        'imblearn==0.10.1',
        'cobra==0.24.0'
    ],
    license='Apache License 3.0',
)
