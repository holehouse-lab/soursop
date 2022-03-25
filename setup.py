"""
soursop
Analysis package for all-atom simulations of proteins, with a specific focus on intrinsically disordered proteins.
"""
import setuptools
from setuptools import setup
import versioneer

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:]),


setup(
    # Self-descriptive entries which should always be present
    name='soursop',
    author='Alex Holehouse',
    author_email='alex.holehouse@wustl.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    #scripts=['scripts/<something>'],
    license='LGPLv3',
    url='https://github.com/holehouse-lab/soursop',

    # Which Python importable modules should be included when your package is installed
    packages=['soursop', "soursop.tests"],

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    # > line below modified from default

    # line here means we include the files defined in MANIFEST.in - the package_data
    # line below is not sufficient to do this
    include_package_data=True,
    package_data={'soursop': ['data/'], 'soursop': ['data/test_data/']},

    # dependencies soursop requires
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.5.0",
        "cython",
        "mdtraj>=1.9.5",
        "pandas>=0.23.0",
        "threadpoolctl>=2.2.0",
    ],


    # Additional entries you may want simply uncomment the lines you want and fill in the data
    #author_email='alex.holehouse@wustl.edu',      # Author email
    #url='http://www.my_package.com',              # Website

    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
