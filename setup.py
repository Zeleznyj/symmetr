import setuptools
from symmetr.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="symmetr",
    #version=__version__,
    version='0.7.1',
    author="Jakub Zelezny",
    author_email="jakub.zelezny@gmail.com",
    description="Package for determining symmetry properties of crystals.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/zeleznyj/linear-response-symmetry",
    packages=setuptools.find_packages(),
    package_data={'symmetr': ['findsym/*']},
    scripts = ['exec/symmetr'],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Operating System :: POSIX :: Linux",
    ],
    install_requires=[
        'future',
        'six',
        'numpy',
        'sympy',
        'prettytable']
)
