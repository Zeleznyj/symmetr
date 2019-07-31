import setuptools

with open("README_pypi.md", "r") as fh:
    long_description = fh.read()

def get_version():
    with open("symmetr/version.py") as f:
        lines = f.readlines()
    for line in lines:
        if 'version' in line:
            version = line.split('=')[1].strip().lstrip('\'').rstrip('\'')
            return version

setuptools.setup(
    name="symmetr",
    version=get_version(),
    author="Jakub Zelezny",
    author_email="jakub.zelezny@gmail.com",
    description="Package for determining symmetry properties of crystals.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/zeleznyj/linear-response-symmetry",
    packages=setuptools.find_packages(),
    package_data={'symmetr': ['findsym/*']},
    data_files=[('',['README_pypi.md'])],
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
        'numpy; python_version>="3"',
        'numpy==1.16.4; python_version<"3"',
        'sympy',
        'prettytable']
)
