from setuptools import setup, find_packages

setup(
    name="evolvepro",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'biopython',
        # add any other dependencies here
    ],
)