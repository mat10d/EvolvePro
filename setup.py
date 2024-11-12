from setuptools import setup, find_packages

setup(
    name="evolvepro",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scikit-learn',
        'scikit-learn-extra',
        'xgboost',
        'matplotlib',
        'seaborn',
        'biopython',
        'scipy'
    ]
)