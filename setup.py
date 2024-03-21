from setuptools import find_packages, setup

setup(
    name='scat',
    packages=find_packages(),
    package_data={
        'scat': ['bin/browseDefinitions.xlsx'],
    },
    version='0.0',
    description='GUI to analyze spectral image cubes.',
    author='Michael S. Phillips',
    # license='CC-BY-NC-SA-4.0',
)
