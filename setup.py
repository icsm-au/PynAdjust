from setuptools import setup, find_packages

setup(name='pynadjust',
      version='0.0.1',
      description='DynAdjust Python Package',
      long_description='A toolkit for working with DynAdjust inputs and outputs'
                       ' in Python',
      url='https://github.com/icsm-au/PynAdjust',
      author='ICSM Australia',
      author_email='geodesy@ga.gov.au',
      license='Apache License 2.0',
      packages=find_packages(),
      install_requires=['geodepy', 'pyshp', 'numpy'])
