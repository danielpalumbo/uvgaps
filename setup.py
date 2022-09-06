from setuptools import setup

setup(name='uvgaps',
      version='1.0',
      description='Package for characterizing circular gaps in u-v coverage.',
      url='https://github.com/danielpalumbo/uvgaps',
      author='danielpalumbo',
      author_email='daniel.palumbo@cfa.harvard.edu',
      license='GPLv3',
      packages=['uvgaps'],
      install_requires=['numpy','scipy','time','matplotlib','ehtim'])
