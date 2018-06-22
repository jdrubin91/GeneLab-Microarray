from setuptools import setup


setup(name='GeneLab-Microarray',version='0.1', 
    description='Standardized processing pipeline for microarray data on GeneLab.', 
    url='https://github.com/jdrubin91/GeneLab-Microarray', 
    author='Jonathan Rubin & Daniel Mattox', 
    author_email='jonathan.d.rubin@nasa.gov', 
    license='GPL-3.0',
    packages=['GeneLab-Microarray'], 
    install_requires=['scipy', 'mpld3','matplotlib'],
    zip_safe=False,
    scripts=['bin/GeneLab-Microarray'])
