from setuptools import setup,find_packages

setup(
    name='phonefleet_codes',
    version='1.0',
    description='Control a smartphone fleet with Python',
#      url='needs a URL',
    author='Stéphane Perrard',
    author_email='stephane.perrard@espci.fr',
    license='GNU',
    packages=find_packages(),
    zip_safe=False,
#      package_data={'tangle': ['cl_src/*.cl']})
)
