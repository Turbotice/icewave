from setuptools import setup,find_packages

setup(
    name='icewave_codes',
    version='1.0',
    description='Measureing waves in sea ice',
#      url='needs a URL',
    author='Baptiste Auvity & Stéphane Perrard',
    author_email='stephane.perrard@espci.fr',
    license='GNU',
    packages=find_packages(),
    zip_safe=False,
#      package_data={'tangle': ['cl_src/*.cl']})
)
