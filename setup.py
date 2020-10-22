from setuptools import setup

with open('README.rst') as f:
	readme = f.read()

with open('LICENSE') as f:
	license = f.read()

setup(name='mutation_gatherer',
	version='1.0',
	description='Gather variant information from NCI GDC API, ClinVar API, and specific CSV/TSV and make a variant lollipop plot',
	url='https://github.com/samaco-lab/mutation_gatherer',
	long_description=readme,
	author='Rocio Dominguez Vidana, PhD',
	author_email='rocio.vidana@alumni.bcm.edu',
	license='MIT',
	packages=['mutation_gatherer']
)
