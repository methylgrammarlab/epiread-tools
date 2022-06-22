import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='epiread-tools',
    version='0.0.1',
    author='Irene Unterman',
    author_email='irene.guberman@mail.huji.ac.il',
    description='Tools for handling epiread files',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/methylgrammarlab/epiread-tools',
    project_urls = {
        "Bug Tracker": "https://github.com/methylgrammarlab/epiread-tools/issues"
    },
    license='MIT',
    packages=['epiread-tools'],
    install_requires=['argparse', 'numpy', 'pandas', 'scipy', 'pysam'],
)