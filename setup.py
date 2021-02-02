from setuptools import setup

setup(name='bayespairing',
      version='2.0',
      description='RNA 3D module prediction',
      url='http://csb.cs.mcgill.ca/BayesPairing2/',
      author='Roman Sarrazin-Gendron',
      author_email='roman.sarrazingendron@mail.mcgill.ca',
      license='McGill',
      packages=["bayespairing","bayespairing.src"],
      install_requires=[
          "networkx", "matplotlib","seaborn","weblogo","ghostscript","biopython","wrapt","anytree", "treedecomp"
      ],
      zip_safe=False)
