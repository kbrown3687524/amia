from setuptools import setup
setup(
    name='amia',
    version='0.1.0',    
    description='Automated Mutation Introduction and Analysis',
    author='Keaghan Brown',
    author_email='3687524@myuwc,.ac.za',
    packages=['amia'],
    install_requires=['biopython>=1.79',
                      'pymol>=2.3.0',
                      'MDAnalysis[analysis]>=2.4.0',
                      'MDAnalysisTests',
                      'pandas>=1.5.1',
		      'numpy>=1.24.2',
	              'scipy>=1.10.1',
		      'matplotlib>=3.7.1',
		      'networkx>=3.0',
		      'pandas>=1.5.3'
                      ],

    classifiers=[
        'Development Status :: 1 - Development',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: ',  
        'Operating System :: Linux',        
        'Programming Language :: Python :: 3',
	  'Programming Language :: Python :: 3.9'
    ],
)
