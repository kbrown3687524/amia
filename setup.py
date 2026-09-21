from setuptools import setup

setup(
    name='amia',
    version='1.1.3',    
    description='Automated Mutation Introduction and Analysis',
    author='Keaghan Brown',
    author_email='3687524@myuwc.ac.za',
    packages=['amia'],
    install_requires=[
        'click>=8.1',
        'PyYAML>=6.0',
        'biopython>=1.79',
        'MDAnalysis[analysis]>=2.4.0',
        'MDAnalysisTests',
        'pandas>=1.5.3',
        'numpy>=1.24.2',
        'scipy>=1.10.1',
        'matplotlib>=3.7.1',
        'networkx>=3.0',
    ],
    entry_points={
        "console_scripts": [
            "amia=amia.AMIA:run_pipeline",  # <-- matches your @click.command run_pipeline()
        ],
    },
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9'
    ],
)
