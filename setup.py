from setuptools import setup, find_packages

setup(
    name='pyGeSyn',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'matplotlib>=3.5',
        'numpy>=1.20',
    ],
    entry_points={
        'console_scripts': [
            'pyGeSyn=pyGeSyn.cli:main',
        ],
    },
    python_requires='>=3.12',
    author='pyGeSyn',
    description='Synteny visualization for genomic regions',
    classifiers=[
        'Programming Language :: Python :: 3.12',
    ],
)
