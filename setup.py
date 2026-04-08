from setuptools import setup, find_packages

# Read requirements from requirements.txt
def read_requirements():
    with open('requirements.txt', 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name="immunex",
    version="0.1.0",
    description="Immunex - Comprehensive Immunological MD Analysis Toolkit",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    author="Research Team",
    author_email="research@example.com",
    url="https://github.com/your-username/Immunex",
    packages=find_packages(),
    install_requires=read_requirements(),
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="immunology, molecular dynamics, TCR-pMHC, MD analysis, PBC processing, bioinformatics",
    entry_points={
        'console_scripts': [
            # Main unified CLI
            'imn=immunex.cli.main:main',

            # Legacy commands (backward compatibility)
            'immunex=immunex.cli.main:main',
            'immunex-batch=immunex.core.batch_process:main',
            'immunex-slurm=immunex.cluster.slurm_generator:main',
            'immunex-batch-pdb=immunex.cli.batch_pdb:main',
        ],
    },
    extras_require={
        'dev': [
            'pytest>=7.0',
            'pytest-cov>=4.0',
            'black>=23.0',
            'flake8>=6.0',
        ],
        'docs': [
            'sphinx>=5.0',
            'sphinx-rtd-theme>=1.0',
        ],
    },
)
