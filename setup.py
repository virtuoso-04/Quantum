"""
Setup script for Quantum Drug Discovery package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="quantum-drug-discovery",
    version="1.0.0",
    author="Quantum Chemistry Research Team",
    author_email="your.email@example.com",
    description="Quantum computing pipeline for drug discovery and molecular binding calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/virtuoso-04/Quantum",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "pandas>=1.3.0",
        "tqdm>=4.62.0",
    ],
    extras_require={
        "quantum": [
            "pyscf>=2.0.0",
            "qiskit>=0.45.0",
            "qiskit-nature>=0.7.0",
            "qiskit-algorithms>=0.2.0",
        ],
        "viz": [
            "seaborn>=0.12.0",
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "ipywidgets>=8.0.0",
            "notebook>=6.4.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
        ],
        "docs": [
            "sphinx>=4.5.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "quantum-drug-discovery=examples.full_pipeline_demo:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yaml", "*.json", "*.txt"],
    },
    zip_safe=False,
)
