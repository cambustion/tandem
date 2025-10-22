#!/usr/bin/env python3
"""
Setup script for Aerosol PSD Analysis Tool
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="aerosol-psd-analysis",
    version="1.0.0",
    author="Morteza Kiasadegh",
    author_email="kiasadeg@ualberta.ca",
    description="Multi-lognormal fitting tool for tandem particle size distribution analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MortezaKiasadegh/aerosol-psd-analysis",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "aerosol-psd-fit=aerosol_psd_fitting_tool:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

