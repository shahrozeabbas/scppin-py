"""Setup script for scppin package (backward compatibility)."""

from setuptools import setup, find_packages

# For modern Python packaging, most configuration is in pyproject.toml
# This file is kept for backward compatibility

setup(
    packages=find_packages(),
    include_package_data=True,
)

