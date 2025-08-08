from setuptools import setup, find_packages

setup(
    name="editable_app",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "regex",
        "biopython",
        "requests",
        "beautifulsoup4"
    ],
)
