import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "dnaapler",
            "VERSION",
        )
    ) as f:
        return f.readline().strip()
    

def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="dnaapler",
    packages=find_packages(),
    url="",
    python_requires=">=3.7",
    description="Re-orients  bacterial chromosome seqeunces to begin with dnaA gene",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="George Bouras",
    author_email="george.bouras@adelaide.edu.au",
    data_files=get_data_files(),
    py_modules=["dnaapler"],
    install_requires=[
        "pyyaml>=6.0",
        "Click>=8.1.3",
        "pandas>=1.4.2",
        "loguru>=0.5.3", 
        "biopython==1.79",
        "pyrodigal>=2.0.0",
        "pytest>=6.2.5",
        "pytest-cov>=3.0.0"
    ],
    entry_points={
        "console_scripts": [
            "dnaapler=dnaapler.__main__:main" 
        ]
    },
    include_package_data=True,
)
