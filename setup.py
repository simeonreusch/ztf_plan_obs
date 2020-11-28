DESCRIPTION = "Plan observations with the Zwicky Transient Facility"
LONG_DESCRIPTION = """Plan observations with the Zwicky Transient Facility. Automated parsing of GCN is only implemented for IceCube at the moment."""

DISTNAME = "ztf_plan_obs"
AUTHOR = "Simeon Reusch"
MAINTAINER = "Simeon Reusch"
MAINTAINER_EMAIL = "simeon.reusch@desy.de"
URL = "https://github.com/simeonreusch/ztf_plan_obs/"
LICENSE = "BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/simeonreusch/ztf_plan_obs/archive/v0.32.tar.gz"
VERSION = "0.32"

try:
    from setuptools import setup, find_packages

except ImportError:
    from distutils.core import setup

    raise Exception("Please install python3 setuptools")


if __name__ == "__main__":

    install_requires = [
        "astropy",
        "astroquery",
        "numpy",
        "astroplan>=0.7",
        "pandas",
        "matplotlib",
        "bs4",
        "flask",
        "ztfquery",
        "requests",
        "lxml",
        "html5lib",
    ]

    setup(
        name=DISTNAME,
        author=AUTHOR,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        packages=find_packages(),
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3.8",
            "License :: OSI Approved :: BSD License",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
    )
