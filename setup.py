DESCRIPTION = "Plan observation with ZTF"
LONG_DESCRIPTION = """Plan observation with ZTF. Currently only neutrino alert parsing is implemented"""

DISTNAME = "ztf_plan_obs"
AUTHOR = "Simeon Reusch"
MAINTAINER = "Simeon Reusch"
MAINTAINER_EMAIL = "simeon.reusch@desy.de"
URL = "https://github.com/simeonreusch/ztf_plan_obs/"
LICENSE = "BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/simeonreusch/ztf_plan_obs/archive/v0.24.tar.gz"
VERSION = "v0.24"

try:
    from setuptools import setup, find_packages

    _has_setuptools = True
except ImportError:
    from distutils.core import setup

    _has_setuptools = False


if __name__ == "__main__":

    install_requires = [
        "astropy==4.0.0",
        "numpy",
        "astroplan",
        "pandas",
        "matplotlib",
        "bs4",
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
        packages=packages,
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
