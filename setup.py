import setuptools
import re
import os
import codecs

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

PKG = "svtool"

version = find_version(PKG, "__version__.py")

setuptools.setup(
    name=PKG,
    version=version,
    author="Huanchang Chen, Quanhu Sheng",
    author_email="hua-chang.chen@vumc.org, quanhu.sheng.1@vumc.org",
    description="Struature Variant Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zeissmania/sv_tool",
    download_url="https://github.com/zeissmania/sv_tool/archive/v" + version + ".tar.gz",
    entry_points = {
        'console_scripts': ['svtool=svtool.__main__:main'],
    },
    packages=setuptools.find_packages(exclude=["tests", "tests.*", "scripts", "debug.py", ".project", ".pydevproject"]),
    package_data={'': ['svtool/*.r', 'svtool/*.rmd', 'svtool/data/*.bed', 'svtool/data/*.txt', 'svtool/data/*.count', 'svtool/data/*.zip']},
    install_requires=['argparse', 'pathlib', 'pysam', 'pytabix', 'biopython' ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    zip_safe=False
)

