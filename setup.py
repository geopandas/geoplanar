import os
import versioneer
from distutils.command.build_py import build_py
from setuptools import setup

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists("MANIFEST"):
    os.remove("MANIFEST")


with open("requirements.txt") as f:
    tests_require = f.readlines()
install_requires = [t.strip() for t in tests_require]

with open("README.md") as f:
    long_description = f.read()


setup(
    name="geoplanar",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass({"build_py": build_py}),
    description="Geographic planar enforcement of polygon geoseries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sjsrey/geoplanar",
    author="Serge Rey",
    author_email="sjsrey@gmail.com",
    license="3-Clause BSD",
    packages=["geoplanar"],
    package_data={"": ["requirements.txt"]},
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.10",
    install_requires=install_requires,
    zip_safe=False,
)
