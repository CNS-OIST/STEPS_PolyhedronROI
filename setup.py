import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="polyhedronROI",
    version="2.0.0",
    author="Weiliang Chen",
    author_email="w.chen@oist.jp",
    description="Create Regions of Interest (ROI)s in STEPS Tetmesh object using Polyhedral surface boundaries.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CNS-OIST/STEPS_PolyhedronROI",
    packages=setuptools.find_packages(),
    install_requires=["trimesh", "numpy", "rtree", "typing"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
