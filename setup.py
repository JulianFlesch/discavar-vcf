import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="discavar-vcf",
    version="0.1",
    author="Julian Flesch",
    author_email="",
    description="Variant Cohort Analysis from VCFs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JulianFlesch/discavar-vcf",
    packages=["discavar"],
    install_requires=[
        "numpy",
        "pandas",
        "scikit-allel",
        "scipy",
        "cython",
        "cyvcf2"
        ],
    scripts=[
        "bin/discavar"],
    classifiers=[
        "Topics :: Bioinformatics ::  NGS Analysis",
        ""
        "Environment :: Console",
        "Programming Language :: Python :: 3",
        "License :: OSI APPROVED :: MIT License",
        "Operating System :: LINUX",
    ],
)
