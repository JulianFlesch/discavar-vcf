import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="discavar-vcf",
    version="0.1",
    author="Julian Flesch",
    author_email="",
    keywords="bioinformatics ngs-analysis vcf cohort disease causing",
    description="Variant Cohort Analysis from VCFs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JulianFlesch/discavar-vcf",
    project_urls={
        "Documentation": "https://packaging.python.org/tutorials/distributing-packages/",
        "Source": "https://github.com/JulianFlesch/discavar-vcf"
        },
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
