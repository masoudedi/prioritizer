import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="genome",
    version="0.3.1",
    author="Masoud-Edizadeh",
    author_email="edizadeh.masoud@yahoo.com",
    description="A comprehensive pipeline for analyzing WES",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject1",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues1",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)