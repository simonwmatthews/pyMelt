import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyMelt",
    version="1.92",
    author="Simon Matthews, Kevin Wong, Matthew Gleeson",
    author_email="simonwmatthews@gmail.com",
    description=("A python library for calculating the melting behaviour of Earth's mantle."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/simonwmatthews/pyMelt",
    packages=setuptools.find_packages(),
    install_requires=[
            'pandas',
            'numpy',
            'matplotlib',
            'scipy',
            ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
