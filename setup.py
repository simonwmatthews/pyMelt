import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyMelt",
    version="3.dev0", # REMEMBER TO UPDATE __INIT__.PY
    author="Simon Matthews, Kevin Wong, Matthew Gleeson",
    author_email="simonwmatthews@gmail.com",
    description=("A python library for calculating the melting behaviour of Earth's mantle."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/simonwmatthews/pyMelt",
    packages=setuptools.find_packages(),
    package_data={'pyMelt':['phaseDiagrams/*', 'phaseDiagrams/build/*']},
    install_requires=[
            'pandas==1.4.4',
            'numpy==1.21.5',
            'matplotlib==3.5.2',
            'scipy==1.9.1',
            ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
