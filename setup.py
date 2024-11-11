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
    package_data={'pyMelt':['phaseDiagrams/depma_holland2018/make_depma_holland2018.ipynb', 
                            'phaseDiagrams/depma_holland2018/table_depma_holland2018.csv',
                            'phaseDiagrams/kg1_holland2018/make_kg1_holland2018.ipynb', 
                            'phaseDiagrams/kg1_holland2018/table_kg1_holland2018.csv',
                            'phaseDiagrams/g2_holland2018/make_g2_holland2018.ipynb', 
                            'phaseDiagrams/g2_holland2018/table_g2_holland2018.csv',
                            'phaseDiagrams/klb1_holland2018/make_klb1_holland2018.ipynb', 
                            'phaseDiagrams/klb1_holland2018/table_klb1_holland2018.csv',]},
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
