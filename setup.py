import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="peepdemo", 
    version="1.0.0",
    author="Seri Kitada",
    author_email="kitada.seri.77w@st.kyoto-u.ac.jp",
    description="Designs paired pegRNAs for deletions and evaluates efficiencies",
    url='git@github.com:serikitada/PEEP_demo/tree/main',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
	'numpy'==1.22.3,
	'pandas'==1.4.1,
	'tensorflow'==2.8.0,
	'biopython'==1.74,
	'viennarna'
    ]
)
