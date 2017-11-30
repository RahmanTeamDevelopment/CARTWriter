from setuptools import setup

setup(
    name = 'CARTWriter',
    version = '0.2.0',
    description = 'A simple tool that outputs CARTs data in various formats',
    url = 'https://github.com/RahmanTeamDevelopment/CARTWriter',
    author = 'RahmanTeam',
    author_email = 'rahmanlab@icr.ac.uk',
    license = 'MIT',
    packages=["cartwriter"],
    scripts=["bin/CARTWriter.py", "bin/cartwriter"],
    zip_safe=False
)
