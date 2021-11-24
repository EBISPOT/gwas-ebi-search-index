from setuptools import setup

setup(
    name='gwas-ebi-search-index',
    description='GWAS search index document generator',
    version='0.0.1',
    packages=['.'],
    include_package_data=True,
    license='Apache License, Version 2.0',
    entry_points={
        "console_scripts": ['create-ebi-search-index-data = create_ebi_search_index_data:main']
    },
    url='https://github.com/EBISPOT/gwas-ebi-search-index'
)