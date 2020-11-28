#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.

from setuptools import setup
import glob

setup_args = {
    'name': "utility_packages",
    'description': "miscellaneous utility methods/scripts/etc",
    'license': "BSD",
    'author': "David DeBoer",
    'author_email': "david.r.deboer@gmail.edu",
    'version': '0.1',
    'scripts': glob.glob('scripts/*'),
    'packages': ['my_utils', 'my_ephem'],
    'package_data': {"my_utils": ["data/*"]}
}

if __name__ == '__main__':
    setup(**setup_args)
