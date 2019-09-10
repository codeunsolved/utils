#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : test
# AUTHOR  : codeunsolved@gmail.com
# CREATED : September 7 2019
# VERSION : v0.0.1

import os
import re
import sys
import random
import argparse
import unittest

from test_ac import TestAnnCoordinates

__VERSION__ = 'v0.0.1'


def parse_config(path_config):
    if path_config and os.path.isfile(path_config):
        try:
            config_dir = os.path.dirname(path_config)
            config_name = os.path.basename(path_config)
            config_name = re.sub(r'\.py$', '', config_name)
            sys.path.append(config_dir)
            config = __import__(config_name)
        except Exception as e:
            print("Can't import config: {}".format(e))
        else:
            sys.path.pop()
            return config
    else:
        print("Invalid config file: {}".format(
            path_config))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="test",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description="test {}\n".format(__VERSION__))

    parser.add_argument('--config', type=str, default=None,
                        help="Specify config file(.py)")
    parser.add_argument('--mysql_user', type=str, default=None,
                        help="Specify MYSQL username")
    parser.add_argument('--mysql_pass', type=str, default=None,
                        help="Specify MYSQL password")
    parser.add_argument('--mysql_host', type=str, default='127.0.0.1',
                        help="Specify MYSQL host(default: 127.0.0.1)")
    parser.add_argument('--mysql_port', type=int, default=3306,
                        help="Specify MYSQL port(default: 3306)")
    parser.add_argument('--mysql_database', type=str, default=None,
                        help="Specify MYSQL databse name")
    parser.add_argument('--random_seed', type=int, default=None,
                        help="Specify random seed")

    args = parser.parse_args()

    config = parse_config(args.config)

    if config:
        mysql_config = getattr(config, 'MYSQL', {})
        TestAnnCoordinates.MYSQL_USER = mysql_config.get('user')
        TestAnnCoordinates.MYSQL_PASS = mysql_config.get('password')
        TestAnnCoordinates.MYSQL_HOST = mysql_config.get('host') or '127.0.0.1'
        TestAnnCoordinates.MYSQL_PORT = mysql_config.get('port') or 3306
        TestAnnCoordinates.MYSQL_DATABASE = mysql_config.get('database')
    else:
        TestAnnCoordinates.MYSQL_USER = args.mysql_user
        TestAnnCoordinates.MYSQL_PASS = args.mysql_pass
        TestAnnCoordinates.MYSQL_HOST = args.mysql_host
        TestAnnCoordinates.MYSQL_PORT = args.mysql_port
        TestAnnCoordinates.MYSQL_DATABASE = args.mysql_database

    if args.random_seed is None:
        TestAnnCoordinates.RANDOM_SEED = random.randint(1, 1000)
    else:
        TestAnnCoordinates.RANDOM_SEED = args.random_seed

    unittest.main(argv=['test.py'], verbosity=2)
