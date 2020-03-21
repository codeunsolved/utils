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
import logging
import argparse
import unittest

from test_ac import TestAnnCoordinates_GENCODE, TestAnnCoordinates_REFSEQ
from test_mysql_cnn import TestMysqlConnector

__VERSION__ = 'v0.0.1'


def parse_config(path_config):
    if path_config and os.path.isfile(path_config):
        try:
            config_dir = os.path.dirname(path_config)
            config_name = os.path.basename(path_config)
            config_name = re.sub(r'\.py$', '', config_name)
            sys.path.insert(0, config_dir)  # Make this priority
            config = __import__(config_name)
        except Exception as e:
            print("Can't import config: {}".format(e))
        else:
            sys.path.pop(0)
            return config
    else:
        print("Invalid config file: {}".format(
            path_config))


def set_config(args):
    def check_if_mysql_config_is_set():
        if args.mysql_user and args.mysql_pass and \
           args.mysql_host and args.mysql_port and args.database:
            return True
        else:
            return False

    config = parse_config(args.config)

    # Set MySQL attributes
    for cls in [TestAnnCoordinates_GENCODE, TestMysqlConnector]:
        if check_if_mysql_config_is_set():
            cls.MYSQL_USER = args.mysql_user
            cls.MYSQL_PASS = args.mysql_pass
            cls.MYSQL_HOST = args.mysql_host
            cls.MYSQL_PORT = args.mysql_port
            cls.MYSQL_DATABASE = args.mysql_database
        elif config:
            mysql_config = getattr(config, 'MYSQL', {})
            cls.MYSQL_USER = mysql_config.get('user')
            cls.MYSQL_PASS = mysql_config.get('password')
            cls.MYSQL_HOST = mysql_config.get('host') or '127.0.0.1'
            cls.MYSQL_PORT = mysql_config.get('port') or 3306
            cls.MYSQL_DATABASE = mysql_config.get('database')
    # Set AnnCoordinates RefSeq attributes
    if args.refseq_files:
        TestAnnCoordinates_REFSEQ.REFSEQ_FILES = args.refseq_files
    elif config:
        refseq_files = getattr(config, 'REFSEQ_FILES', None)
        TestAnnCoordinates_REFSEQ.REFSEQ_FILES = refseq_files


def set_random_seed(random_seed):
    for cls in [TestAnnCoordinates_GENCODE, TestAnnCoordinates_REFSEQ]:
        if random_seed is None:
            cls.RANDOM_SEED = random.randint(1, 1000)
        else:
            cls.RANDOM_SEED = random_seed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="test",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description="test {}\n".format(__VERSION__))

    parser.add_argument('-m', '--module', type=str, default='__main__',
                        help="Specify test module(default: all)")
    parser.add_argument('--config', type=str, default=None,
                        help="Specify config file(.py), "
                             "Its priority is lower than other parameters below")
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
    parser.add_argument('--refseq_files', type=str, nargs='+',
                        help="Specify path of RefSeq file(s)")
    parser.add_argument('--random_seed', type=int, default=None,
                        help="Specify random seed")

    args = parser.parse_args()

    logger = logging.basicConfig(level=logging.CRITICAL)  # Kepp silence

    set_config(args)

    set_random_seed(args.random_seed)

    unittest.main(module=args.module, argv=['test.py'], verbosity=2)
