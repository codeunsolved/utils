#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : base
# AUTHOR  : codeunsolved@gmail.com
# CREATED : August 10 2016
# VERSION : v0.0.1a4
# UPDATE  : [v0.0.1a1] February 13 2017
# 1. add `color_term()`, `execute_cmd()` and log related function/class;
# 2. add :FileHandlerFormatter: to remove ANSI color format when :SetupLogger: log to `FileHandler`;
# [ToDo] Not real time output with `subprocess.Popen()` in `execute_cmd()`;
# UPDATE  : [v0.0.1a2] May 18 2018
# 1. add `print()` action to `color_term()`;
# 2. optimize :StreamToLogger: to support seperate line;
# 3. optimize :SetupLogger: by adding `log()`, etc;
# UPDATE  : [v0.0.1a3] May 22 2018
# 1. seperate `colour()` from `color_term()`;
# 2. optimize :SetupLogger: by adding more specific log level functions;
# UPDATE  : [v0.0.1a4] November 13 2018
# 1. complete log levels in `colour` and make it consistent with :SetupLogger:;

import os
import re
import sys
import shlex
import logging
import subprocess


def colour(string, color='blue', bold=True):
    colors = {
        'RESET': '\033[0m',      # IFO - INFO; DBG - DEBUG

        'grey': '\033[0;30m',
        'red': '\033[0;31m',     # ERR - ERROR; CRT - CRITICAL
        'green': '\033[0;32m',
        'yellow': '\033[0;33m',  # WRN - WARNING
        'blue': '\033[0;34m',
        'magenta': '\033[0;35m',
        'cyan': '\033[0;36m',
        'white': '\033[0;37m',

        'END': '\033[0m',
    }
    # Log levels
    colors['IFO'] = colors['INFO'] = colors['RESET']
    colors['DBG'] = colors['DEBUG'] = colors['RESET']
    colors['WRN'] = colors['WARNING'] = colors['yellow']
    colors['ERR'] = colors['ERROR'] = colors['red']
    colors['CRT'] = colors['CRITICAL'] = colors['red']

    color_start = colors[color].replace('[0', '[1') if bold else colors[color]
    color_string = color_start + str(string) + colors['END']

    return color_string


def color_term(string, color='blue', bold=True,
               on_stream=True, end='\n', exit_code=None):
    color_string = colour(string, color=color, bold=bold)

    if on_stream:
        print(color_string, end=end)

    if exit_code is not None:
        exit(exit_code)

    return color_string


def execute_cmd(cmd):
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    for line in iter(p.stdout.readline, b''):
        sys.stdout.write(line.decode())

    error = p.stderr.read().decode()
    if error:
        raise Exception(error)


class FileHandlerFormatter(logging.Formatter):

    def format(self, record):
        msg = super(FileHandlerFormatter, self).format(record)
        return re.sub('\\033\[[\d;]+m', '', msg)


class SetupLogger(object):
    def __init__(self, name, path=None, level=logging.DEBUG, on_file=True, on_stream=True, log_mode='a',
                 format_fh='%(asctime)s | %(filename)s - line:%(lineno)-4d | %(levelname)s | %(message)s',
                 format_sh='[%(levelname)s] %(message)s',
                 format_date='[%b-%d-%Y] %H:%M:%S',
                 verbose_level=0):
        self.name = name
        self.path = self.handle_path(path)

        self.on_file = on_file
        self.on_stream = on_stream
        self.log_mode = log_mode

        self.format_fh = format_fh
        self.format_sh = format_sh
        self.format_date = format_date

        self.verbose_level = verbose_level

        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        if self.on_file:
            self.add_filehandler()
        if self.on_stream:
            self.add_streamhandler()

    def handle_path(self, path):
        if path is None:
            current_dir = os.path.abspath(os.path.dirname(__file__))
            return os.path.join(current_dir, 'log')
        else:
            return os.path.abspath(path)

    def add_filehandler(self):
        # Handle log directory
        dir_name = os.path.dirname(self.path)
        if not os.path.exists(dir_name):
            color_term("[SetupLogger] Log directory doesn't exist ", 'WRN', end='')
            color_term('CREATE ', 'grey', end='')
            try:
                os.makedirs(dir_name)
            except Exception as e:
                raise Exception(e)
            else:
                color_term('OK!', 'green')

        file_handler = logging.FileHandler(self.path, mode=self.log_mode)
        file_handler.setFormatter(FileHandlerFormatter(self.format_fh, self.format_date))
        self.logger.addHandler(file_handler)

    def add_streamhandler(self):
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter(self.format_sh))
        self.logger.addHandler(stream_handler)

    def remove_handler(self, class_=None, name=None):
        i = 0
        for handler in self.logger.handlers:
            if isinstance(handler, class_):
                self.logger.removeHandler(handler)
                print("{} removed".format(handler))
                i += 1
        if i == 0:
            print("No '{}' found!".format(name))

    def remove_filehandler(self):
        self.remove_handler(class_=logging.FileHandler, name='FileHandler')

    def remove_streamhandler(self):
        self.remove_handler(class_=logging.StreamHandler, name='StreamHandler')

    def log(self, msg, level='INFO', verbose=0, exit_code=None):
        level = str(level).upper()
        if level not in ['CRT', 'CRITICAL',
                         'ERR', 'ERROR',
                         'WRN', 'WARNING',
                         'IFO', 'INFO',
                         'DBG', 'DEBUG']:
            color_term("[SetupLogger] Unrecognized log level: {}, set to 'INFO'".format(level), 'ERR')
            level = 'INFO'

        if self.verbose_level >= verbose:
            if level in ['CRT', 'CRITICAL']:
                self.logger.critical(colour(msg, 'ERR'))
            elif level in ['ERR', 'ERROR']:
                self.logger.error(colour(msg, 'ERR'))
            elif level in ['WRN', 'WARNING']:
                self.logger.warning(colour(msg, 'WRN'))
            elif level in ['IFO', 'INFO']:
                self.logger.info(msg)
            elif level in ['DBG', 'DEBUG']:
                self.logger.debug(msg)

        if exit_code is not None:
            exit(exit_code)

    def critical(self, msg, verbose=0, exit_code=None):
        self.log(msg, level='CRT', verbose=verbose, exit_code=exit_code)

    def error(self, msg, verbose=0, exit_code=None):
        self.log(msg, level='ERR', verbose=verbose, exit_code=exit_code)

    def warning(self, msg, verbose=0, exit_code=None):
        self.log(msg, level='WRN', verbose=verbose, exit_code=exit_code)

    def info(self, msg, verbose=0, exit_code=None):
        self.log(msg, level='IFO', verbose=verbose, exit_code=exit_code)

    def debug(self, msg, verbose=0, exit_code=None):
        self.log(msg, level='DBG', verbose=verbose, exit_code=exit_code)


class StreamToLogger(object):
    """Fake file-like stream object that redirects writes to a logger instance.

    :Refer:
    1.  [Redirect stdout and stderr to a logger in Python]
        (https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/)
    2.  [How to redirect stdout and stderr to logger in Python]
        (https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python)
    """

    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, message):
        if message != '\n':
            self.linebuf += message
        else:
            self.logger.log(self.log_level, self.linebuf)
            self.linebuf = ''

    def flush(self):
        pass
