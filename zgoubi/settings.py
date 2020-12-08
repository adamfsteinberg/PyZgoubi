#!/usr/bin/env python
"Read settings from ~/.pyzgoubi"
from __future__ import division, print_function

import configparser
import os
import tempfile

home = os.path.expanduser('~')
config_dir = os.path.join(home, ".pyzgoubi")

if not os.path.exists(config_dir):
	os.mkdir(config_dir)


config_path = os.path.join(config_dir, "settings.ini")

config = configparser.RawConfigParser()
config.add_section('pyzgoubi')
#default values
config.set('pyzgoubi', 'tmp_dir', tempfile.gettempdir())
config.set('pyzgoubi', 'extra_defs_files', '')
config.set('pyzgoubi', 'zgoubi_path', "zgoubi")
config.set('pyzgoubi', 'log_level', "warn")
config.set('pyzgoubi', 'max_label_size', 20)



if os.path.exists(config_path):
	config.read(config_path)
else:
	config.write(open(config_path, 'w'))


zgoubi_settings = {}
zgoubi_settings['tmp_dir'] = os.path.expanduser(config.get('pyzgoubi', 'tmp_dir'))
#zgoubi_settings['pyzgoubi_dir'] = config.get('pyzgoubi', 'pyzgoubi_dir')
#zgoubi_settings['defs_dir'] = config.get('pyzgoubi', 'defs_dir')

zgoubi_settings['zgoubi_path'] = os.path.expanduser(config.get('pyzgoubi', 'zgoubi_path'))

zgoubi_settings['extra_defs_files'] = [os.path.expanduser(x) for x in config.get('pyzgoubi', 'extra_defs_files').split(',')]
if zgoubi_settings['extra_defs_files'] == ['']:
	zgoubi_settings['extra_defs_files'] = []

zgoubi_settings['log_level'] = config.get('pyzgoubi','log_level').upper()
if not zgoubi_settings['log_level'] in ['WARN', 'ERROR', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL']:
	print("invalid setting for log_level '%s' in settings.ini"%zgoubi_settings['log_level'])
	print("should be one of 'WARN', 'ERROR', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL'")
	exit(2)

zgoubi_settings['max_label_size'] = int(config.get('pyzgoubi', 'max_label_size'))

# create example defs file
example_defs_path = os.path.join(config_dir, "user_elements.defs")

if not os.path.exists(example_defs_path):
	example_defs = open(example_defs_path, 'w')
	example_defs.write("# define any elements you need here, the will be compiled to ~/.pyzgoubi/defs_cache\n")

