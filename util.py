import sys
import os
import ConfigParser
from os.path import expanduser
import logging


class UTIL:
	name2path = {}
	PATH = os.path.join(expanduser("~"), "UTIL.ini")
	config = ConfigParser.ConfigParser()
	#config.optionxform = str
	config.read(PATH)
	for section in config.sections():
		options = config.options(section)
		for option in options:
			try:
				logging.debug("[UTIL] %s = %s" %(option,config.get(section, option)))
				name2path[option] = config.get(section, option)
				if name2path[option.lower()] == -1:
					continue
			except:
				name2path[option] = None

	@staticmethod
	def get(name):
		if not os.path.exists(UTIL.name2path[name]):
			print UTIL.name2path[name], "doesn't exist!"
		return UTIL.name2path[name]

	def __init__(self):
		pass

		
	def __del__(self):
		pass
		
