import sys
import os
import logging
import ConfigParser
from os.path import expanduser


class REF:
	name2path = {}
	PATH = os.path.join(expanduser("~"), "REF.ini")
	config = ConfigParser.ConfigParser()
	#config.optionxform = str
	config.read(PATH)
	for section in config.sections():
		options = config.options(section)
		for option in options:
			try:
				logging.debug("[REF] %s = %s" %(option,config.get(section, option)))
				name2path[option.lower()] = config.get(section, option)
				if name2path[option] == -1:
					continue
			except:
				name2path[option] = None

	@staticmethod
	def get(name):
		if not os.path.exists(REF.name2path[name]):
			print REF.name2path[name], "doesn't exist!"
		return REF.name2path[name]

	def __init__(self):
		pass
		
			
		
	def __del__(self):
		pass
		
