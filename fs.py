#To use this module, you need to include the following in your code:
#import sys
#if not "~/python_lib" in sys.path: sys.path.append("~/python_lib")
#from fs import FS
import os
#sed cheat sheet
#http://sed.sourceforge.net/sed1line.txt


class FS:
	def __init__(self):
		self.filename2handle = {}

	def file_exist(self, filename):
		if not os.path.exists(filename):
			print "[ERROR] File NOT exists:", filename
			return False
		else:
			return True

	def get_lines_between(self, infile, outfile, line_from, line_to): # line numbers inclusive
	#print section of file based on line numbers (lines 8-12, inclusive)
	#sed -n '8,12p' filename	# method 1
	#sed '8,12!d'			# method 2
		pass


	def read2list(self, filename):
	#use this function to read file with small to medium size
	#Example:
	#fs = FS()
	#mylist = fs.read(r"E:\python_lib\FS_example.txt"):
	#print "The file is read into a list", mylist
		content = []
		fh = open(filename)
		#for i, line in enumerate(fh.readlines()):
		#	print i, line
		for line in fh:
			#print line
			line = line.rstrip('\n')
			line = line.rstrip('\r')
			content.append(line)
		return content

	
	def read2list_skip_header(self, filename):
	#use this function to read file with small to medium size
	#Example:
	#fs = FS()
	#mylist = fs.read(r"E:\python_lib\FS_example.txt"):
	#print "The file is read into a list", mylist
		content = []
		fh = open(filename)
		#for i, line in enumerate(fh.readlines()):
		#	print i, line
		fh.next()
		for line in fh:
			#print line
			line = line.rstrip('\n')
			line = line.rstrip('\r')
			content.append(line)
		return content



	def read(self, filename):
	#use this function to read file line by line
	#Example:
	#fs = FS()
	#for line in fs.read(r"E:\python_lib\FS_example.txt"):
	#	print line.rstrip()
		fh = open(filename)
		if filename in self.filename2handle:
			fh = self.filename2handle[filename]
		else:
			fh = open(filename)
			self.filename2handle[filename] = fh
		return fh


	def read_header(self, filename):
		fh = open(filename)
		if filename in self.filename2handle:
			fh = self.filename2handle[filename]
		else:
			fh = open(filename)
			self.filename2handle[filename] = fh
		return fh.next()
	
	def read_skip_header(self, filename):
		if filename in self.filename2handle:
			fh = self.filename2handle[filename]
		else:
			fh = open(filename)
			fh.next()
			self.filename2handle[filename] = fh
		return fh
	#To read large file line by line
	#fh = open("infile.txt")
	#header = fh.next()
	#for line in fh:
	#	sample_list.append(line[:10].rstrip())	
	#fh.close()

	def write(self, filename, str_buffer):
		if filename in self.filename2handle:
			#print "existing file"
			fh = self.filename2handle[filename]
		else:
			#print "new file"
			fh = open(filename, "wt")
			self.filename2handle[filename] = fh
		fh.write(str_buffer)

	def close(self):
		for filename, fh in self.filename2handle.items():
			fh.close()
			del self.filename2handle[filename]
			print "[CLOSE]", filename
			
	def close_file(self, filename):
		if filename in self.filename2handle:
			self.filename2handle[filename].close()
			del self.filename2handle[filename]
			print "[CLOSE]", filename
		else:
			print "No file handle found for", filename




def main():

	fs = FS()
	for line in fs.read(r"E:\python_lib\FS_example.txt"):
		print line.rstrip()
	fs.close()
	#print "\n*****************************\n"
	#print "header:", fs.read_header(r"E:\python_lib\Excel.py")
	#for line in fs.readSkipHeader(r"E:\python_lib\Excel.py"):
	#	print line
	#	break
	#fs.close()
	#fs.write(r"E:\python_lib\FS_example.txt", "test")
	#fs.close()


if __name__ == '__main__':
	main()


