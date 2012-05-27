import sys

for f in sys.argv[1:] : 
    status = 0
    for line in open(f) :
      if line [0:3] == "/**" : status = 1
      if line.rstrip() [-3:] == "**/" : status = 3
      if status ==2 : print line.rstrip('\n') 
      if status ==3 : 
          print "/*************** %s *******************/"% f.upper()
          status = 2
    print "\f"
