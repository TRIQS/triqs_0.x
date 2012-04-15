#!/usr/bin/env python
import sys,os,re

dic =  [ 
( r'\s*boost::fusion' , 'bf') ,
(r',\s*bf::void_' , ''),
(r',\s*void_' , ''),
(r',\s*void' , ''),
(r'triqs::arrays::' , 'tqa::'),
(r'\s*boost::proto' , 'proto'),
(r'\s*proto::tagns_::tag' , 'p_tag'),
#r"proto::exprns_::basic_expr<(.*)>" , r"\1"), 
(r"proto::exprns_::basic_expr" , r""), 
(r"proto::argsns_::list2",""),
(r"proto::argsns_::term",""),
(r",\s*mpl_::na",""),
]


for line in sys.stdin.readlines():
  # find the [0x adress line...]
  m = re.search('\[(0x[0-9A-Fa-f]*)\]', line)
  if m : 
      executable = line.split()[0]
      addr = m.group(1) 
      command = "addr2line -C -f -e %s %s " %(executable,addr)
      cin,cout,cerr = os.popen3(command)  # execute in a unix pipe.
      res = cout.readlines()
      #print res
      #s =  line + "---> %s \n---> %s\n"%(res[0].strip(),res[1].strip())
      s =  line + "---> %s\n"%(res[1].strip())
      for i,r in dic : 
         s = re.sub(i,r,s)
      print s 
  else : print line
