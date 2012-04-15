#!/usr/bin/env python
import sys,re

dic =  [ 
( r'\s*boost::fusion' , 'bf') ,
(r',\s*bf::void_' , ''),
(r',\s*void_' , ''),
(r'\s*boost::proto' , 'proto'),
(r'\s*proto::tagns_::tag' , 'p_tag'),
#r"proto::exprns_::basic_expr<(.*)>" , r"\1"), 
(r"proto::exprns_::basic_expr" , r""), 
(r"proto::argsns_::list2",""),
(r"proto::argsns_::term",""),
(r",\s*mpl_::na",""),
]


#s = open(sys.argv[1]).read()
s = sys.stdin.read()

for i,r in dic : 
  s = re.sub(i,r,s)
  #print i,r ,s

print s

