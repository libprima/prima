This folder contains the intersection-format version of the Fortran source
files. The files in this folder are generated automatically by the interform.m
script, and they are NOT intended to be readable.

This project is coded in the free format, yet some platforms accept only
fixed-format Fortran code, for example, the MATLAB MEX on Windows. The code
in this folder can serve such a purpose.

In the intersection format, each continued line has an ampersand at
column 73, and each continuation line has an ampersand at column 6.
A Fortran file in such a format can be compiled both as a fixed-format
file and as a free-format file.
See http://fortranwiki.org/fortran/show/Continuation+lines for details.

Zaikun Zhang (www.zhangzk.net), 03-Aug-2020.