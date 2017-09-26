# garbageIn

Any program that takes GFF files as input is totally not going to work. You
might think it is working, but as soon as new GFF file from a new organism
is fed into your program, well, end of game.

There are plenty of programs that can read a GFF file, the files are only
TAB-delimited tables after all. Should be easy right?

The fault of every other GFF reader I have ever seen is that they either
underestimate the wiles of GFFs or pass them off as the user's problem.

`garbageIn` is my attempt to build a single program that can deal reasonably
with any GFF input. It must properly handle evertying that is legal in the GFF
spec (an unmaintained rag that everyone reads and no one follows). It should
also handle the popular deviants (e.g. AUGUSTUS) who so shamelessly butcher the
format.  
