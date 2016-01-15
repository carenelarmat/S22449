#!/bin/tcsh -f

#New compilation flags
set NEWF=" -axAVX,SSE4.2 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -diag-disable 6477 -implicitnone -gen-interfaces -warn all -O3 -check nobounds -DFORCE_VECTORIZATION"
set NEWC="-axAVX,SSE4.2  "'$(CPPFLAGS)'
set files=`find . -name 'Makefile'`
foreach file ($files)
  echo $file
  set index=2
  eval "sed  '/^FLAGS_CHECK/s/[^=]*/$NEWF/'"'$index '$file > toto
  eval "sed  '/^CFLAGS/s/[^=]*/$NEWC/'"'$index 'toto > $file
end
