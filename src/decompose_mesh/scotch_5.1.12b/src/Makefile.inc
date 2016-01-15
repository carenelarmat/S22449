EXE		=
LIB		= .a
OBJ		= .o

MAKE		= make
AR		= ar
ARFLAGS		= -ruv
CAT		= cat
CCS		= /usr/projects/hpcsoft/toss2.1/common/intel/2011.11.339/composer_xe_2011.11.339/bin/intel64/icc
CCP		= mpicc
CCD		= /usr/projects/hpcsoft/toss2.1/common/intel/2011.11.339/composer_xe_2011.11.339/bin/intel64/icc
CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -DSCOTCH_PTHREAD -Drestrict=__restrict -pthread
CLIBFLAGS	=
LDFLAGS		= -lz -lm -lrt 
CP		= cp
LEX		= flex -Pscotchyy -olex.yy.c
LN		= ln
MKDIR		= mkdir
MV		= mv
RANLIB		= ranlib
YACC		= bison -y -pscotchyy -b y
