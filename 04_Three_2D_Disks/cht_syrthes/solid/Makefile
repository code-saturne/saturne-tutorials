SHELL = /bin/sh

#*****************************************************
#***     Makefile pour la géneration 
#***     d'un executable pour SYRTHES
#*****************************************************

include $(SYRTHES_HOME)/share/syrthes/Makefile.in


OBJ=$(patsubst %.c, %.o, $(wildcard *.c))


exe : $(OBJ)
	ar xv $(LIB_SYRTHES) mainsyrthes.o
	$(CC) -o syrthes   $(COPTIM) $(CFLAG) -D_FILE_OFFSET_BITS=64 \
	$(INC_SYRTHES_LIBMAT) $(INC_SYRTHES) $(INC_SYRTHES_STARDIS)  $(INC) $(LDFLAG) *.o \
	$(LIB_SYRTHES) -lm  $(LIB_STARDIS) $(LIB) $(LIB_PLE) 

	@rm *.o

%.o : %.c
	$(CC)  -c  $(COPTIM) $(CFLAG) -D_FILE_OFFSET_BITS=64 $(INC_SYRTHES_LIBMAT) $(INC_SYRTHES)  $(INC_SYRTHES_STARDIS) $(INC_STARDIS)  $(INC) $? -o $@


