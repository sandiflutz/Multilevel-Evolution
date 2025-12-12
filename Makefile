# -------------------------------------------------------
#  MACRO default name                                    
# -------------------------------------------------------
MACRO = NOTHING

# -------------------------------------------------------
#  Compiler and flags
# -------------------------------------------------------
CC     = gcc
CFLAGS = -O3 -D$(MACRO) -fsanitize=address -Wall -I includes/
LDFLAGS_GSL = -lgsl -lgslcblas -lm -fsanitize=address
LDFLAGS = -lm -fsanitize=address

# -------------------------------------------------------
#  Source files and object files
# -------------------------------------------------------
CORE = core/init.c core/evo.c core/measures.c core/tools.c core/randgen_ufrgs.c globals.c

EULER_SRCS = main_euler.c euler/bac_euler.c $(CORE)
GSL_SRCS = main_gsl.c gsl/bac_gsl.c $(CORE)

EULER_OBJS = $(EULER_SRCS:.c=.o)
GSL_OBJS = $(GSL_SRCS:.c=.o)

# -------------------------------------------------------
#  Final executable name
# -------------------------------------------------------
EXEC = exec.out

# -------------------------------------------------------
#  Default rule:buid the gsl version
# -------------------------------------------------------
all: 
	@echo "Choose one of: make euler, make gsl"
	@exit 1

# -------------------------------------------------------
#  Euler
# -------------------------------------------------------
euler:	$(EULER_OBJS)
	$(CC) $(EULER_OBJS) -o $(EXEC) $(LDFLAGS)
	@echo "Built EULER executable: $(EXEC)"
# -------------------------------------------------------
#  GSL:
# -------------------------------------------------------
gsl:	$(GSL_OBJS)
	$(CC) $(GSL_OBJS) -o $(EXEC) $(LDFLAGS_GSL)
	@echo "Built GSL executable: $(EXEC)"
# -------------------------------------------------------
#  Compilation step for each .c file
# -------------------------------------------------------
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# -------------------------------------------------------
#  Clean euler object files
# -------------------------------------------------------
cleaneuler:
	rm -f *.o core/*.o euler/*.o 
# -------------------------------------------------------
#  Clean gsl object files
# -------------------------------------------------------
cleangsl:
	rm -f *.o core/*.o gsl/*.o  
# -------------------------------------------------------
#  Phony targets
# -------------------------------------------------------
.PHONY: all clean
