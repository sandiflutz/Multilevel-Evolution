# -------------------------------------------------------
#  MACRO default name                                    
# -------------------------------------------------------
MACRO = TMEAS
# -------------------------------------------------------
#  List if available macros (for choosing what to measure)                                    
# -------------------------------------------------------
MEASURE_MACROS = \
	INV_DIST \
	AVERINVxT \
	DENSb1xT \
	SAVE_CONFIG \
# -------------------------------------------------------
#  Compiler and flags
# -------------------------------------------------------
CC     = gcc
CFLAGS = -O3 -fsanitize=address -Wall -I includes/
LDFLAGS_GSL = -lgsl -lgslcblas -lm -fsanitize=address
LDFLAGS_EULER = -lm -fsanitize=address

# -------------------------------------------------------
#  Source files and object files
# -------------------------------------------------------
CORE = core/init.c core/evo.c core/measures.c core/tools.c core/randgen_ufrgs.c globals.c

EULER_SRCS = main_euler.c euler/bac_euler.c $(CORE)
GSL_SRCS = main_gsl.c gsl/bac_gsl.c $(CORE)
DEBUG_SRCS = debug.c $(CORE)

EULER_OBJS = $(EULER_SRCS:.c=.o)
GSL_OBJS = $(GSL_SRCS:.c=.o)

DEBUG_OBJS = $(DEBUG_SRCS:.c=.o)

# -------------------------------------------------------
#  Final executable name
# -------------------------------------------------------
EXEC = exec_$(MACRO).out

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
	$(CC) $(EULER_OBJS) -o $(EXEC) $(LDFLAGS_EULER)
	@echo "Built EULER executable: $(EXEC)"
# -------------------------------------------------------
#  GSL:
# -------------------------------------------------------
gsl:	$(GSL_OBJS)
	$(CC) $(GSL_OBJS) -o $(EXEC) $(LDFLAGS_GSL)
	@echo "Built GSL executable: $(EXEC)"
# -------------------------------------------------------
#  DEBUG:
# -------------------------------------------------------
debug:	$(DEBUG_OBJS)
	$(CC) $(DEBUG_OBJS) -o $(EXEC) $(LDFLAGS_GSL)
	@echo "Built DEBUG executable: $(EXEC)"
# -------------------------------------------------------
#  Compilation step for each .c file
# -------------------------------------------------------
%.o: %.c
	$(CC) $(CFLAGS) -D$(MACRO) -c $< -o $@
# -------------------------------------------------------
#  Clean object files
# -------------------------------------------------------
cleanobj:
	find . -type f -name '*.o' -exec rm -f {} +
# -------------------------------------------------------
#  Clean executables
# -------------------------------------------------------
cleanexec:
	find . -type f -name '*.out' -exec rm -f {} +
# -------------------------------------------------------
#  Clean object and executables files
# -------------------------------------------------------
cleanall:
	find . -type f \( -name '*.o' -o -name '*.out' \) -exec rm -f {} +
# -------------------------------------------------------
#  Show the Available Macros on the .c files, used to
#  choose what to measure 
# -------------------------------------------------------
help:
	@echo ""
	@echo "Available measurement MACRO options:"
	@$(foreach m,$(MEASURE_MACROS),echo "  - $(m)";)
	@echo ""
	@echo "Other options:"
	@echo "clean just objects: make cleanobj"
	@echo "clean just executables: make cleanexec"
	@echo "clean all objects and executables: make cleanall"
	@echo ""
	@echo "Example (for the gsl version): make gsl MACRO=AVERINVxT EXEC=invXt_gsl.out"
	@echo "Example (for the euler version): make euler MACRO=AVERINVxT EXEC=invXt_eu.out"
	@echo ""
# -------------------------------------------------------
#  Phony targets
# -------------------------------------------------------
.PHONY: all clean
