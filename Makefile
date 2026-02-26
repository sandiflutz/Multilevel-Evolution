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
	MEANBFRACxT \
	NUMHEVENTSxT \
	DENSb1xT \
	SAVE_CONFIG \
# -------------------------------------------------------
#  Compiler and flags
# -------------------------------------------------------
CC     = gcc
CFLAGS = -O3 -fsanitize=address -Wall -I includes/
LDFLAGS = -lm -fsanitize=address

# -------------------------------------------------------
#  Source files and object files
# -------------------------------------------------------
CORE = core/init.c core/evo.c core/measures.c core/tools.c core/randgen_ufrgs.c globals.c

SRCS = main_euler.c euler/bac_euler.c $(CORE)

OBJS = $(SRCS:.c=.o)

# -------------------------------------------------------
#  Final executable name
# -------------------------------------------------------
EXEC = exec_$(MACRO).out

# -------------------------------------------------------
#  Default rule
# -------------------------------------------------------
all: $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(LDFLAGS)
	@echo "Built executable: $(EXEC)"
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
	@echo "Example:"
	@echo "	make MACRO=AVERINVxT EXEC=invXt.out"
	@echo ""
	@echo "Other options:"
	@echo "  - clean objects: make cleanobj"
	@echo "  - clean executables: make cleanexec"
	@echo "  - clean all objects and executables: make cleanall"
	@echo ""
# -------------------------------------------------------
#  Phony targets
# -------------------------------------------------------
.PHONY: all clean
