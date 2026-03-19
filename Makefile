# -------------------------------------------------------
#  TERMINAL TEXT STYLE                                    
# -------------------------------------------------------
BOLD := \033[1m
RESET := \033[0m
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
	CORRxT \
	GENTIME \
	DIFBACOMPxT \
	DENSb1xT \
	SAVE_CONFIG \
	AVINVxRH \

TESTING_MACROS = \
	DEBUG_SELECT_EV \
# -------------------------------------------------------
#  Compiler and flags
# -------------------------------------------------------
CC     = gcc
CFLAGS = -O3 -fsanitize=address -Wall -I includes/
LDFLAGS = -lm -fsanitize=address

# -------------------------------------------------------
#  Source files and object files
# -------------------------------------------------------
CORE = core/init.c core/evo.c core/measures.c core/tools.c core/randgen_ufrgs.c core/bac.c globals.c

SRCS = main.c $(CORE)

OBJS = $(SRCS:.c=.o)

# -------------------------------------------------------
#  Final executable name
# -------------------------------------------------------
EXEC = exec_$(MACRO).out
# -------------------------------------------------------
#  Default rule
# -------------------------------------------------------
all: 
	@echo "To see options and syntax:\n make help\n"	
# -------------------------------------------------------
#  Main program
# -------------------------------------------------------
main: $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(LDFLAGS)
	@echo "Built executable: $(EXEC)"
# -------------------------------------------------------
#  DEBUG: for testing routines
# -------------------------------------------------------
SRCS_DB = debug.c $(CORE)

OBJS_DB = $(SRCS_DB:.c=.o)

debug: $(OBJS_DB)
	$(CC) $(OBJS_DB) -o $(EXEC) $(LDFLAGS)
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
	@echo "\nAvailable measurement MACRO options:"
	@$(foreach m,$(MEASURE_MACROS),echo "  - $(m)";)
	@echo ""
	@echo "Example:"
	@echo "	make main MACRO=AVERINVxT EXEC=invXt.out\n"
	@echo "Available testing MACRO options:"
	@$(foreach m,$(TESTING_MACROS),echo "  - $(m)";)
	@echo ""
	@echo "Example:"
	@echo "	make debug MACRO=DEBUG_SELECT_EV EXEC=db_selectev.out\n"
	@echo "Other options:"
	@echo "  - clean objects: make cleanobj"
	@echo "  - clean executables: make cleanexec"
	@echo "  - clean all objects and executables: make cleanall"
	@echo "\nImportant:"
	@echo "  - Choose executable names with $(BOLD).out$(RESET) as extention to make sure the options $(BOLD)cleanobj$(RESET) and $(BOLD)cleanall$(RESET) work properly"
	@echo "  - Remember to erase old objects before compiling the system.\n"
# -------------------------------------------------------
#  Phony targets
# -------------------------------------------------------
.PHONY: all clean
