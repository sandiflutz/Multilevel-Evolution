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
	AVERINVxT \
	INV_DIST \
	NUMHEVENTSxT \
	CORRxT \
	CLUSTERSxT \
	CLUSTERS_DISTxT \
	BESTCLUSTER_TIMES \
	RVNxTxCORRBAC \
	DIFBACOMPxT \
	SAVE_CONFIG \
	AVINVxRH \
	AVINVxGH \
	AVINVxMB \
	AVINVxMH \
	AVINVxCOST \
	RHxMHxAVINV \
	COSTxMBxAVINV \
# -------------------------------------------------------
#  Compiler and flags
# -------------------------------------------------------
CC     = gcc
CFLAGS_RELEASE = -O3 -Wall -I includes/
LDFLAGS_RELEASE = -lm

CFLAGS_DEBUG = -O3 -fsanitize=address -Wall -I includes/
LDFLAGS_DEBUG = -lm -fsanitize=address

# -------------------------------------------------------
#  Default MODE
# -------------------------------------------------------
Mode ?= release
# -------------------------------------------------------
#  If a MODE is specified
# -------------------------------------------------------
ifeq ($(MODE),debug)
	CFLAGS = $(CFLAGS_DEBUG)
	LDFLAGS = $(LDFLAGS_DEBUG)
else
	CFLAGS = $(CFLAGS_RELEASE)
	LDFLAGS = $(LDFLAGS_RELEASE)
endif
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
all: $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(LDFLAGS)
	@echo "Built executable: $(EXEC)"
# -------------------------------------------------------
#  Compilation step for each .c file
# -------------------------------------------------------
%.o: %.c
	$(CC) $(CFLAGS) -D$(MACRO) -c $< -o $@
# -------------------------------------------------------
#  Debug with valgrind
# -------------------------------------------------------
valgrind: all
	valgrind --leak-check=full --track-origins=yes ./$(EXEC)
# -------------------------------------------------------
#  Clean object files
# -------------------------------------------------------
cleanobj:
	find . -type f -name '*.o' -exec rm -fv {} +
# -------------------------------------------------------
#  Clean executables
# -------------------------------------------------------
cleanexec:
	find . -type f -name '*.out' -exec rm -f {} +
# -------------------------------------------------------
#  Clean object and executables files
# -------------------------------------------------------
cleanall:
	find . -type f \( -name '*.o' -o -name '*.out' \) -exec rm -fv {} +
# -------------------------------------------------------
#  Show the Available Macros on the .c files, used to
#  choose what to measure 
# -------------------------------------------------------
help:
	@echo "\nAvailable measurement MACRO options:"
	@$(foreach m,$(MEASURE_MACROS),echo "  - $(m)";)
	@echo ""
	@echo "Example:"
	@echo "	make MACRO=AVERINVxT EXEC=invXt.out\n"
	@echo "For debugging while executing:"
	@echo "	make MODE=debug MACRO=AVERINVxT EXEC=invXt.out\n"
	@echo "For debugging and executing with valgring:"
	@echo "	make valgrind MODE=DEBUG MACRO=AVERINVxT EXEC=invXt.out\n"
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
