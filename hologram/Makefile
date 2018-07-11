
### Programe name
NAME = hologram
BIN  = bin/$(NAME)
DATE = $(shell date +'%y.%m.%d-%H.%M.%S')

####### Compiler, tools and options
CC            = g++
DEFINES       =
CFLAGS        = -m64 -pipe -std=c++0x -O3 -Wall -W $(DEFINES)
INCPATH       = -I./src -I/usr/local/include
LINK          = g++
LFLAGS        = -L/usr/local/lib -Wl,-R /usr/local/lib '-Wl,-R $$ORIGIN'
LIBS          =

####### Directories
# Dossier
SRCD = src/
OBJD = src/

# Liste des fichiers
SRCS=$(wildcard $(SRCD)*.cpp $(SRCD)*/*.cpp $(SRCD)*/*/*.cpp)



OBJS = $(SRCS:$(SRCD)%.cpp=$(OBJD)%.o)
DEPS = $(SRCS:$(SRCD)%.cpp=$(OBJD)%.d)

########## Regles ##########
default:
	@mkdir -p bin
	@make -s c
	@make -s all

all:
	@echo "\033[44m\033[30m=========== Compilation =============\033[m\033[m"
	@echo "Will put bin to $(BIN)"
	@make $(BIN)
	@echo "\033[42m\033[30m=============== OK ==================\033[m\033[m"

# Executable
$(BIN): $(OBJS)
	@$(CC) -o $@ $^ $(LFLAGS) $(CFLAGS) $(LIBS)

# Objets
$(OBJD)%.o: $(SRCD)%.cpp
	@echo "\033[43m\033[30m -> $@\033[m\033[m"
	@$(CC) $(INCPATH) -o $@ -c $< $(CFLAGS)

# Dependances
$(OBJD)%.d: $(SRCD)%.cpp
	@echo "\033[43m\033[30m -> $@\033[m\033[m"
	@$(CC) $(CFLAGS) -MM -MD -o $@ $<

#-include $(DEPS)

clean:
	@echo "\033[41m\033[30m============= CLEANING ==============\033[m\033[m"
	@echo "Will remove :"
	@echo "- $(OBJS)"
	@echo "- $(DEPS)"
	@rm -rf $(OBJS) $(DEPS)

# Execute la commande "clean", puis supprime les executables et la doc.
veryclean: clean
	@echo "- $(BIN)/$(NAME)"
	@rm -f $(BIN)/$(NAME)

cleanQt: clean
	@echo "- *Debug"
	@echo "- *Release"
	@echo "- src/*.o"
	@echo "- src/*.d"
	@rm -rf *Debug
	@rm -rf *Release
	@rm -rf src/*.o src/*.d

verycleanQt: cleanQt
	@echo "- src/*.user"
	@echo "- $(BIN)/$(NAME)"
	@rm -f $(BIN)/$(NAME)
	@rm -f src/*.user

archive: verycleanQt
	@make
	@make cleanQt
	@echo "\033[45m\033[30m============= ARCHIVE ===============\033[m\033[m"
	@mkdir -p archive
	@zip -r archive/$(NAME).$(DATE).zip $(SRCD) $(BIN) src/Makefile src/$(NAME).pro Makefile --exclude src/$(NAME).pro.user

.PHONY: clean c doc archiver
