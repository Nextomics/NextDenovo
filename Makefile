CC=gcc

TOP_DIR=$(PWD)
BIN_DIR=$(TOP_DIR)/bin
SUB_DIR=minimap2 util lib

TARGET_MINIMAP=minimap2-nd
TARGER_UTIL=seq_dump ovl_sort seq_stat ovl_cvt seq_bit bam_sort nextgraph

all:make_bin_dir
	$(foreach N, $(SUB_DIR), make -C $(N);)
	cp $(TOP_DIR)/minimap2/$(TARGET_MINIMAP) $(TOP_DIR)/bin/
	$(foreach N, $(TARGER_UTIL), cp $(TOP_DIR)/util/$(N) $(TOP_DIR)/bin/;)

make_bin_dir:
	@if test ! -d $(BIN_DIR) ; \
	then \
		mkdir $(BIN_DIR) ; \
	fi

clean:
	rm -rf $(BIN_DIR)
	$(foreach N, $(SUB_DIR), make -C $(N) clean;)
