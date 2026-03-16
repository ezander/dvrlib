SRCS = $(shell find src test demo -name '*.h' -o -name '*.cc')

.PHONY: all demo test coverage coverage-html check clean


BUILD     = build
BUILD_COV = build-cov
COV_OUT   = coverage-html

all: $(BUILD)/CMakeCache.txt
	cmake --build $(BUILD) -- -s

demo: all
	$(BUILD)/dvrlib_demo

test: all
	@$(BUILD)/dvrlib_test

coverage: $(BUILD_COV)/CMakeCache.txt
	cmake --build $(BUILD_COV) -- -s
	$(BUILD_COV)/dvrlib_test
	@cd $(BUILD_COV)/CMakeFiles/dvrlib.dir/src && \
	    gcov *.gcno 2>/dev/null \
	        | awk '/File.*\/src\/[^/]+\.cc/{f=1; print} f && /Lines executed/{print; f=0}'; \
	    rm -f *.gcov 2>/dev/null || true

coverage-html: $(BUILD_COV)/CMakeCache.txt
	@which lcov > /dev/null 2>&1 || { echo "lcov not found. Install with: sudo apt install lcov"; exit 1; }
	cmake --build $(BUILD_COV) -- -s
	$(BUILD_COV)/dvrlib_test
	lcov --capture --directory $(BUILD_COV) --output-file $(BUILD_COV)/coverage.info
	lcov --remove $(BUILD_COV)/coverage.info '*/build-cov/_deps/*' '/usr/*' --output-file $(BUILD_COV)/coverage.info
	genhtml $(BUILD_COV)/coverage.info --output-directory $(COV_OUT)
	@echo "Report: $(COV_OUT)/index.html"

$(BUILD)/CMakeCache.txt:
	cmake -B $(BUILD) -S . -DCMAKE_BUILD_TYPE=Debug

$(BUILD_COV)/CMakeCache.txt:
	cmake -B $(BUILD_COV) -S . -DCMAKE_BUILD_TYPE=Debug \
	    -DCMAKE_CXX_FLAGS="--coverage" \
	    -DCMAKE_EXE_LINKER_FLAGS="--coverage"

check: $(BUILD)/CMakeCache.txt
	# Commented out, clang-tidy reports too much (sometimes dangerous) nonsense
	#	@echo "--- clang-tidy ---"
	#	@find src test demo -name '*.cc' | xargs clang-tidy -p $(BUILD)
	@echo "--- clang-format ---"
	@clang-format --dry-run --Werror $(SRCS) && echo "format OK"

format_source:
	@echo "--- clang-format ---"
	@clang-format $(SRCS) -i --verbose

clean:
	rm -rf $(BUILD) $(BUILD_COV) $(COV_OUT)
