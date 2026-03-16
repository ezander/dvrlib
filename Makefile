.PHONY: all demo test coverage clean

BUILD     = build
BUILD_COV = build-cov
COV_OUT   = coverage-html

all: $(BUILD)/CMakeCache.txt
	cmake --build $(BUILD) -- -s

demo: all
	$(BUILD)/dvrlib_demo

test: all
	@$(BUILD)/dvrlib_main -d yes

coverage: $(BUILD_COV)/CMakeCache.txt
	cmake --build $(BUILD_COV) -- -s
	$(BUILD_COV)/dvrlib_main
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

clean:
	rm -rf $(BUILD) $(BUILD_COV) $(COV_OUT)
