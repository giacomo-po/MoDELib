# Add all header symlinks as dependencies to this target
header_symlinks::

MODEL_DIR := $(CURDIR)/../../../..
depend_dirs := $(MODEL_DIR)/model

# header files
include_dirs	:= $(shell find $(depend_dirs) -type d)
include_files	:= $(shell find $(depend_dirs) -regex "[^\#~]*\.h")
include_files	+= $(shell find $(depend_dirs) -regex "[^\#~]*\.hpp")

# all_header_directory
all_header_dir := $(MODEL_DIR)/header_symlinks

# header file links
link_names := $(foreach i, $(include_files), $(all_header_dir)/$(notdir $(i)))

#
# header symlinks
#
define all_header_dir_rule
$(1):
	@echo Rebuilding symlinks in $$@
	@$$(shell mkdir -p $$@)
endef

# Create a rule for one symlink for one header file
# Args
# 1: the header file
# 2: the symlink to create
define symlink_rule
$(2): $(1)
	@ln -sf $$< $$@
endef

# Create a rule for a symlink for each header file
# Args:
# 1: all_header_dir
# 2: list of header files
define symlink_rules
$(foreach i, $(2), $(eval $(call symlink_rule, $(i), $(1)/$(notdir $(i)))))
endef

$(eval $(call all_header_dir_rule, $(all_header_dir)))
$(call symlink_rules, $(all_header_dir), $(include_files))

header_symlinks:: $(all_header_dir) $(link_names)
