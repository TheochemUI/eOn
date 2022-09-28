#!/usr/bin/env -S ruby -w
# frozen_string_literal: true

# Wraps eonclient in a friendly helper. Generates timestamped output folders,
# and a prompt driven customization routine.
#
# Usage examples:
#
# To run a standard GPRD calculation:
#
#   eonc.rb # Then select the right options
#
# Importantly, this is meant to be run with a "fat" configuration file, along with the three required files:
# [1] post.con
# [2] displacement.con
# [3] direction.dat
#
# The "fat" configuration file essentially must be called config.ini
# Furthermore, the following must be present **as written**!!
# - use_prune = false
# - start_prune_at
# - nprune_vals
# - prune_threshold
#
# This is because regular expressions are used to ensure the correct configuration is set up.
#
# TODO: Add a config reader, replace regexps
# TODO: Add templates to remove the requirement for the config
# TODO: Add a classic CLI with option parsing as well
# TODO: Generate a README within the directory
# TODO: Add more logging
# TODO: Support potentials

require 'tty-command' # For the Command
require 'tty-prompt' # For the minmodemethod choice
require 'tty-file' # For the replacement
require 'solid_assert' # Need assert

def replace_value_conf(input_key, new_value)
  # TODO: Super fragile, update
  # This expects that the named but discarded capture group is present, and in
  # particular will break if spaces are not present around the = symbol
  re = TTY::File.replace_in_file 'config.ini', /(?<=#{input_key} = ).*$/, new_value.to_s
  assert re == true # Can't fail
end

prompt = TTY::Prompt.new

# Populate prompt
minmodemethod = prompt.select('Choose the minmodemethod', %w[GPRDimer Dimer])
if minmodemethod.downcase == 'gprdimer'
  prune_status = prompt.yes?('With pruning?', default: false)
  if prune_status == true
    prune_defaults = prompt.yes?(
      "Use the defaults?\n This is the removal of 4 elements greater than 0.5 starting from size 10", default: true
    )
    prune_sched = if prune_defaults == true
                    { 'start_prune_at' => 10, 'nprune_vals' => 4, 'prune_threshold' => 0.5 }
                  else # User needs to provide the values
                    prompt.collect do
                      key('start_prune_at').ask('Size to start pruning at?', convert: :int)
                      key('nprune_vals').ask('How many values are to be dropped per thinning round?', convert: :int)
                      key('prune_threshold').ask('Initial highest allowed force value?', convert: :float)
                    end
                  end
  end
end
run_eonclient = prompt.yes?('Run eonclient?')

# Prepare output string
out_root = if minmodemethod.downcase == 'gprdimer'
             "#{minmodemethod.downcase}_prune_#{prune_status}"
           else
             minmodemethod.downcase.to_s
           end

# We would like to ideally run every one of these in tandem, all at once without
# having to worry about overwrites, this means timestamping a folder and moving
# to it. We use the seconds since epoch for better granularity.
# TODO: This way we can even make a pretty README
# HACK: Might not be the most elegant approach. Consider when multiple of these
# are started concurrently
puts 'Moving one level down'
mytime = Time.now.to_i
orig_dir = Dir.pwd
lowered_dir = "#{orig_dir}/#{out_root}_#{mytime}"
Dir.mkdir(lowered_dir)
FileUtils.cp %w[direction.dat displacement.con pos.con config.ini], lowered_dir.to_s

# HACK: I do not like this, the actual shift is jarring
Dir.chdir(lowered_dir)
puts Dir.pwd

# Unconditional
replace_value_conf('min_mode_method', minmodemethod.downcase)
replace_value_conf('use_prune', prune_status)
prune_sched.each { |key, value| replace_value_conf(key, value) } if prune_status == true

# Run
if run_eonclient == false
  puts "Writing to #{out_root}_stdout and #{out_root}_stderr"
  cmd = TTY::Command.new(printer: :pretty, dry_run: true)
  print TTY::File.diff_files("#{orig_dir}/config.ini", "#{lowered_dir}/config.ini")
  cmd.run('eonclient')
  # Revert directories
  puts 'Deleting directories'
  Dir.chdir(orig_dir)
  FileUtils.remove_dir(lowered_dir)
else
  cmd = TTY::Command.new(printer: :pretty)
  out, err = cmd.run('eonclient')
  File.write("./#{out_root}_stdout", out)
  File.write("./#{out_root}_stderr", err)

  # Cleanup
  # Move generated files around
  File.rename 'client.log', "./#{out_root}_client.log"
  File.rename 'mode.dat', "./#{out_root}_mode.dat"
  File.rename 'results.dat', "./#{out_root}_results.dat"
  File.rename 'saddle.con', "./#{out_root}_saddle.con"

  # TODO: Fixup
  File.delete('runAMS.sh') if File.exist?('runAMS.sh')
  File.delete('updCoord.sh') if File.exist?('updCoord.sh')
  File.delete('myrestart.in.sh') if File.exist?('myrestart.in.sh')

  FileUtils.remove_dir('firstRun.results') if File.directory?('firstRun.results')
  FileUtils.remove_dir('secondRun.results') if File.directory?('secondRun.results')
  FileUtils.remove_dir('ams.results') if File.directory?('ams.results')
end
