#!/usr/bin/env ruby

require "tty-command" # For the Command
require "tty-prompt" # For the minmodemethod choice
require "tty-file" # For the replacement
require "solid_assert" # Need assert

def replace_value_conf(input_key, new_value)
  # TODO: Super fragile, update
  re = TTY::File.replace_in_file "config.ini", /(?<=#{input_key} = ).*$/, "#{new_value}"
  assert re == true # Can't fail
end

prompt = TTY::Prompt.new

# Populate prompt
minmodemethod = prompt.select("Choose the minmodemethod", %w(GPRDimer Dimer))
if minmodemethod.downcase == "gprdimer"
  prune_status = prompt.yes?("With pruning?", default: false)
  if prune_status == true
    prune_defaults = prompt.yes?("Use the defaults?\n This is the removal of 4 elements greater than 0.5 starting from size 10", default: true)
    if prune_defaults == true
      prune_sched = {'start_dropout_at' => 10, 'ndropout_vals' => 4, 'dropout_threshold' => 0.5}
    else # User needs to provide the values
      prune_sched = prompt.collect do
        key('start_dropout_at').ask("Initial dropout size at?", convert: :int)
        key('ndropout_vals').ask("How many values are to be dropped per thinning round?", convert: :int)
        key('dropout_threshold').ask("Initial highest allowed force value?", convert: :float)
    end
    end
  end
end
run_eonclient = prompt.yes?("Run eonclient?")

# Unconditional
replace_value_conf("min_mode_method", minmodemethod.downcase)
replace_value_conf("use_dropout", prune_status)
if prune_status == true
  prune_sched.each {|key, value| replace_value_conf(key, value) }
end

# Generate output string
if  minmodemethod.downcase == "gprdimer"
  out_root = "./#{minmodemethod.downcase}_prune_#{prune_status}"
else
  out_root = "./#{minmodemethod.downcase}"
end

# Run
if run_eonclient == false
  puts "#{out_root}_stdout"
  puts "#{out_root}_stderr"
  cmd = TTY::Command.new(printer: :pretty)
  TTY::File.tail_file("config.ini", 15) do |line|
    puts line
  end
else
  out, err = cmd.run("eonclient")
  File.write("#{out_root}_stdout", out)
  File.write("#{out_root}_stderr", err)

  # Cleanup
  # Move generated files around
  File.rename "client.log", "#{out_root}_client.log"
  File.rename "mode.dat", "#{out_root}_mode.dat"
  File.rename "results.dat", "#{out_root}_results.dat"
  File.rename "saddle.con", "#{out_root}_saddle.con"

  # TODO: Fixup
  File.delete("runAMS.sh") if File.exist?("runAMS.sh")
  File.delete("updCoord.sh") if File.exist?("updCoord.sh")
  File.delete("myrestart.in.sh") if File.exist?("myrestart.in.sh")

  FileUtils.remove_dir("firstRun.results") if File.directory?("firstRun.results")
  FileUtils.remove_dir("secondRun.results") if File.directory?("secondRun.results")
  FileUtils.remove_dir("ams.results") if File.directory?("ams.results")
end

