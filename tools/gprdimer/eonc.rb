#!/usr/bin/env ruby

require "tty-command" # For the Command
require "tty-prompt" # For the method choice
require "tty-file" # For the replacement

prompt = TTY::Prompt.new
a = prompt.select("Choose the method", %w(GPRDimer Dimer))

# TODO: Super fragile, update
re = TTY::File.replace_in_file "config.ini", /(?<=min_mode_method = ).*$/, a.downcase
if re==false
  # Couldn't replace, exit
  exit
end

# Run eonclient
cmd = TTY::Command.new(printer: :pretty)
out, err = cmd.run("eonclient")
File.write("./%s_stdout" % [a.downcase], out)
File.write("./%s_stderr" % [a.downcase], err)

# Cleanup
# Move generated files around
File.rename "client.log", "./%s_client.log" % [a.downcase]
File.rename "mode.dat", "./%s_mode.dat" % [a.downcase]
File.rename "results.dat", "./%s_results.dat" % [a.downcase]
File.rename "saddle.con", "./%s_saddle.con" % [a.downcase]

# TODO: Fixup
File.delete("runAMS.sh") if File.exist?("runAMS.sh")
File.delete("updCoord.sh") if File.exist?("updCoord.sh")
File.delete("myrestart.in.sh") if File.exist?("myrestart.in.sh")

FileUtils.remove_dir("firstRun.results") if File.directory?("firstRun.results")
FileUtils.remove_dir("secondRun.results") if File.directory?("secondRun.results")
FileUtils.remove_dir("ams.results") if File.directory?("ams.results")
