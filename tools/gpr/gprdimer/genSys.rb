#!/usr/bin/env -S ruby -w
# frozen_string_literal: true

# Script to facilitate generating multiple configuration runs from one "fat"
# config
#  TODO: Document
#  TODO: Add cli, optparse
#  TODO: Clean up into classes

require 'fileutils' # move things
# require "liquid" # nicer documentation
require 'tty-file' # for diffs
require 'solid_assert' # need assert
require 'inifile' # read a fat configuration

def condel(x)
  FileUtils.remove_dir(x) if Dir.exist?(x)
  File.delete(x) if File.exist?(x)
end

# HACK: This has the rather major side effect of moving into a folder
def lower_dir(out_root, orig_dir)
  puts 'Moving one level down'
  mytime = Time.now.to_i
  lowered_dir = "#{orig_dir}/#{out_root}_#{mytime}"
  Dir.mkdir(lowered_dir)
  # XXX: Does not generate configuration file
  FileUtils.cp %w[direction.dat displacement.con pos.con full.sample.ini], lowered_dir.to_s
  lowered_dir
end

# TODO: Bundle this into a class ==> LowerDir
def relative_diff(orig_dir, lowered_dir)
  # Must be run after the configuration is written out
  # Assumes that the relative root only has the full.sample.ini file
  reldiff = TTY::File.diff_files("#{orig_dir}/full.sample.ini", "#{lowered_dir}/config.ini")
  Dir.chdir(lowered_dir)
  File.write("#{lowered_dir}/rel.diff", reldiff)
end

# Main stuff

# Full config, needed only once
fconf = IniFile.load("#{Dir.pwd}/full.sample.ini")

# First lowering to make an "experiment"
ldir = lower_dir('exp1', Dir.pwd)
Dir.chdir(ldir)

# Make a "bash" executor
bashfsc = File.new('exec.sh', 'w+')

%w[gprdimer dimer].each do |method|
  mdir = lower_dir(method, Dir.pwd) # LOWERING!
  rconf = IniFile.new
  rconf = fconf
  rconf['Saddle Search']['min_mode_method'] = method
  rconf.filename = "#{mdir}/config.ini"
  rconf.write
  relative_diff(ldir, mdir) # Now is in the lowered directory
  # Not needed
  # bashfsc.write("\n cd #{mdir}/\n eonclient 2> error.txt 1> output.txt &")
  # End
  # Now do the potentials
  { 'Hybrid B3LYP' => 'ADF', 'PM3' => 'MOPAC' }.each do |model, engine|
    pdir = lower_dir(model, mdir) # LOWERING!
    if model == 'PM3'
      rconf['Potential']['potential'] = 'ams'
      rconf['AMS'] = {
        'engine' => engine,
        'model' => model
      }
    else # ADF
      rconf['Potential']['potential'] = 'ams'
      rconf['AMS'] = {
        'engine' => engine,
        'xc' => model,
        'basis' => "DZ"
      }
    end
    rconf.filename = "#{pdir}/config.ini"
    rconf.write
    relative_diff(ldir, pdir)
    # Write path for eonclient
    bashfsc.write("\n cd \"#{pdir}/\"\n eonclient 2> error.txt 1> output.txt &")
    # End
    # Parameter scans
    if method == 'gprdimer'
      # Divisor T Dimer ==> 10 default
      [1, 10, 100, 1000].each do |divTdimer|
        dtddir = lower_dir("ddtd_#{divTdimer}", pdir) # LOWERING
        rconf['GPR Dimer']['divisor_t_dimer'] = divTdimer
        rconf.filename = "#{dtddir}/config.ini"
        rconf.write
        relative_diff(ldir, dtddir)
        # Write path for eonclient
        bashfsc.write("\n cd \"#{dtddir}/\"\n eonclient 2> error.txt 1> output.txt &")
        # End
        # Max step
        [0.1, 0.05, 0.2, 0.001].each do |maxstep|
          mxssdir = lower_dir("dmxss_#{maxstep}", dtddir) # LOWERING
          rconf['GPR Dimer']['max_step_size'] = maxstep
          rconf.filename = "#{mxssdir}/config.ini"
          rconf.write
          relative_diff(ldir, mxssdir)
          # Write path for eonclient
          bashfsc.write("\n cd \"#{mxssdir}/\"\n eonclient 2> error.txt 1> output.txt &")
          # End
          Dir.chdir(mxssdir) # max_step
        end
        Dir.chdir(dtddir) # divisor_t_dimer
      end
      # End parameter scan
    end
    Dir.chdir(pdir) # Potentials
  end
  Dir.chdir(ldir) # Methods
end

# Cleanup
['direction.dat', 'displacement.con', 'pos.con', 'full.sample.ini', 'rel.diff', 'config.ini'].each do |item|
  condel(item)
end
bashfsc.close

# # Need the following inputs
# reqInp = Set['config.ini', 'displacement.con', 'direction.dat', 'pos.con']
# allSys = ['ally_vinyl_ether_re', 'oxadiazole', 'sulfolene']

# puts "We are running from #{Dir.pwd}"
# puts "Set #{ENV["EONHOME"]}"
# FileUtils.mkdir_p(Dir.pwd + '/stuff')
# allSys.each do |i|
#    puts "Generating #{i}"
#    FileUtils.cp_r "#{Dir.pwd}/client/gtests/data/systems/#{i}/.", "#{Dir.pwd}/stuff"
# end
# puts reqInp == Dir.children("#{Dir.pwd}/stuff/").to_set
# puts reqInp ^ Dir.children("#{Dir.pwd}/stuff/")
# puts Dir.children("#{Dir.pwd}/stuff/").to_set
# # puts Dir["#{Dir.pwd}/stuff/*"]
