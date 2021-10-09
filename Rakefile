# -* coding: utf-8 -*-
#require 'rake/testtask'
require 'rake/extensiontask'
require 'rake/packagetask'
begin
  require 'bundler/gem_helper'  # instead of 'bundler/gem_tasks' -> need manual
                                # calls of install_tasks (see below)
rescue LoadError
  puts 'If you want to create gem, You must install Bundler'
end

Bundler::GemHelper.install_tasks(name: "vttrac")

require './lib/numru/vttrac_version.rb'
def version
  NumRu::VTTrac::VERSION
end

#task :default => :test
#task :test => :compile
#Rake::TestTask.new do |t|
#  t.libs << 'lib' << 'test'
#  t.test_files = FileList['test/test_*.rb'].exclude('test/test_assoccoords.rb')
#end

task :default => :compile
Rake::ExtensionTask.new do |ext|
  ext.name = 'vttrac_core'
  ext.ext_dir = 'ext/numru'
  ext.lib_dir = 'lib/numru'
end

Rake::PackageTask.new('vttrac', "#{version}") do |t|
  t.need_tar_gz = true
  t.package_files.include `git ls-files`.split("\n")
end
