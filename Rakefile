# encoding: utf-8

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'juwelier'
Juwelier::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://guides.rubygems.org/specification-reference/ for more options
  gem.name = "bio-pangenome"
  gem.homepage = "http://github.com/Uauy-Lab/bioruby-pangenome"
  gem.license = "MIT"
  gem.summary = %Q{Scripts to analyse pangenomes.}
  gem.description = %Q{Tools to find similarity between pangenomes.}
  gem.email = "ricardo.ramirez-gonzalez@jic.ac.uk"
  gem.authors = ["Ricardo H. Ramirez-Gonzalez"]
  # dependencies defined in Gemfile
end

Juwelier::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  #test.libs << 'lib' << 'test'
  #test.pattern = 'test/**/test_*.rb'
  test.test_files = FileList['test/enable_coverage.rb', 'test/test*.rb']
  test.verbose = true
  test.warning = false
end

desc "Code coverage detail"
task :simplecov do
  ENV['COVERAGE'] = "true"
  Rake::Task['test'].execute
end

task :default => :test

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-pangenome #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
