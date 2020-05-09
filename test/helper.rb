# require 'simplecov'

# module SimpleCov::Configuration
#   def clean_filters
#     @filters = []
#   end
# end

# SimpleCov.configure do
#   clean_filters
#   load_adapter 'test_frameworks'
# end

# ENV["COVERAGE"] && SimpleCov.start do
#   add_filter "/.rvm/"
# end
# require 'rubygems'
# require 'bundler'
# begin
#   Bundler.setup(:default, :development)
# rescue Bundler::BundlerError => e
#   $stderr.puts e.message
#   $stderr.puts "Run `bundle install` to install missing gems"
#   exit e.status_code
# end
# require 'test/unit'
# require 'shoulda'

# $LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
# $LOAD_PATH.unshift(File.dirname(__FILE__))
# require 'bio-pangenome'

# class Test::Unit::TestCase
# end



# require 'coveralls'
# Coveralls.wear!

# require 'simplecov'
# SimpleCov.start do
#   add_filter %r{^/test/}
#   if ENV['CI']
#     formatter SimpleCov::Formatter::SimpleFormatter
#   else
#     formatter SimpleCov::Formatter::MultiFormatter.new([
#       SimpleCov::Formatter::SimpleFormatter,
#       SimpleCov::Formatter::HTMLFormatter
#     ])
#   end

#   track_files "**/*.rb"
# end

