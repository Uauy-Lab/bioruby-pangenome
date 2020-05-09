SimpleCov.start do
   add_filter %r{^/test/}
   if ENV['CI']
     formatter SimpleCov::Formatter::SimpleFormatter
   else
     formatter SimpleCov::Formatter::MultiFormatter.new([
       SimpleCov::Formatter::SimpleFormatter,
       SimpleCov::Formatter::HTMLFormatter
     ])
   end
end

