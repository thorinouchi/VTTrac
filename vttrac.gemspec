# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'numru/vttrac_version'

Gem::Specification.new do |spec|
  spec.name          = "vttrac"
  spec.version       = NumRu::VTTrac::VERSION
  spec.authors           = ["Takeshi Horinouchi"]

  #if spec.respond_to?(:metadata)
  #  spec.metadata['allowed_push_host'] = "TODO: Set to 'http://mygemserver.com' to prevent pushes to rubygems.org, or delete to allow pushes to any server."
  #end

  spec.summary          = %q{Velocimetry by Template Tracking}
  spec.description      = %q{This library provides the core C-language implementation and its Ruby interface. The algorithm used in this library is the simple template matching of PIV (particle image velocimetry) for monochromatic image-like data, but the matching is conducted multiple times in a Lagrangian manner as in PTV (particle tracking velocimetry) over a specified number of times.}
  spec.homepage         = 'https://github.com/thorinouchi/VTTrac'
  spec.licenses       = ["BSD-2-Clause"]

  #spec.files         = Dir['{ext/numru/*[cb],lib/numru/*.rb'].sort + %w(ChangeLog LICENSE README.md Rakefile vttrac.gemspec Gemfile)
  spec.files         = Dir['{ext/numru/*[cb],lib/numru/*.rb'].sort + %w(LICENSE README.md Rakefile vttrac.gemspec Gemfile)
  spec.require_paths = ["ext","lib"]
  spec.extensions << "ext/numru/extconf.rb"

  spec.post_install_message = "Thanks for installing!"

  spec.required_ruby_version = Gem::Requirement.new(">= 2.0")

  spec.add_runtime_dependency(%q<narray>, [">= 0.5.7"])
  spec.add_runtime_dependency(%q<narray_miss>, [">= 1.2.4"])
end
