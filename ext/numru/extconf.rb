require 'mkmf'
gem_path = nil

if Gem.respond_to?(:find_files) and Gem.find_files("narray.h").length > 0
  require "rbconfig"
  so = RbConfig::CONFIG["DLEXT"]
  narray_include = File.expand_path(File.dirname(Gem.find_files("narray.h")[0]))
  narray_lib = File.expand_path(File.dirname(Gem.find_files("narray." + so)[0]))
else
  gem_home=(`gem environment GEM_HOME`).chomp
  narray_dir = Dir.glob("#{gem_home}/gems/narray-*/ext/narray").sort[-1]
  if narray_dir
    narray_include = narray_lib = narray_dir
  else
    narray_include = narray_lib = [ $sitearchdir, $vendorarchdir]
  end
end
dir_config('narray', narray_include, narray_lib)

require "narray"
if NArray.const_defined?(:SUPPORT_BIGMEM) && NArray::SUPPORT_BIGMEM

  case RbConfig::CONFIG["CC"]
  when "gcc"
    omp_opt = "-fopenmp"
  else
    omp_opt = nil
  end
  omp_opt = arg_config("--openmp", omp_opt)

  omp_opt = nil if omp_opt.to_s.empty?

  if omp_opt
    $CFLAGS << " " << omp_opt
    $DLDFLAGS << " " << omp_opt
    warn "OpenMP support: ON"
  else
    warn "OpenMP support: OFF"
    warn "if you want to enable openmp, set --openmp=compile_option"
  end
end

create_makefile('numru/vttrac_core')
