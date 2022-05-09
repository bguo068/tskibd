project('tskibd', ['c', 'cpp'],
    version: files('VERSION.txt'),    
    default_options: ['c_std=c99', 'cpp_std=c++17']#, 'warning_level=3', 'werror=true']
)

# KASTORE
kastore_inc = include_directories('kastore/c')
kastore = static_library('kastore', 'kastore/c/kastore.c')
kastore_dep = declare_dependency(link_with : kastore, include_directories: kastore_inc)

# TSKIT
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)
tskit_inc = [kastore_inc, include_directories('tskit/c')]
# need tskit header as dep as we are compliling from outside
lib_deps = [m_dep, kastore_dep, declare_dependency(include_directories:tskit_inc)]
lib_sources = [
    'tskit/c/tskit/core.c', 'tskit/c/tskit/tables.c', 'tskit/c/tskit/trees.c',
    'tskit/c/tskit/genotypes.c', 'tskit/c/tskit/stats.c', 'tskit/c/tskit/convert.c', 'tskit/c/tskit/haplotype_matching.c']
tskit_lib = static_library('tskit',
    sources: lib_sources, dependencies: lib_deps)
tskit_dep = declare_dependency(include_directories:tskit_inc, link_with: tskit_lib)

# TSKIBD
executable('tskibd', 'src/tskibd.cpp', dependencies: [tskit_dep])
executable('newickibd', 'src/newickibd.cpp')