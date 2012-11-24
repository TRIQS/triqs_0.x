
# Should we use MathJax....
option(Sphinx_Math_Generator_MathJax " " OFF)
if (Sphinx_Math_Generator_MathJax)
 SET(SPHINX_MATH_GENERATOR "mathjax")
else (Sphinx_Math_Generator_MathJax)
 SET(SPHINX_MATH_GENERATOR "pngmath")
endif (Sphinx_Math_Generator_MathJax)

add_all_subdirectories_with_cmakelist()

# the local extensions....
file(COPY ../foreignlibs/breathe ../foreignlibs/doxylink ../foreignlibs/mpl ../foreignlibs/autocompile DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


